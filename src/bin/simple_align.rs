use std::fs::File;
use std::io::{BufReader, BufRead};
use lib_wfa2::affine_wavefront::{AffineWavefronts, MemoryMode, AlignmentStatus};

#[derive(Debug, Clone)]
struct Sequence {
    id: String,
    data: Vec<u8>,
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b,
        })
        .collect()
}

fn load_fasta(path: &str) -> Result<Vec<Sequence>, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut sequences = Vec::new();
    let mut current_id = String::new();
    let mut current_data = Vec::new();
    
    for line in reader.lines() {
        let line = line?;
        if let Some(stripped) = line.strip_prefix('>') {
            if !current_id.is_empty() {
                sequences.push(Sequence {
                    id: current_id.clone(),
                    data: current_data.clone(),
                });
                current_data.clear();
            }
            current_id = stripped.split_whitespace().next().unwrap_or("").to_string();
        } else {
            current_data.extend(line.trim().bytes());
        }
    }
    
    if !current_id.is_empty() {
        sequences.push(Sequence {
            id: current_id,
            data: current_data,
        });
    }
    
    Ok(sequences)
}

fn invert_cigar(cigar: &str) -> String {
    let mut result = String::new();
    let mut count = 0;
    
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            count = count * 10 + (ch as usize - '0' as usize);
        } else {
            if count == 0 { count = 1; }
            
            let op = match ch {
                'I' => 'D',
                'D' => 'I',
                _ => ch,
            };
            
            result.push_str(&count.to_string());
            result.push(op);
            count = 0;
        }
    }
    
    result
}

fn convert_to_eqx(cigar: &str, query: &[u8], target: &[u8]) -> String {
    let mut result = String::new();
    let mut count = 0;
    let mut qpos = 0;
    let mut tpos = 0;
    
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            count = count * 10 + (ch as usize - '0' as usize);
        } else {
            if count == 0 { count = 1; }
            
            match ch {
                'M' => {
                    // Convert M to =/X
                    for _ in 0..count {
                        if qpos < query.len() && tpos < target.len() && query[qpos] == target[tpos] {
                            result.push_str("1=");
                        } else {
                            result.push_str("1X");
                        }
                        qpos += 1;
                        tpos += 1;
                    }
                }
                'I' => {
                    result.push_str(&count.to_string());
                    result.push('I');
                    qpos += count;
                }
                'D' => {
                    result.push_str(&count.to_string());
                    result.push('D');
                    tpos += count;
                }
                '=' | 'X' => {
                    result.push_str(&count.to_string());
                    result.push(ch);
                    if ch == '=' || ch == 'X' {
                        qpos += count;
                        tpos += count;
                    }
                }
                _ => {}
            }
            count = 0;
        }
    }
    
    // Compact the CIGAR
    compact_cigar(&result)
}

fn compact_cigar(cigar: &str) -> String {
    let mut result = String::new();
    let mut count = 0;
    let mut current_op: Option<char> = None;
    let mut accumulated = 0;
    
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            count = count * 10 + (ch as usize - '0' as usize);
        } else {
            if count == 0 { count = 1; }
            
            if Some(ch) == current_op {
                accumulated += count;
            } else {
                if let Some(op) = current_op {
                    result.push_str(&accumulated.to_string());
                    result.push(op);
                }
                current_op = Some(ch);
                accumulated = count;
            }
            count = 0;
        }
    }
    
    if let Some(op) = current_op {
        result.push_str(&accumulated.to_string());
        result.push(op);
    }
    
    result
}

fn parse_cigar_stats(cigar: &str) -> (usize, usize, usize) {
    let mut matches = 0;
    let mut qlen = 0;
    let mut tlen = 0;
    let mut count = 0;
    
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            count = count * 10 + (ch as usize - '0' as usize);
        } else {
            if count == 0 { count = 1; }
            
            match ch {
                '=' => {
                    matches += count;
                    qlen += count;
                    tlen += count;
                }
                'X' | 'M' => {
                    qlen += count;
                    tlen += count;
                }
                'I' => {
                    qlen += count;
                }
                'D' => {
                    tlen += count;
                }
                _ => {}
            }
            count = 0;
        }
    }
    
    (matches, qlen, tlen)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 {
        eprintln!("Usage: {} <fasta_file>", args[0]);
        std::process::exit(1);
    }
    
    let sequences = load_fasta(&args[1])?;
    eprintln!("Loaded {} sequences", sequences.len());
    
    // Alignment parameters matching wfmash defaults
    let match_score = 0;
    let mismatch_penalty = 5;
    let gap_open1 = 8;
    let gap_extend1 = 2;
    let gap_open2 = 24;
    let gap_extend2 = 1;
    
    // All vs all alignment
    for i in 0..sequences.len() {
        for j in 0..sequences.len() {
            let query = &sequences[i];
            let target = &sequences[j];
            
            // Try both orientations
            for (is_reverse, target_seq) in [(false, target.data.clone()), (true, reverse_complement(&target.data))] {
                // Create WFA2 aligner
                let mut wf = AffineWavefronts::with_penalties_affine2p_and_memory_mode(
                    match_score,
                    mismatch_penalty,
                    gap_open1,
                    gap_extend1,
                    gap_open2,
                    gap_extend2,
                    MemoryMode::Ultralow
                );
                
                // Set alignment type - maybe this matters?
                wf.set_alignment_span(lib_wfa2::affine_wavefront::AlignmentSpan::End2End);
                
                // Align query to target
                // WFA2 expects (pattern, text) = (query, target)
                let status = wf.align(&query.data, &target_seq);
                
                if matches!(status, AlignmentStatus::Completed) {
                    let score = wf.score().abs();
                    let cigar_bytes = wf.cigar();
                    let wfa_cigar = String::from_utf8_lossy(cigar_bytes);
                    
                    // Debug specific alignment
                    if query.id.contains("157734152") && target.id.contains("568815592") && !is_reverse {
                        eprintln!("DEBUG: Aligning {} vs {} (forward strand)", query.id, target.id);
                        eprintln!("  WFA2 CIGAR (first 200): {}", &wfa_cigar.chars().take(200).collect::<String>());
                        eprintln!("  Score: {}", score);
                        eprintln!("  Matches from CIGAR after EQX conversion: (will compute below)");
                        
                        // Check first 100bp similarity
                        let mut matches_100 = 0;
                        for i in 0..100.min(query.data.len()).min(target.data.len()) {
                            if query.data[i] == target.data[i] {
                                matches_100 += 1;
                            }
                        }
                        eprintln!("  First 100bp similarity: {}/100", matches_100);
                        
                        // Let's also try the other way
                        let mut wf2 = AffineWavefronts::with_penalties_affine2p_and_memory_mode(
                            match_score,
                            mismatch_penalty,
                            gap_open1,
                            gap_extend1,
                            gap_open2,
                            gap_extend2,
                            MemoryMode::Ultralow
                        );
                        wf2.set_alignment_span(lib_wfa2::affine_wavefront::AlignmentSpan::End2End);
                        let status2 = wf2.align(&target.data, &query.data);
                        if matches!(status2, AlignmentStatus::Completed) {
                            let cigar2 = wf2.cigar();
                            let wfa_cigar2 = String::from_utf8_lossy(cigar2);
                            let inv_cigar2 = invert_cigar(&wfa_cigar2);
                            let eqx_cigar2 = convert_to_eqx(&inv_cigar2, &query.data, &target.data);
                            let (matches2, _, _) = parse_cigar_stats(&eqx_cigar2);
                            eprintln!("  REVERSE TEST: target->query alignment gives {} matches", matches2);
                        }
                    }
                    
                    // Debug: let's see what raw WFA2 gives us
                    if query.id.contains("568815592") && target.id.contains("157734152") && !is_reverse {
                        eprintln!("\nRAW WFA2 output:");
                        eprintln!("  Query: {} (len {})", query.id, query.data.len());
                        eprintln!("  Target: {} (len {})", target.id, target.data.len());
                        eprintln!("  Raw CIGAR (first 100): {}", &wfa_cigar.chars().take(100).collect::<String>());
                    }
                    
                    // WFA2 already gives us query->target CIGAR since we passed (query, target)
                    let paf_cigar = wfa_cigar.to_string();
                    
                    // Convert to EQX format
                    let eqx_cigar = if is_reverse {
                        let query_rc = reverse_complement(&query.data);
                        convert_to_eqx(&paf_cigar, &query_rc, &target.data)
                    } else {
                        convert_to_eqx(&paf_cigar, &query.data, &target.data)
                    };
                    
                    // Get stats
                    let (matches, qlen, tlen) = parse_cigar_stats(&eqx_cigar);
                    let block_len = qlen.max(tlen);
                    
                    // Debug: show exact comparison for specific alignment
                    if query.id.contains("568815592") && target.id.contains("157734152") && !is_reverse {
                        eprintln!("\nDEBUG CIGAR comparison:");
                        eprintln!("  Our CIGAR (first 100): {}", &eqx_cigar.chars().take(100).collect::<String>());
                        eprintln!("  Expected from wfmash: 159=1X75=1X32=2X1=1X41=1X85=1X73=1X52=2X32=1I6=11I3=1X...");
                        eprintln!("  Matches: {} (expected 3252)", matches);
                    }
                    
                    // Output PAF
                    println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tNM:i:{}\tAS:i:{}\tcg:Z:{}",
                        query.id,
                        query.data.len(),
                        0,  // query start
                        query.data.len(),  // query end
                        if is_reverse { '-' } else { '+' },
                        target.id,
                        target.data.len(),
                        0,  // target start
                        target.data.len(),  // target end
                        matches,
                        block_len,
                        60,  // mapping quality
                        score,  // edit distance
                        -score,  // alignment score
                        eqx_cigar
                    );
                }
            }
        }
    }
    
    Ok(())
}