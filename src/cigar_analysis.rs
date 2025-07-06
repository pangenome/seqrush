/// CIGAR analysis for detecting gaps that might indicate inversions

#[derive(Debug, Clone)]
pub struct Gap {
    pub query_start: usize,
    pub query_end: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub gap_type: GapType,
}

#[derive(Debug, Clone, PartialEq)]
pub enum GapType {
    /// Both query and target have unaligned regions
    Divergent,
    /// Only query has unaligned region
    QueryOnly,
    /// Only target has unaligned region
    TargetOnly,
}

/// Parse CIGAR string to find large gaps that might be inversions
pub fn find_potential_inversion_sites(
    cigar: &str,
    min_gap_size: usize,
) -> Vec<Gap> {
    let mut gaps = Vec::new();
    
    // Parse CIGAR into operations
    let mut operations = Vec::new();
    let mut count = 0;
    
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            count = count * 10 + (ch as usize - '0' as usize);
        } else {
            if count == 0 {
                count = 1;
            }
            operations.push((ch, count));
            count = 0;
        }
    }
    
    // Track positions
    let mut query_pos = 0;
    let mut target_pos = 0;
    
    // Look for patterns like: matches -> large indels -> matches
    let mut i = 0;
    while i < operations.len() {
        let (op, count) = operations[i];
        
        match op {
            'M' | '=' => {
                // Look ahead for large gaps
                if i + 1 < operations.len() {
                    let mut j = i + 1;
                    let mut query_gap = 0;
                    let mut target_gap = 0;
                    let gap_start_query = query_pos + count;
                    let gap_start_target = target_pos + count;
                    
                    // Accumulate gap size
                    while j < operations.len() {
                        let (next_op, next_count) = operations[j];
                        match next_op {
                            'I' => target_gap += next_count,
                            'D' => query_gap += next_count,
                            'X' => {
                                query_gap += next_count;
                                target_gap += next_count;
                            }
                            'M' | '=' => {
                                // End of gap region
                                break;
                            }
                            _ => {}
                        }
                        j += 1;
                    }
                    
                    // Check if this is a significant gap
                    if query_gap >= min_gap_size && target_gap >= min_gap_size {
                        gaps.push(Gap {
                            query_start: gap_start_query,
                            query_end: gap_start_query + query_gap,
                            target_start: gap_start_target,
                            target_end: gap_start_target + target_gap,
                            gap_type: GapType::Divergent,
                        });
                    } else if query_gap >= min_gap_size {
                        gaps.push(Gap {
                            query_start: gap_start_query,
                            query_end: gap_start_query + query_gap,
                            target_start: gap_start_target,
                            target_end: gap_start_target,
                            gap_type: GapType::QueryOnly,
                        });
                    } else if target_gap >= min_gap_size {
                        gaps.push(Gap {
                            query_start: gap_start_query,
                            query_end: gap_start_query,
                            target_start: gap_start_target,
                            target_end: gap_start_target + target_gap,
                            gap_type: GapType::TargetOnly,
                        });
                    }
                }
                
                query_pos += count;
                target_pos += count;
            }
            'X' => {
                query_pos += count;
                target_pos += count;
            }
            'I' => {
                target_pos += count;
            }
            'D' => {
                query_pos += count;
            }
            _ => {}
        }
        
        i += 1;
    }
    
    gaps
}

/// Check if a gap is large enough to potentially contain an inversion
pub fn is_potential_inversion(gap: &Gap, min_inversion_size: usize) -> bool {
    match gap.gap_type {
        GapType::Divergent => {
            let query_size = gap.query_end - gap.query_start;
            let target_size = gap.target_end - gap.target_start;
            
            // For inversions, the sizes should be similar
            let size_ratio = query_size.max(target_size) as f64 / query_size.min(target_size) as f64;
            
            query_size >= min_inversion_size && 
            target_size >= min_inversion_size && 
            size_ratio <= 1.5  // Allow up to 50% size difference
        }
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_simple_gap_detection() {
        // Simple case: 10M 5D 5I 10M
        let gaps = find_potential_inversion_sites("10M5D5I10M", 3);
        
        assert_eq!(gaps.len(), 1);
        assert_eq!(gaps[0].query_start, 10);
        assert_eq!(gaps[0].query_end, 15);
        assert_eq!(gaps[0].target_start, 10);
        assert_eq!(gaps[0].target_end, 15);
    }
    
    #[test]
    fn test_large_divergent_gap() {
        // Large divergent region that could be an inversion
        let gaps = find_potential_inversion_sites("20M50X20M", 30);
        
        assert_eq!(gaps.len(), 1);
        assert_eq!(gaps[0].gap_type, GapType::Divergent);
        assert_eq!(gaps[0].query_end - gaps[0].query_start, 50);
        assert_eq!(gaps[0].target_end - gaps[0].target_start, 50);
    }
    
    #[test]
    fn test_no_gaps() {
        let gaps = find_potential_inversion_sites("100M", 10);
        assert_eq!(gaps.len(), 0);
    }
    
    #[test]
    fn test_multiple_gaps() {
        let gaps = find_potential_inversion_sites("10M20D20I10M30X10M", 15);
        
        // Should find two gaps
        assert_eq!(gaps.len(), 2);
    }
    
    #[test]
    fn test_is_potential_inversion() {
        let gap = Gap {
            query_start: 100,
            query_end: 200,
            target_start: 100,
            target_end: 195,
            gap_type: GapType::Divergent,
        };
        
        assert!(is_potential_inversion(&gap, 50));
        assert!(!is_potential_inversion(&gap, 150)); // Too large minimum
        
        let gap2 = Gap {
            query_start: 100,
            query_end: 200,
            target_start: 100,
            target_end: 100,
            gap_type: GapType::QueryOnly,
        };
        
        assert!(!is_potential_inversion(&gap2, 50)); // Not divergent
    }
}