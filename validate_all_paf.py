#!/usr/bin/env python3
import sys

def load_fasta(filename):
    """Load FASTA sequences into a dictionary"""
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                # Extract ID as first word after >
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences

def parse_cigar(cigar):
    """Parse CIGAR string into operations"""
    ops = []
    count = ""
    for c in cigar:
        if c.isdigit():
            count += c
        else:
            ops.append((int(count) if count else 1, c))
            count = ""
    return ops

def validate_alignment(query_seq, target_seq, cigar_str, is_reverse):
    """Validate CIGAR against sequences"""
    if is_reverse:
        # Reverse complement target
        rc = []
        for b in reversed(target_seq):
            if b == 'A' or b == 'a': rc.append('T')
            elif b == 'T' or b == 't': rc.append('A')
            elif b == 'C' or b == 'c': rc.append('G')
            elif b == 'G' or b == 'g': rc.append('C')
            else: rc.append(b)
        target_seq = ''.join(rc)
    
    ops = parse_cigar(cigar_str)
    q_pos = 0
    t_pos = 0
    errors = []
    
    for count, op in ops:
        if op == '=':
            # Matches
            for i in range(count):
                if q_pos >= len(query_seq) or t_pos >= len(target_seq):
                    errors.append(f"Out of bounds at {op}{count}: q_pos={q_pos}, t_pos={t_pos}")
                    return errors
                if query_seq[q_pos] != target_seq[t_pos]:
                    errors.append(f"Match claimed but mismatch at: query[{q_pos}]={query_seq[q_pos]} != target[{t_pos}]={target_seq[t_pos]}")
                q_pos += 1
                t_pos += 1
        elif op == 'X':
            # Mismatch - just advance
            q_pos += count
            t_pos += count
        elif op == 'I':
            # Insertion in target
            t_pos += count
        elif op == 'D':
            # Deletion from target
            q_pos += count
    
    if q_pos != len(query_seq):
        errors.append(f"Query not fully consumed: {q_pos}/{len(query_seq)}")
    if t_pos != len(target_seq):
        errors.append(f"Target not fully consumed: {t_pos}/{len(target_seq)}")
    
    return errors

# Main validation
if len(sys.argv) < 3:
    print("Usage: validate_all_paf.py <fasta> <paf>")
    sys.exit(1)

sequences = load_fasta(sys.argv[1])
print(f"Loaded {len(sequences)} sequences")

total_alignments = 0
valid_alignments = 0
error_alignments = 0

with open(sys.argv[2], 'r') as f:
    for line_num, line in enumerate(f, 1):
        if not line.strip():
            continue
        
        fields = line.strip().split('\t')
        if len(fields) < 12:
            continue
        
        total_alignments += 1
        
        query_name = fields[0]
        target_name = fields[5]
        strand = fields[4]
        
        # Extract CIGAR
        cigar = None
        for field in fields[12:]:
            if field.startswith('cg:Z:'):
                cigar = field[5:]
                break
        
        if query_name in sequences and target_name in sequences and cigar:
            query_seq = sequences[query_name]
            target_seq = sequences[target_name]
            
            errors = validate_alignment(query_seq, target_seq, cigar, strand == '-')
            
            if errors:
                error_alignments += 1
                if error_alignments <= 5:  # Show first 5 problematic alignments
                    print(f"\nLine {line_num}: {query_name} vs {target_name} ({strand})")
                    print(f"  CIGAR: {cigar}")
                    print(f"  Errors: {errors[0]}")
            else:
                valid_alignments += 1

print(f"\nSummary:")
print(f"Total alignments: {total_alignments}")
print(f"Valid alignments: {valid_alignments}")
print(f"Error alignments: {error_alignments}")
print(f"Success rate: {valid_alignments/total_alignments*100:.1f}%")