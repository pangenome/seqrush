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
        if op == '=' or op == 'M':
            for i in range(count):
                if q_pos >= len(query_seq) or t_pos >= len(target_seq):
                    errors.append(f"Out of bounds at {op}{count}: q_pos={q_pos}, t_pos={t_pos}")
                    return errors
                if query_seq[q_pos] != target_seq[t_pos]:
                    errors.append(f"Mismatch at {op} operation: query[{q_pos}]={query_seq[q_pos]} != target[{t_pos}]={target_seq[t_pos]}")
                q_pos += 1
                t_pos += 1
        elif op == 'X':
            # Mismatch - advance both
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

# Test with first PAF line
if len(sys.argv) < 3:
    print("Usage: validate_paf_simple.py <fasta> <paf>")
    sys.exit(1)

sequences = load_fasta(sys.argv[1])
print(f"Loaded {len(sequences)} sequences")

# Read first PAF line
with open(sys.argv[2], 'r') as f:
    line = f.readline().strip()
    fields = line.split('\t')
    
    query_name = fields[0]
    query_len = int(fields[1])
    target_name = fields[5]
    target_len = int(fields[6])
    strand = fields[4]
    
    # Extract CIGAR from cg:Z: tag
    cigar = None
    for field in fields[12:]:
        if field.startswith('cg:Z:'):
            cigar = field[5:]
            break
    
    print(f"\nValidating alignment:")
    print(f"Query: {query_name} (len={query_len})")
    print(f"Target: {target_name} (len={target_len})")
    print(f"Strand: {strand}")
    print(f"CIGAR: {cigar}")
    
    if query_name in sequences and target_name in sequences:
        query_seq = sequences[query_name]
        target_seq = sequences[target_name]
        
        errors = validate_alignment(query_seq, target_seq, cigar, strand == '-')
        
        if errors:
            print(f"\nErrors found ({len(errors)}):")
            for i, err in enumerate(errors[:10]):  # Show first 10 errors
                print(f"  {i+1}. {err}")
            if len(errors) > 10:
                print(f"  ... and {len(errors) - 10} more errors")
        else:
            print("\nAlignment validated successfully!")
    else:
        print(f"\nError: Sequences not found in FASTA")
        if query_name not in sequences:
            print(f"  Query '{query_name}' not found")
        if target_name not in sequences:
            print(f"  Target '{target_name}' not found")