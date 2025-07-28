#!/usr/bin/env python3

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
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences

def reverse_complement(seq):
    """Get reverse complement of sequence"""
    rc = []
    for b in reversed(seq):
        if b == 'A' or b == 'a': rc.append('T')
        elif b == 'T' or b == 't': rc.append('A')
        elif b == 'C' or b == 'c': rc.append('G')
        elif b == 'G' or b == 'g': rc.append('C')
        else: rc.append(b)
    return ''.join(rc)

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

def debug_alignment(query_seq, target_seq, cigar_str, is_reverse):
    """Debug the first few operations of an alignment"""
    print(f"Query length: {len(query_seq)}")
    print(f"Target length: {len(target_seq)}")
    print(f"Reverse: {is_reverse}")
    print(f"CIGAR: {cigar_str}")
    print()
    
    if is_reverse:
        target_seq = reverse_complement(target_seq)
        print(f"Target RC length: {len(target_seq)}")
    
    ops = parse_cigar(cigar_str)
    q_pos = 0
    t_pos = 0
    
    print("First 10 operations:")
    for i, (count, op) in enumerate(ops[:10]):
        print(f"  {i+1}. {count}{op}")
        
        if op == '=':
            print(f"    Should match at q_pos={q_pos}, t_pos={t_pos}")
            for j in range(min(count, 5)):  # Show first 5 of this operation
                if q_pos + j < len(query_seq) and t_pos + j < len(target_seq):
                    qc = query_seq[q_pos + j]
                    tc = target_seq[t_pos + j]
                    match = "✓" if qc == tc else "✗"
                    print(f"      [{q_pos + j},{t_pos + j}]: {qc} vs {tc} {match}")
                else:
                    print(f"      [{q_pos + j},{t_pos + j}]: OUT OF BOUNDS")
            q_pos += count
            t_pos += count
        elif op == 'X':
            print(f"    Should mismatch at q_pos={q_pos}, t_pos={t_pos}")
            q_pos += count
            t_pos += count
        elif op == 'I':
            print(f"    Insertion in target at t_pos={t_pos}")
            t_pos += count
        elif op == 'D':
            print(f"    Deletion from target at q_pos={q_pos}")
            q_pos += count
    
    print(f"\nAfter first 10 operations: q_pos={q_pos}, t_pos={t_pos}")
    
    # Calculate total consumption
    total_q = 0
    total_t = 0
    for count, op in ops:
        if op in ['=', 'X']:
            total_q += count
            total_t += count
        elif op == 'D':
            total_q += count
        elif op == 'I':
            total_t += count
    
    print(f"Total CIGAR consumption: query={total_q}, target={total_t}")

# Load sequences and test first PAF line
sequences = load_fasta('HLA-zoo/seqs/B-3106.fa')

# Read first PAF line that has errors (line 3 according to pafcheck)
with open('d.paf', 'r') as f:
    lines = f.readlines()
    
# Look at line 3 (pafcheck uses 1-based indexing)
line = lines[2].strip()  # 0-based indexing
fields = line.split('\t')

query_name = fields[0]
target_name = fields[5]
strand = fields[4]

# Extract CIGAR
cigar = None
for field in fields[12:]:
    if field.startswith('cg:Z:'):
        cigar = field[5:]
        break

print(f"Debugging PAF line 3:")
print(f"Query: {query_name}")
print(f"Target: {target_name}")
print(f"Strand: {strand}")
print()

if query_name in sequences and target_name in sequences and cigar:
    query_seq = sequences[query_name]
    target_seq = sequences[target_name]
    debug_alignment(query_seq, target_seq, cigar, strand == '-')
else:
    print("Sequences not found in FASTA!")