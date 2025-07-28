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

# Load sequences
sequences = load_fasta('HLA-zoo/seqs/B-3106.fa')

# Check the problematic line (line 3)
with open('d.paf', 'r') as f:
    lines = f.readlines()

line = lines[2].strip()  # Line 3 (0-based)
fields = line.split('\t')

query_name = fields[0]
query_len = int(fields[1])
query_start = int(fields[2])
query_end = int(fields[3])
strand = fields[4]
target_name = fields[5]
target_len = int(fields[6])
target_start = int(fields[7])
target_end = int(fields[8])

print(f"PAF line 3 coordinates:")
print(f"Query: {query_name} len={query_len} start={query_start} end={query_end}")
print(f"Target: {target_name} len={target_len} start={target_start} end={target_end}")
print(f"Strand: {strand}")
print()

# Check actual sequence lengths
if query_name in sequences and target_name in sequences:
    actual_query_len = len(sequences[query_name])
    actual_target_len = len(sequences[target_name])
    
    print(f"Actual sequence lengths:")
    print(f"Query: {actual_query_len}")
    print(f"Target: {actual_target_len}")
    print()
    
    if actual_query_len != query_len:
        print(f"ERROR: Query length mismatch! PAF={query_len}, Actual={actual_query_len}")
    if actual_target_len != target_len:
        print(f"ERROR: Target length mismatch! PAF={target_len}, Actual={actual_target_len}")
    
    # Check if we're reporting global alignment incorrectly
    if query_start != 0 or query_end != actual_query_len:
        print(f"WARNING: Not reporting full query range: {query_start}-{query_end} vs 0-{actual_query_len}")
    if target_start != 0 or target_end != actual_target_len:
        print(f"WARNING: Not reporting full target range: {target_start}-{target_end} vs 0-{actual_target_len}")

print("\nChecking some other lines for pattern...")

# Check a few more lines
for i in [0, 1, 4, 7]:  # Lines 1, 2, 5, 8
    if i < len(lines):
        line = lines[i].strip()
        fields = line.split('\t')
        
        query_name = fields[0]
        query_len = int(fields[1])
        target_name = fields[5]
        target_len = int(fields[6])
        
        if query_name in sequences and target_name in sequences:
            actual_query_len = len(sequences[query_name])
            actual_target_len = len(sequences[target_name])
            
            print(f"Line {i+1}: Query PAF={query_len} Actual={actual_query_len}, Target PAF={target_len} Actual={actual_target_len}")
            
            if actual_query_len != query_len or actual_target_len != target_len:
                print(f"  MISMATCH!")