#!/usr/bin/env python3
import sys

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

def apply_cigar(query_seq, target_seq, cigar_str):
    """Apply CIGAR to sequences and show alignment"""
    ops = parse_cigar(cigar_str)
    
    query_aligned = ""
    target_aligned = ""
    match_line = ""
    
    q_pos = 0
    t_pos = 0
    
    for count, op in ops:
        if op == '=' or op == 'M':
            # Match
            for i in range(count):
                if q_pos < len(query_seq) and t_pos < len(target_seq):
                    query_aligned += query_seq[q_pos]
                    target_aligned += target_seq[t_pos]
                    if query_seq[q_pos] == target_seq[t_pos]:
                        match_line += "|"
                    else:
                        match_line += "X"
                    q_pos += 1
                    t_pos += 1
                else:
                    print(f"ERROR: Out of bounds at {op} operation")
                    return
        elif op == 'X':
            # Mismatch
            for i in range(count):
                if q_pos < len(query_seq) and t_pos < len(target_seq):
                    query_aligned += query_seq[q_pos]
                    target_aligned += target_seq[t_pos]
                    match_line += "X"
                    q_pos += 1
                    t_pos += 1
                else:
                    print(f"ERROR: Out of bounds at {op} operation")
                    return
        elif op == 'I':
            # Insertion in target
            for i in range(count):
                if t_pos < len(target_seq):
                    query_aligned += "-"
                    target_aligned += target_seq[t_pos]
                    match_line += " "
                    t_pos += 1
                else:
                    print(f"ERROR: Out of bounds at {op} operation")
                    return
        elif op == 'D':
            # Deletion from target (present in query)
            for i in range(count):
                if q_pos < len(query_seq):
                    query_aligned += query_seq[q_pos]
                    target_aligned += "-"
                    match_line += " "
                    q_pos += 1
                else:
                    print(f"ERROR: Out of bounds at {op} operation")
                    return
    
    print(f"Query consumed: {q_pos}/{len(query_seq)}")
    print(f"Target consumed: {t_pos}/{len(target_seq)}")
    print()
    
    # Print alignment in chunks
    chunk_size = 80
    for i in range(0, len(query_aligned), chunk_size):
        print(f"Query:  {query_aligned[i:i+chunk_size]}")
        print(f"        {match_line[i:i+chunk_size]}")
        print(f"Target: {target_aligned[i:i+chunk_size]}")
        print()

if __name__ == "__main__":
    # Test with our simple example
    query = "ATCGATCGATCG"
    target = "ATCGATTTGATCG"
    cigar = "6=1I1X5="
    
    print("Test case:")
    print(f"Query:  {query} (len={len(query)})")
    print(f"Target: {target} (len={len(target)})")
    print(f"CIGAR:  {cigar}")
    print()
    
    apply_cigar(query, target, cigar)