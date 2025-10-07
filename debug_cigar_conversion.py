#!/usr/bin/env python3

def trace_conversion(cigar, seq1, seq2):
    """Trace through CIGAR conversion to understand the issue"""
    print(f"Input CIGAR: {cigar}")
    print(f"Seq1: {seq1}")
    print(f"Seq2: {seq2}")
    print()
    
    pos1 = 0
    pos2 = 0
    result = []
    
    count = ""
    for ch in cigar:
        if ch.isdigit():
            count += ch
        else:
            n = int(count) if count else 1
            
            if ch == 'M':
                print(f"Processing {n}M:")
                for i in range(n):
                    if pos1 < len(seq1) and pos2 < len(seq2):
                        if seq1[pos1] == seq2[pos2]:
                            print(f"  pos1={pos1} ({seq1[pos1]}) == pos2={pos2} ({seq2[pos2]}) -> =")
                            result.append('=')
                        else:
                            print(f"  pos1={pos1} ({seq1[pos1]}) != pos2={pos2} ({seq2[pos2]}) -> X")
                            result.append('X')
                        pos1 += 1
                        pos2 += 1
            elif ch == 'I':
                print(f"Processing {n}I: skip {n} positions in target")
                for i in range(n):
                    result.append('I')
                    pos2 += 1
            elif ch == 'D':
                print(f"Processing {n}D: skip {n} positions in query")
                for i in range(n):
                    result.append('D')
                    pos1 += 1
                    
            count = ""
    
    print()
    print(f"Result: {''.join(result)}")
    print(f"Final positions: pos1={pos1}/{len(seq1)}, pos2={pos2}/{len(seq2)}")

# Test case
trace_conversion("6MI1X5M", "ATCGATCGATCG", "ATCGATTTGATCG")