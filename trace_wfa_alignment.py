#!/usr/bin/env python3

def trace_alignment(cigar, seq1, seq2):
    """Trace WFA2 alignment"""
    print(f"CIGAR: {cigar}")
    print(f"Seq1:  {seq1}")
    print(f"Seq2:  {seq2}")
    print()
    
    aligned1 = []
    aligned2 = []
    match_line = []
    
    pos1 = 0
    pos2 = 0
    
    count = ""
    for ch in cigar:
        if ch.isdigit():
            count += ch
        else:
            n = int(count) if count else 1
            
            if ch == 'M':
                for _ in range(n):
                    aligned1.append(seq1[pos1])
                    aligned2.append(seq2[pos2])
                    if seq1[pos1] == seq2[pos2]:
                        match_line.append('|')
                    else:
                        match_line.append('X')
                    pos1 += 1
                    pos2 += 1
            elif ch == 'I':
                for _ in range(n):
                    aligned1.append('-')
                    aligned2.append(seq2[pos2])
                    match_line.append(' ')
                    pos2 += 1
            elif ch == 'D':
                for _ in range(n):
                    aligned1.append(seq1[pos1])
                    aligned2.append('-')
                    match_line.append(' ')
                    pos1 += 1
            elif ch == 'X':
                for _ in range(n):
                    aligned1.append(seq1[pos1])
                    aligned2.append(seq2[pos2])
                    match_line.append('X')
                    pos1 += 1
                    pos2 += 1
                    
            count = ""
    
    print("Alignment:")
    print(f"Seq1: {''.join(aligned1)}")
    print(f"      {''.join(match_line)}")
    print(f"Seq2: {''.join(aligned2)}")
    print()
    print(f"Consumed: seq1={pos1}/{len(seq1)}, seq2={pos2}/{len(seq2)}")

# Trace the actual WFA2 CIGAR
trace_alignment("6MI1X5M", "ATCGATCGATCG", "ATCGATTTGATCG")