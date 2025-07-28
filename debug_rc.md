# Debug RC Alignment Logic

## Example:
- Seq1: ACGT (offset 0)
- Seq2: ACGT (offset 10)
- Alignment: seq1 forward vs seq2 RC

## What allwave does:
- Compares seq1 "ACGT" with RC(seq2) = "ACGT" -> "TGCA"
- Finds no match (ACGT != TGCA)

## What should happen:
- If seq1 aligns to seq2 in RC orientation:
  - seq1[0]=A matches complement of seq2[3]=T (which is A)
  - seq1[1]=C matches complement of seq2[2]=G (which is C)
  - etc.

## Current unite_matching_region logic:
For RC alignment (seq2_is_rc = true):
1. For each position i in the match:
   - pos1_fwd = seq1[i] forward
   - rc_local_pos = i (position in RC sequence)
   - forward_local_pos = seq2_len - 1 - i (corresponding position in original seq2)
   - pos2_rev = seq2[forward_local_pos] reverse
   - Unite: pos1_fwd with pos2_rev

This seems correct! seq1[i] forward is united with seq2[len-1-i] reverse.

## The real issue:
When we build the graph, we only consider forward positions. But if a position was only ever aligned in reverse orientation, it won't be connected to the rest of the graph.

## Solution:
We need to ensure that forward and reverse orientations of the same position are properly connected when they represent the same sequence content.