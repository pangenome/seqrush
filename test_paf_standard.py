#!/usr/bin/env python3

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

# Test case: what should PAF look like for reverse alignment?
query = "ATCGATCG"
target = "CGATCGAT"  # This is RC of query

print("Test: Query aligns to reverse complement of target")
print(f"Query:      {query}")
print(f"Target:     {target}")
print(f"Target RC:  {reverse_complement(target)}")
print()

print("WFA2 aligns query to target_RC:")
print(f"Query:      {query}")
print(f"Target RC:  {reverse_complement(target)}")
print("This gives perfect match: 8M -> 8=")
print()

print("But in PAF format for reverse strand:")
print("- Query should be: query")  
print("- Target should be: target (forward strand)")
print("- Strand should be: -")
print("- CIGAR should represent: query vs target (not target_RC)")
print()

print("So the CIGAR should represent alignment of:")
print(f"Query:  {query}")
print(f"Target: {target}")
print("Which would be 8X (all mismatches), not 8=")
print()

print("OR, if PAF follows SAM convention:")
print("- The CIGAR represents alignment of RC(query) vs target")
print(f"RC(Query): {reverse_complement(query)}")
print(f"Target:    {target}")
print("Which would be 8= (perfect match)")

# Let's check what the actual expectation is
print("\nActual alignment that makes biological sense:")
print(f"Query:        {query}")
print(f"Target:       {target}")
print(f"RC(Query):    {reverse_complement(query)}")
print(f"RC(Target):   {reverse_complement(target)}")
print()
print("Query matches RC(Target), so in PAF:")
print("- We report the match between query and target")
print("- With strand = '-' to indicate target is reverse complemented")
print("- CIGAR should be 8= since query matches RC(target)")
print("- But coordinates should be relative to forward target")