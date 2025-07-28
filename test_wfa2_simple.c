#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "wavefront/wfa.h"

int main() {
    // Test sequences
    char* query = "ATCGATCGATCG";
    char* target = "ATCGATTTGATCG";
    
    // Configure alignment
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.distance_metric = gap_affine;
    attributes.affine_penalties.mismatch = 4;
    attributes.affine_penalties.gap_opening = 6;
    attributes.affine_penalties.gap_extension = 2;
    
    // Initialize aligner
    wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
    
    // Align
    wavefront_align(wf_aligner, query, strlen(query), target, strlen(target));
    
    // Extract CIGAR
    fprintf(stdout, "Query:  %s (len=%zu)\n", query, strlen(query));
    fprintf(stdout, "Target: %s (len=%zu)\n", target, strlen(target));
    fprintf(stdout, "Score: %d\n", wf_aligner->cigar->score);
    
    // Print compact CIGAR
    cigar_print_pretty(stdout, wf_aligner->cigar, query, strlen(query), target, strlen(target));
    
    wavefront_aligner_delete(wf_aligner);
    
    return 0;
}