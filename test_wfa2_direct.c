#include <stdio.h>
#include <string.h>
#include "wavefront/wfa.h"

int main() {
    // Test sequences
    char* query = "ATCGATCGATCG";
    char* target = "ATCGATTTGATCG";
    
    // Configure alignment
    wavefront_aligner_attributes_t attributes = wavefront_aligner_attr_default;
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
    
    // Get CIGAR string
    char* buffer = malloc(1000);
    memset(buffer, 0, 1000);
    int pos = 0;
    
    for (int i = wf_aligner->cigar->begin_offset; i < wf_aligner->cigar->end_offset; i++) {
        buffer[pos++] = wf_aligner->cigar->operations[i];
    }
    
    fprintf(stdout, "CIGAR operations: %s\n", buffer);
    
    // Count operations
    int count_m = 0, count_x = 0, count_i = 0, count_d = 0;
    for (int i = 0; i < pos; i++) {
        switch(buffer[i]) {
            case 'M': count_m++; break;
            case 'X': count_x++; break;
            case 'I': count_i++; break;
            case 'D': count_d++; break;
        }
    }
    
    fprintf(stdout, "Operation counts: M=%d, X=%d, I=%d, D=%d\n", count_m, count_x, count_i, count_d);
    
    free(buffer);
    wavefront_aligner_delete(wf_aligner);
    
    return 0;
}