// Instrumentation patch for seqwish to debug RC alignments
// This adds debug logging to transclosure.cpp

#include <iostream>
#include <sstream>

// Debug function to log position uniting operations
void debug_unite_positions(uint64_t j, uint64_t offset_p, 
                          const pos_t& p, const range_t& r,
                          const seqindex_t& seqidx) {
    // Only log for our test sequences
    std::string seq1_name = seqidx.seq_name(r.query_pos);
    std::string seq2_name = seqidx.seq_name(p);
    
    if ((seq1_name.find("seq1") != std::string::npos || 
         seq1_name.find("seq2") != std::string::npos) &&
        (seq2_name.find("seq1") != std::string::npos || 
         seq2_name.find("seq2") != std::string::npos)) {
        
        std::cerr << "[SEQWISH_DEBUG] Unite: "
                  << seq1_name << " pos " << j - r.query_pos.offset 
                  << " (global " << j << ")"
                  << " <-> "
                  << seq2_name << " pos " << offset_p - offset(p) 
                  << " (global " << offset_p << ")"
                  << " is_rev=" << is_rev(p)
                  << std::endl;
                  
        if (j < r.query_pos.offset + 5) {  // Log first few positions
            std::cerr << "  Query rank: " << q_curr_rank(j) 
                      << " Target rank: " << q_curr_rank(offset_p) 
                      << std::endl;
        }
    }
}