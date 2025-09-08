use crate::pos::{make_pos, Pos};
use std::sync::Arc;
use uf_rush::UFRush;

/// Bidirected union-find that properly handles orientation-aware positions
/// This wrapper ensures that positions with different orientations are treated correctly
pub struct BidirectedUnionFind {
    /// Underlying union-find structure
    /// We use raw u64 values from Pos as the element IDs
    uf: Arc<UFRush>,
}

impl BidirectedUnionFind {
    /// Create a new bidirected union-find for positions up to max_offset
    /// We need to handle the encoded position values which include orientation bit
    pub fn new(max_offset: usize) -> Self {
        // Positions are encoded as (offset << 1) | orientation
        // So the maximum value is approximately (max_offset << 1) + 1
        // We need extra capacity to handle all possible encoded values
        let capacity = (max_offset << 1) + 2;
        Self {
            uf: Arc::new(UFRush::new(capacity)),
        }
    }

    /// Find the representative position for a given position
    pub fn find(&self, pos: Pos) -> Pos {
        let raw_pos = pos;
        let representative = self.uf.find(raw_pos as usize);
        representative as u64
    }

    /// Unite two positions in the union-find
    /// This preserves the orientation information in the union
    pub fn unite(&self, pos1: Pos, pos2: Pos) {
        let raw1 = pos1;
        let raw2 = pos2;

        // Only unite if they're different positions
        if raw1 != raw2 {
            self.uf.unite(raw1 as usize, raw2 as usize);
        }
    }

    /// Check if two positions are in the same connected component
    pub fn same(&self, pos1: Pos, pos2: Pos) -> bool {
        if pos1 == pos2 {
            return true;
        }

        let rep1 = self.uf.find(pos1 as usize);
        let rep2 = self.uf.find(pos2 as usize);
        rep1 == rep2
    }

    /// Unite positions from a matching region, handling reverse complement properly
    /// seq1_is_rc indicates if seq1 (query) was reverse complemented for alignment
    /// seq1_len is the length of seq1 (used for RC coordinate transformation)
    #[allow(clippy::too_many_arguments)]
    pub fn unite_matching_region(
        &self,
        seq1_offset: usize,
        seq2_offset: usize,
        seq1_local_start: usize,
        seq2_local_start: usize,
        match_length: usize,
        seq1_is_rc: bool,
        seq1_len: usize,
    ) {
        // Debug for specific case
        let debug = seq1_is_rc && (match_length > 100 || match_length == 12 || match_length == 3);

        for i in 0..match_length {
            if seq1_is_rc {
                // seq1 (query) was reverse complemented for alignment
                // So the alignment coordinates are in RC space
                let rc_local_pos = seq1_local_start + i;
                // Convert back to forward strand coordinates
                let forward_local_pos = seq1_len - 1 - rc_local_pos;
                let seq1_global_offset = seq1_offset + forward_local_pos;

                // Create reverse orientation for seq1
                let pos1_rev = make_pos(seq1_global_offset, true);

                // seq2 is always forward in this case
                let pos2_fwd = make_pos(seq2_offset + seq2_local_start + i, false);

                if debug || (seq1_offset == 100 && seq2_offset == 200 && match_length == 3) {
                    eprintln!("  [BIDIRECTED_UF_DEBUG] RC transform: rc_local[{}] = {} -> forward_local[{}] = {} -> global[{}] = {}",
                        i, rc_local_pos, i, forward_local_pos, i, seq1_global_offset);
                    eprintln!("    Unite: pos1_rev={} (make_pos({}, true)) <-> pos2_fwd={} (make_pos({}, false))", 
                        pos1_rev, seq1_global_offset, pos2_fwd, seq2_offset + seq2_local_start + i);
                }

                // Unite seq1 reverse with seq2 forward
                self.unite(pos1_rev, pos2_fwd);
            } else {
                // Normal forward-forward alignment
                let pos1_fwd = make_pos(seq1_offset + seq1_local_start + i, false);
                let pos2_fwd = make_pos(seq2_offset + seq2_local_start + i, false);
                self.unite(pos1_fwd, pos2_fwd);
            }
        }
    }

    /// Unite positions when seq2 is reverse complemented (for test compatibility)
    /// This is used by test_graph_construction.rs which has different semantics
    #[allow(clippy::too_many_arguments)]
    pub fn unite_matching_region_seq2_rc(
        &self,
        seq1_offset: usize,
        seq2_offset: usize,
        seq1_local_start: usize,
        seq2_local_start: usize,
        match_length: usize,
        seq2_is_rc: bool,
        seq2_len: usize,
    ) {
        for i in 0..match_length {
            if seq2_is_rc {
                // seq1[i] matches seq2[len-1-i] in reverse orientation
                let pos1_fwd = make_pos(seq1_offset + seq1_local_start + i, false);
                let seq2_rc_pos = seq2_len - 1 - (seq2_local_start + i);
                let pos2_rev = make_pos(seq2_offset + seq2_rc_pos, true);

                // Unite seq1 forward with seq2 reverse
                self.unite(pos1_fwd, pos2_rev);
            } else {
                // Normal forward-forward alignment
                let pos1_fwd = make_pos(seq1_offset + seq1_local_start + i, false);
                let pos2_fwd = make_pos(seq2_offset + seq2_local_start + i, false);
                self.unite(pos1_fwd, pos2_fwd);
            }
        }
    }

    /// Get the underlying UFRush for compatibility
    pub fn inner(&self) -> &Arc<UFRush> {
        &self.uf
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_operations() {
        let uf = BidirectedUnionFind::new(1000);

        let pos1 = make_pos(100, false);
        let pos2 = make_pos(200, false);
        let pos3 = make_pos(100, true); // Same offset, different orientation

        // Initially, positions should not be connected
        assert!(!uf.same(pos1, pos2));
        assert!(!uf.same(pos1, pos3));

        // Unite forward positions
        uf.unite(pos1, pos2);
        assert!(uf.same(pos1, pos2));
        assert!(!uf.same(pos1, pos3)); // Different orientation still separate

        // Unite different orientations of same offset
        uf.unite(pos1, pos3);
        assert!(uf.same(pos1, pos3));
        assert!(uf.same(pos2, pos3)); // Transitivity
    }

    #[test]
    fn test_simple_rc_unite() {
        // Need bigger capacity for the test values
        let uf = BidirectedUnionFind::new(1000);

        // Simple test: unite two positions
        let p1 = make_pos(139, true); // 279
        let p2 = make_pos(215, false); // 430

        println!("Before unite: p1={}, p2={}", p1, p2);
        println!(
            "UFRush capacity should handle max({}, {}) = {}",
            p1,
            p2,
            std::cmp::max(p1, p2)
        );

        uf.unite(p1, p2);

        let r1 = uf.find(p1);
        let r2 = uf.find(p2);
        println!("After unite: find(p1)={}, find(p2)={}", r1, r2);

        assert_eq!(r1, r2, "Representatives should be equal");
        assert!(uf.same(p1, p2), "Positions should be united");
    }

    #[test]
    fn test_unite_matching_region_forward() {
        let uf = BidirectedUnionFind::new(1000);

        // Forward-forward alignment
        uf.unite_matching_region(
            100,   // seq1_offset
            200,   // seq2_offset
            10,    // seq1_local_start
            15,    // seq2_local_start
            5,     // match_length
            false, // seq1_is_rc
            100,   // seq1_len (unused for forward)
        );

        // Check that corresponding positions are united
        let pos1_start = make_pos(110, false); // 100 + 10
        let pos2_start = make_pos(215, false); // 200 + 15
        assert!(uf.same(pos1_start, pos2_start));

        let pos1_end = make_pos(114, false); // 100 + 10 + 4
        let pos2_end = make_pos(219, false); // 200 + 15 + 4
        assert!(uf.same(pos1_end, pos2_end));
    }

    #[test]
    fn test_unite_matching_region_reverse() {
        let uf = BidirectedUnionFind::new(1000);

        // Debug: test direct unite first
        let pos1 = make_pos(139, true);
        let pos2 = make_pos(215, false);
        println!("Direct unite test:");
        println!("pos1: {} (raw)", pos1);
        println!("pos2: {} (raw)", pos2);
        uf.unite(pos1, pos2);
        println!("After unite - same? {}", uf.same(pos1, pos2));
        assert!(uf.same(pos1, pos2), "Direct unite failed");

        // Now test the full method
        let uf2 = BidirectedUnionFind::new(1000);

        // Reverse-forward alignment (seq1/query was reverse complemented)
        uf2.unite_matching_region(
            100,  // seq1_offset
            200,  // seq2_offset
            10,   // seq1_local_start (in RC coordinates)
            15,   // seq2_local_start
            3,    // match_length
            true, // seq1_is_rc
            50,   // seq_len (seq1_len since seq1_is_rc is true)
        );

        // seq1 RC position 10 corresponds to forward position 50-1-10 = 39
        // So seq1 position 139 (100+39) with reverse orientation
        let pos1_start = make_pos(139, true); // 100 + (50-1-10) = 139, reverse
        let pos2_start = make_pos(215, false); // 200 + 15 = 215, forward
        assert!(
            uf2.same(pos1_start, pos2_start),
            "First position not united"
        );

        // seq1 RC position 12 corresponds to forward position 50-1-12 = 37
        let pos1_end = make_pos(137, true); // 100 + (50-1-12) = 137, reverse
        let pos2_end = make_pos(217, false); // 200 + 15 + 2 = 217, forward
        assert!(uf2.same(pos1_end, pos2_end), "Last position not united");
    }
}
