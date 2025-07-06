/// Position type that encodes both offset and orientation
/// The least significant bit (LSB) indicates orientation:
/// - 0 = forward strand (+)
/// - 1 = reverse strand (-)
/// The remaining bits store the offset
pub type Pos = u64;

/// Create a position with the given offset and orientation
#[inline]
pub fn make_pos(offset: usize, is_reverse: bool) -> Pos {
    ((offset as u64) << 1) | (is_reverse as u64)
}

/// Check if a position is on the reverse strand
#[inline]
pub fn is_rev(pos: Pos) -> bool {
    (pos & 1) == 1
}

/// Get the offset from a position
#[inline]
pub fn offset(pos: Pos) -> usize {
    (pos >> 1) as usize
}

/// Increment a position (respecting orientation)
#[inline]
pub fn incr_pos(pos: Pos) -> Pos {
    if is_rev(pos) {
        // On reverse strand, incrementing means going backward
        let off = offset(pos);
        if off > 0 {
            make_pos(off - 1, true)
        } else {
            pos // Can't go below 0
        }
    } else {
        // On forward strand, increment normally
        make_pos(offset(pos) + 1, false)
    }
}

/// Decrement a position (respecting orientation)
#[inline]
pub fn decr_pos(pos: Pos) -> Pos {
    if is_rev(pos) {
        // On reverse strand, decrementing means going forward
        make_pos(offset(pos) + 1, true)
    } else {
        // On forward strand, decrement normally
        let off = offset(pos);
        if off > 0 {
            make_pos(off - 1, false)
        } else {
            pos // Can't go below 0
        }
    }
}

/// Flip the orientation of a position
#[inline]
pub fn flip_orientation(pos: Pos) -> Pos {
    pos ^ 1
}

/// Get orientation as a character ('+' or '-')
#[inline]
pub fn orientation_char(pos: Pos) -> char {
    if is_rev(pos) { '-' } else { '+' }
}

/// Reverse complement a base
#[inline]
pub fn rc_base(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'N' | b'n' => b'N',
        _ => base,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_position_encoding() {
        // Test forward position
        let pos_fwd = make_pos(100, false);
        assert_eq!(offset(pos_fwd), 100);
        assert!(!is_rev(pos_fwd));
        assert_eq!(orientation_char(pos_fwd), '+');
        
        // Test reverse position
        let pos_rev = make_pos(100, true);
        assert_eq!(offset(pos_rev), 100);
        assert!(is_rev(pos_rev));
        assert_eq!(orientation_char(pos_rev), '-');
    }
    
    #[test]
    fn test_position_increment() {
        // Forward strand increment
        let pos = make_pos(10, false);
        let next = incr_pos(pos);
        assert_eq!(offset(next), 11);
        assert!(!is_rev(next));
        
        // Reverse strand increment (goes backward)
        let pos_rev = make_pos(10, true);
        let next_rev = incr_pos(pos_rev);
        assert_eq!(offset(next_rev), 9);
        assert!(is_rev(next_rev));
    }
    
    #[test]
    fn test_position_decrement() {
        // Forward strand decrement
        let pos = make_pos(10, false);
        let prev = decr_pos(pos);
        assert_eq!(offset(prev), 9);
        assert!(!is_rev(prev));
        
        // Reverse strand decrement (goes forward)
        let pos_rev = make_pos(10, true);
        let prev_rev = decr_pos(pos_rev);
        assert_eq!(offset(prev_rev), 11);
        assert!(is_rev(prev_rev));
    }
    
    #[test]
    fn test_flip_orientation() {
        let pos_fwd = make_pos(50, false);
        let pos_rev = flip_orientation(pos_fwd);
        assert_eq!(offset(pos_rev), 50);
        assert!(is_rev(pos_rev));
        
        let back = flip_orientation(pos_rev);
        assert_eq!(pos_fwd, back);
    }
    
    #[test]
    fn test_boundary_conditions() {
        // Test at offset 0
        let pos_zero = make_pos(0, false);
        let dec = decr_pos(pos_zero);
        assert_eq!(offset(dec), 0); // Can't go negative
        
        let pos_zero_rev = make_pos(0, true);
        let inc = incr_pos(pos_zero_rev);
        assert_eq!(offset(inc), 0); // Can't go negative
    }
}