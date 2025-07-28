fn invert_cigar(cigar: &str) -> String {
    let mut result = String::new();
    let mut count = 0;
    
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            count = count * 10 + (ch as usize - '0' as usize);
        } else {
            if count > 0 {
                result.push_str(&count.to_string());
            }
            let op = match ch {
                'I' => 'D',
                'D' => 'I',
                _ => ch,
            };
            result.push(op);
            count = 0;
        }
    }
    
    result
}

fn main() {
    println\!("2I5M3D -> {}", invert_cigar("2I5M3D"));
    println\!("IIMMMDD -> {}", invert_cigar("IIMMMDD"));
    println\!("10M2I3M -> {}", invert_cigar("10M2I3M"));
}
