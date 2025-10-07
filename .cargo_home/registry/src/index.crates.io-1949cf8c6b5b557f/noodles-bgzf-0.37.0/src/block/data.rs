use std::cmp;

/// An uncompressed block data buffer with a cursor.
#[derive(Debug, Default)]
pub struct Data {
    buf: Vec<u8>,
    pos: usize,
}

impl Data {
    pub fn has_remaining(&self) -> bool {
        self.pos < self.len()
    }

    pub fn position(&self) -> usize {
        self.pos
    }

    /// Sets the cursor position.
    ///
    /// The given position must be <= the length of the buffer to be valid.
    pub fn set_position(&mut self, position: usize) {
        self.pos = position;
    }

    pub fn len(&self) -> usize {
        self.buf.len()
    }

    pub fn resize(&mut self, len: usize) {
        self.buf.resize(len, 0);
    }

    /// Moves the cursor from the current position by `amt` bytes.
    ///
    /// This clamps the amount to the length of the buffer.
    pub fn consume(&mut self, amt: usize) {
        self.pos = cmp::min(self.pos + amt, self.buf.len());
    }
}

impl AsRef<[u8]> for Data {
    fn as_ref(&self) -> &[u8] {
        &self.buf[self.pos..]
    }
}

impl AsMut<[u8]> for Data {
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.buf[self.pos..]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_has_remaining() {
        let data = Data::default();
        assert!(!data.has_remaining());

        let mut data = Data::default();
        data.resize(1);
        assert!(data.has_remaining());
    }

    #[test]
    fn test_position() {
        let data = Data::default();
        assert_eq!(data.position(), 0);

        let mut data = Data::default();
        data.resize(1);
        data.set_position(1);
        assert_eq!(data.position(), 1);
    }

    #[test]
    fn test_consume() {
        let mut data = Data::default();
        data.consume(3);
        assert_eq!(data.position(), 0);

        let mut data = Data::default();
        data.resize(5);
        data.consume(3);
        assert_eq!(data.position(), 3);
    }
}
