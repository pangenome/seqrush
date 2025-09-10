use std::num::NonZeroUsize;

use bytes::{Bytes, BytesMut};
use futures::SinkExt;
use tokio::io::AsyncWrite;
use tokio_util::codec::FramedWrite;

use super::{Deflater, Writer};
use crate::{
    r#async::BlockCodec,
    writer::{CompressionLevel, BGZF_EOF, MAX_BUF_SIZE},
};

/// An async BGZF writer builder.
#[derive(Debug, Default)]
pub struct Builder {
    compression_level: Option<CompressionLevel>,
    worker_count: Option<NonZeroUsize>,
}

impl Builder {
    /// Sets a compression level.
    ///
    /// By default, the compression level is set to level 6.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::{self as bgzf, writer::CompressionLevel};
    ///
    /// let builder = bgzf::r#async::writer::Builder::default()
    ///     .set_compression_level(CompressionLevel::best());
    /// ```
    pub fn set_compression_level(mut self, compression_level: CompressionLevel) -> Self {
        self.compression_level = Some(compression_level);
        self
    }

    /// Sets a worker count.
    ///
    /// By default, the worker count is set to the number of available logical CPUs.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    ///
    /// use noodles_bgzf as bgzf;
    ///
    /// let worker_count = NonZeroUsize::try_from(1)?;
    /// let builder = bgzf::r#async::writer::Builder::default()
    ///     .set_worker_count(worker_count);
    /// # Ok::<_, std::num::TryFromIntError>(())
    /// ```
    pub fn set_worker_count(mut self, worker_count: NonZeroUsize) -> Self {
        self.worker_count = Some(worker_count);
        self
    }

    /// Builds an async BGZF writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use tokio::io;
    /// use noodles_bgzf as bgzf;
    /// let writer = bgzf::r#async::writer::Builder::default()
    ///     .build_with_writer(io::sink());
    /// ```
    pub fn build_with_writer<W>(self, writer: W) -> Writer<W>
    where
        W: AsyncWrite,
    {
        let compression_level = self.compression_level.unwrap_or_default();

        let worker_count = self.worker_count.unwrap_or_else(|| {
            // SAFETY: `num_cpus::get` is guaranteed to be non-zero.
            NonZeroUsize::new(num_cpus::get()).unwrap()
        });

        Writer {
            sink: Deflater::new(FramedWrite::new(writer, BlockCodec)).buffer(worker_count.get()),
            buf: BytesMut::with_capacity(MAX_BUF_SIZE),
            eof_buf: Bytes::from_static(BGZF_EOF),
            compression_level: compression_level.into(),
        }
    }
}
