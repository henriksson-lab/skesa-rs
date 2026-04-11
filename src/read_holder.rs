/// Compact DNA sequence storage using 2-bit encoding.
///
/// Port of SKESA's CReadHolder from common_util.hpp.
///
/// Key design: sequences are stored in REVERSE order (last nucleotide first)
/// within the bit-packed storage. This is intentional for k-mer extraction
/// compatibility — the kmer_iterator reads consecutive bits without reversal.
///
/// Encoding: A=0, C=1, T=2, G=3 (matching LargeInt's bin2NT)
use crate::kmer::Kmer;

const BIN2NT: [char; 4] = ['A', 'C', 'T', 'G'];

#[inline]
fn nt_to_bin(c: char) -> u64 {
    match c {
        'A' | 'a' => 0,
        'C' | 'c' => 1,
        'T' | 't' => 2,
        'G' | 'g' => 3,
        _ => 0,
    }
}

#[derive(Clone, Debug)]
pub struct ReadHolder {
    storage: Vec<u64>,
    read_length: Vec<u32>,
    total_seq: usize,
    contains_paired: bool,
}

impl ReadHolder {
    /// Create a new empty ReadHolder.
    ///
    /// # Example
    /// ```
    /// use skesa_rs::read_holder::ReadHolder;
    /// let mut rh = ReadHolder::new(false);
    /// rh.push_back_str("ACGTACGT");
    /// assert_eq!(rh.read_num(), 1);
    /// assert_eq!(rh.total_seq(), 8);
    /// ```
    pub fn new(contains_paired: bool) -> Self {
        ReadHolder {
            storage: Vec::new(),
            read_length: Vec::new(),
            total_seq: 0,
            contains_paired,
        }
    }

    /// Insert a read from a string (stores in reverse for k-mer compatibility)
    pub fn push_back_str(&mut self, read: &str) {
        let mut shift = (self.total_seq * 2) % 64;
        let mut read_len: u32 = 0;
        // Iterate in reverse (matches C++ `read.rbegin()..read.rend()`)
        for c in read.chars().rev() {
            if shift == 0 {
                self.storage.push(0);
            }
            *self.storage.last_mut().unwrap() += nt_to_bin(c) << shift;
            shift = (shift + 2) % 64;
            read_len += 1;
        }
        self.read_length.push(read_len);
        self.total_seq += read_len as usize;
    }

    /// Insert a read from a slice of chars (stores in reverse)
    pub fn push_back_chars(&mut self, read: &[char]) {
        let mut shift = (self.total_seq * 2) % 64;
        let len = read.len() as u32;
        // Iterate in reverse
        for &c in read.iter().rev() {
            if shift == 0 {
                self.storage.push(0);
            }
            *self.storage.last_mut().unwrap() += nt_to_bin(c) << shift;
            shift = (shift + 2) % 64;
        }
        self.read_length.push(len);
        self.total_seq += len as usize;
    }

    /// Swap contents with another ReadHolder
    pub fn swap(&mut self, other: &mut ReadHolder) {
        std::mem::swap(&mut self.storage, &mut other.storage);
        std::mem::swap(&mut self.read_length, &mut other.read_length);
        std::mem::swap(&mut self.total_seq, &mut other.total_seq);
    }

    /// Delete all sequences and release memory
    pub fn clear(&mut self) {
        self.storage.clear();
        self.storage.shrink_to_fit();
        self.read_length.clear();
        self.read_length.shrink_to_fit();
        self.total_seq = 0;
    }

    /// Total nucleotide count
    pub fn total_seq(&self) -> usize {
        self.total_seq
    }

    /// Maximum length of included sequences
    pub fn max_length(&self) -> usize {
        self.read_length.iter().copied().max().unwrap_or(0) as usize
    }

    /// Number of k-mers of given length that could be generated
    pub fn kmer_num(&self, kmer_len: usize) -> usize {
        let mut num = 0;
        for &l in &self.read_length {
            let l = l as usize;
            if l >= kmer_len {
                num += l - kmer_len + 1;
            }
        }
        num
    }

    /// Total number of sequences
    pub fn read_num(&self) -> usize {
        self.read_length.len()
    }

    /// Whether this holder contains paired reads
    pub fn contains_paired(&self) -> bool {
        self.contains_paired
    }

    /// Memory footprint in bytes
    pub fn memory_footprint(&self) -> usize {
        8 * self.storage.capacity() + 4 * self.read_length.capacity()
    }

    /// Reserve storage
    pub fn reserve(&mut self, seq: usize, num: usize) {
        self.storage.reserve(seq / 32 + 1);
        if num > 0 {
            self.read_length.reserve(num);
        }
    }

    /// Shortest sequence length at xx% of total length (NXX statistic)
    pub fn nxx(&self, xx: f64) -> usize {
        let mut lengths: Vec<u32> = self.read_length.clone();
        lengths.sort();
        let mut nxx = 0usize;
        let mut len = 0usize;
        let threshold = (xx * self.total_seq as f64) as usize;
        for j in (0..lengths.len()).rev() {
            nxx = lengths[j] as usize;
            len += lengths[j] as usize;
            if len >= threshold {
                break;
            }
        }
        nxx
    }

    /// N50 statistic
    pub fn n50(&self) -> usize {
        self.nxx(0.5)
    }

    /// Get read length at index
    pub fn read_length_at(&self, idx: usize) -> u32 {
        self.read_length[idx]
    }

    /// Access raw storage (for CopyBits compatibility)
    pub fn storage(&self) -> &[u64] {
        &self.storage
    }

    /// Access raw storage as a flat byte slice (native endian).
    /// Uses the safe `as_flattened` approach via a helper.
    pub fn storage_as_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(self.storage.len() * 8);
        for &word in &self.storage {
            bytes.extend_from_slice(&word.to_ne_bytes());
        }
        bytes
    }

    /// Access raw storage as bytes (for byte-level k-mer extraction).
    /// This returns a reference to the underlying storage reinterpreted as bytes.
    pub fn storage_bytes(&self) -> &[u8] {
        // Safety: u64 slice can always be viewed as u8 slice on any platform.
        // The alignment of u64 is >= alignment of u8.
        let ptr = self.storage.as_ptr() as *const u8;
        let len = self.storage.len() * 8;
        unsafe { std::slice::from_raw_parts(ptr, len) }
    }

    /// Access raw read lengths
    pub fn read_lengths(&self) -> &[u32] {
        &self.read_length
    }

    /// Copy bits from storage to destination.
    /// bit_from/bit_to: source bit range in storage
    /// destination: target buffer
    /// dest_bit_from: starting bit position in destination
    /// dest_size: number of u64 words used in destination
    pub fn copy_bits(
        &self,
        bit_from: usize,
        bit_to: usize,
        destination: &mut [u64],
        dest_bit_from: usize,
        dest_size: usize,
    ) {
        if bit_to <= bit_from {
            return;
        }

        let mut word = bit_from / 64;
        let last_word = (bit_to - 1) / 64;
        let shift = (bit_from % 64) as u32;
        let mut dest_word = dest_bit_from / 64;
        let mut dest_shift = (dest_bit_from % 64) as u32;

        if shift > 0 {
            let chunk = self.storage[word] >> shift;
            word += 1;
            if dest_shift > 0 {
                destination[dest_word] = destination[dest_word].wrapping_add(chunk << dest_shift);
                if shift <= dest_shift {
                    dest_word += 1;
                }
                if shift < dest_shift && dest_word < dest_size {
                    destination[dest_word] =
                        destination[dest_word].wrapping_add(chunk >> (64 - dest_shift));
                }
            } else {
                destination[dest_word] = chunk;
            }
            dest_shift = (dest_shift + 64 - shift) % 64;
        }

        while word <= last_word {
            if dest_shift > 0 {
                destination[dest_word] =
                    destination[dest_word].wrapping_add(self.storage[word] << dest_shift);
                if dest_word + 1 < dest_size {
                    destination[dest_word + 1] = destination[dest_word + 1]
                        .wrapping_add(self.storage[word] >> (64 - dest_shift));
                }
            } else {
                destination[dest_word] = self.storage[word];
            }
            word += 1;
            dest_word += 1;
        }

        let partial_bits = (dest_bit_from + bit_to - bit_from) % 64;
        if partial_bits > 0 {
            let mask = (1u64 << partial_bits) - 1;
            destination[dest_size - 1] &= mask;
        }
    }

    // ── Iterators ──

    /// Create a k-mer iterator starting from the beginning
    pub fn kmer_iter(&self, kmer_len: usize) -> KmerIterator<'_> {
        let mut iter = KmerIterator {
            holder: self,
            read: 0,
            position: 0,
            kmer_len: kmer_len as u32,
            position_in_read: 0,
        };
        iter.skip_short_reads();
        iter
    }

    /// Create a k-mer end sentinel
    pub fn kmer_end(&self) -> KmerIterator<'_> {
        KmerIterator {
            holder: self,
            read: self.read_length.len(),
            position: 2 * self.total_seq,
            kmer_len: 0,
            position_in_read: 0,
        }
    }

    /// Create a string iterator from the beginning
    pub fn string_iter(&self) -> StringIterator<'_> {
        StringIterator {
            holder: self,
            position: 0,
            read: 0,
        }
    }

    /// Create a string end sentinel
    pub fn string_end(&self) -> StringIterator<'_> {
        StringIterator {
            holder: self,
            position: 2 * self.total_seq,
            read: self.read_length.len(),
        }
    }
}

/// Iterator over k-mers in the read holder
#[derive(Clone)]
pub struct KmerIterator<'a> {
    holder: &'a ReadHolder,
    read: usize,
    position: usize,          // bit position in concatenated storage
    kmer_len: u32,
    position_in_read: u32,    // nucleotide position within current read
}

impl<'a> KmerIterator<'a> {
    /// Dereference: extract the k-mer at current position
    pub fn get(&self) -> Kmer {
        let kmer_len = self.kmer_len as usize;
        // Fast path for precision=1 (kmer_len <= 32): extract single u64 without allocation
        if kmer_len <= 32 {
            let val = self.get_val_p1();
            return Kmer::from_u64(kmer_len, val);
        }
        let num_words = (2 * kmer_len).div_ceil(64);
        let mut buf = vec![0u64; num_words];
        self.holder
            .copy_bits(self.position, self.position + 2 * kmer_len, &mut buf, 0, num_words);
        let mut kmer = Kmer::zero(kmer_len);
        kmer.copy_words_from(&buf);
        kmer
    }

    /// Fast extraction of single-word k-mer value (kmer_len <= 32).
    /// Extracts bits directly from storage without heap allocation.
    #[inline]
    fn get_val_p1(&self) -> u64 {
        let kmer_len = self.kmer_len as usize;
        let bit_len = 2 * kmer_len;
        let word_idx = self.position / 64;
        let bit_offset = self.position % 64;
        let storage = self.holder.storage();

        let mut val = storage[word_idx] >> bit_offset;
        if bit_offset + bit_len > 64 && word_idx + 1 < storage.len() {
            val |= storage[word_idx + 1] << (64 - bit_offset);
        }
        // Mask to kmer_len * 2 bits
        if bit_len < 64 {
            val &= (1u64 << bit_len) - 1;
        }
        val
    }

    /// Advance to next k-mer
    pub fn advance(&mut self) {
        let kmer_len = self.kmer_len as usize;
        if self.position == 2 * (self.holder.total_seq - kmer_len) {
            self.position = 2 * self.holder.total_seq;
            return;
        }

        self.position += 2;
        self.position_in_read += 1;
        let read_len = self.holder.read_length[self.read];
        if self.position_in_read == read_len - self.kmer_len + 1 {
            self.position += 2 * (kmer_len - 1);
            self.read += 1;
            self.position_in_read = 0;
            self.skip_short_reads();
        }
    }

    /// Check if at end
    pub fn at_end(&self) -> bool {
        self.position >= 2 * self.holder.total_seq
    }

    fn skip_short_reads(&mut self) {
        let kmer_len = self.kmer_len as usize;
        while self.position < 2 * self.holder.total_seq
            && self.read < self.holder.read_length.len()
            && (self.holder.read_length[self.read] as usize) < kmer_len
        {
            self.position += 2 * self.holder.read_length[self.read] as usize;
            self.read += 1;
        }
    }
}

impl<'a> PartialEq for KmerIterator<'a> {
    fn eq(&self, other: &Self) -> bool {
        std::ptr::eq(self.holder, other.holder) && self.position == other.position
    }
}

/// Iterator over reads as strings
#[derive(Clone)]
pub struct StringIterator<'a> {
    holder: &'a ReadHolder,
    position: usize,
    read: usize,
}

impl<'a> StringIterator<'a> {
    /// Dereference: extract the read as a string
    pub fn get(&self) -> String {
        let read_length = self.holder.read_length[self.read] as usize;
        let mut read = String::with_capacity(read_length);
        // Stored in reverse, so read from high position backward
        let mut position = self.position + 2 * (read_length - 1);
        for _ in 0..read_length {
            let idx = (self.holder.storage[position / 64] >> (position % 64)) & 3;
            read.push(BIN2NT[idx as usize]);
            position = position.wrapping_sub(2);
        }
        read
    }

    /// Get read length
    pub fn read_len(&self) -> usize {
        self.holder.read_length[self.read] as usize
    }

    /// Advance to next read
    pub fn advance(&mut self) {
        if self.read >= self.holder.read_length.len() {
            return;
        }
        self.position += 2 * self.holder.read_length[self.read] as usize;
        self.read += 1;
    }

    /// Check if at end
    pub fn at_end(&self) -> bool {
        self.read >= self.holder.read_length.len()
    }

    /// Whether the holder has paired reads
    pub fn has_mate(&self) -> bool {
        self.holder.contains_paired
    }

    /// Pair type: 0=single, 1=first mate, 2=second mate
    pub fn pair_type(&self) -> u8 {
        if !self.holder.contains_paired {
            0 // single
        } else if self.read % 2 == 1 {
            2 // second mate (odd index)
        } else {
            1 // first mate (even index)
        }
    }

    /// Get the mate iterator (undefined if not paired)
    pub fn get_mate(&self) -> StringIterator<'a> {
        if self.read % 2 == 1 {
            // odd -> mate is previous
            StringIterator {
                holder: self.holder,
                position: self.position - 2 * self.holder.read_length[self.read - 1] as usize,
                read: self.read - 1,
            }
        } else {
            // even -> mate is next
            StringIterator {
                holder: self.holder,
                position: self.position + 2 * self.holder.read_length[self.read] as usize,
                read: self.read + 1,
            }
        }
    }

    /// Get a k-mer iterator for this read
    pub fn kmers_for_read(&self, kmer_len: usize) -> KmerIterator<'a> {
        if kmer_len <= self.holder.read_length[self.read] as usize {
            let mut iter = KmerIterator {
                holder: self.holder,
                read: self.read,
                position: self.position,
                kmer_len: kmer_len as u32,
                position_in_read: 0,
            };
            iter.skip_short_reads();
            iter
        } else {
            self.holder.kmer_end()
        }
    }
}

impl<'a> PartialEq for StringIterator<'a> {
    fn eq(&self, other: &Self) -> bool {
        std::ptr::eq(self.holder, other.holder) && self.read == other.read
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_push_back_and_read() {
        let mut rh = ReadHolder::new(false);
        rh.push_back_str("ACGT");
        assert_eq!(rh.total_seq(), 4);
        assert_eq!(rh.read_num(), 1);

        let mut si = rh.string_iter();
        assert!(!si.at_end());
        assert_eq!(si.get(), "ACGT");
        si.advance();
        assert!(si.at_end());
    }

    #[test]
    fn test_multiple_reads() {
        let mut rh = ReadHolder::new(false);
        rh.push_back_str("ACGT");
        rh.push_back_str("TTGG");
        rh.push_back_str("CCAA");
        assert_eq!(rh.total_seq(), 12);
        assert_eq!(rh.read_num(), 3);

        let mut si = rh.string_iter();
        assert_eq!(si.get(), "ACGT");
        si.advance();
        assert_eq!(si.get(), "TTGG");
        si.advance();
        assert_eq!(si.get(), "CCAA");
        si.advance();
        assert!(si.at_end());
    }

    #[test]
    fn test_kmer_iterator() {
        let mut rh = ReadHolder::new(false);
        rh.push_back_str("ACGTACGT"); // 8bp read

        let mut ki = rh.kmer_iter(4);
        let mut kmers = Vec::new();
        while !ki.at_end() {
            kmers.push(ki.get().to_kmer_string(4));
            ki.advance();
        }
        // 8 - 4 + 1 = 5 kmers
        // K-mers come in reverse order because reads are stored reversed
        assert_eq!(kmers.len(), 5);
        assert_eq!(kmers[0], "ACGT");
        assert_eq!(kmers[1], "TACG");
        assert_eq!(kmers[2], "GTAC");
        assert_eq!(kmers[3], "CGTA");
        assert_eq!(kmers[4], "ACGT");
    }

    #[test]
    fn test_kmer_across_reads() {
        let mut rh = ReadHolder::new(false);
        rh.push_back_str("ACGT"); // 4bp, only 1 kmer of length 4
        rh.push_back_str("TTGG"); // 4bp, only 1 kmer of length 4

        let mut ki = rh.kmer_iter(4);
        let mut kmers = Vec::new();
        while !ki.at_end() {
            kmers.push(ki.get().to_kmer_string(4));
            ki.advance();
        }
        assert_eq!(kmers.len(), 2);
        assert_eq!(kmers[0], "ACGT");
        assert_eq!(kmers[1], "TTGG");
    }

    #[test]
    fn test_skip_short_reads() {
        let mut rh = ReadHolder::new(false);
        rh.push_back_str("AC"); // too short for kmer_len=4
        rh.push_back_str("ACGTACGT"); // long enough

        let mut ki = rh.kmer_iter(4);
        let mut count = 0;
        while !ki.at_end() {
            ki.advance();
            count += 1;
        }
        assert_eq!(count, 5); // only from the second read
    }

    #[test]
    fn test_max_length() {
        let mut rh = ReadHolder::new(false);
        rh.push_back_str("AC");
        rh.push_back_str("ACGTACGT");
        rh.push_back_str("ACGT");
        assert_eq!(rh.max_length(), 8);
    }

    #[test]
    fn test_kmer_num() {
        let mut rh = ReadHolder::new(false);
        rh.push_back_str("ACGTACGT"); // 5 kmers of length 4
        rh.push_back_str("ACGT");      // 1 kmer
        rh.push_back_str("AC");        // 0 kmers
        assert_eq!(rh.kmer_num(4), 6);
    }

    #[test]
    fn test_paired_reads() {
        let mut rh = ReadHolder::new(true);
        rh.push_back_str("ACGT"); // first mate
        rh.push_back_str("TTGG"); // second mate

        let si = rh.string_iter();
        assert_eq!(si.pair_type(), 1); // first mate
        let mate = si.get_mate();
        assert_eq!(mate.pair_type(), 2); // second mate
        assert_eq!(mate.get(), "TTGG");
    }

    #[test]
    fn test_long_read_crossing_word_boundary() {
        // A read longer than 32 nucleotides crosses u64 word boundaries
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 36bp
        let mut rh = ReadHolder::new(false);
        rh.push_back_str(seq);
        assert_eq!(rh.total_seq(), 36);

        let mut si = rh.string_iter();
        assert_eq!(si.get(), seq);
        si.advance();
        assert!(si.at_end());

        // Verify k-mer extraction (first k-mer from reversed storage = last 21bp of original)
        let ki = rh.kmer_iter(21);
        let first_kmer = ki.get().to_kmer_string(21);
        assert_eq!(first_kmer, &seq[seq.len()-21..]);
    }

    #[test]
    fn test_clear() {
        let mut rh = ReadHolder::new(false);
        rh.push_back_str("ACGT");
        assert_eq!(rh.read_num(), 1);
        rh.clear();
        assert_eq!(rh.read_num(), 0);
        assert_eq!(rh.total_seq(), 0);
    }
}
