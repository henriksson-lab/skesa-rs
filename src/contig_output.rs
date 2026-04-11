/// Contig output formatting matching SKESA's PrintRslt from skesa.cpp.
///
/// Produces FASTA output with contig headers including abundance information,
/// variant annotations, and circular topology markers.
use std::io::Write;

use crate::contig::ContigSequence;
use crate::counter::KmerCount;
use crate::read_holder::ReadHolder;

/// Write contigs in FASTA format matching SKESA's output.
///
/// For each contig above min_contig length:
/// - Header: >Contig_N_abundance [topology=circular]
/// - Sequence: the primary (highest-scoring) variant path
/// - Variants: >Variant_N_for_Contig_M:pos_pos:score
///
/// This simplified version outputs the primary sequence without variant scoring
/// (which requires graph abundance lookups not yet available for the Rust assembler).
pub fn write_contigs<W: Write>(
    out: &mut W,
    contigs: &[ContigSequence],
    min_contig: usize,
) -> std::io::Result<()> {
    let mut num = 0;

    for contig in contigs {
        if contig.len_min() < min_contig {
            continue;
        }

        num += 1;
        let primary = contig.primary_sequence();

        // Simple abundance placeholder (full version needs graph access)
        let abundance = 0.0f64;

        if contig.circular {
            writeln!(out, ">Contig_{}_{}_Circ [topology=circular]", num, abundance)?;
        } else {
            writeln!(out, ">Contig_{}_{}", num, abundance)?;
        }
        writeln!(out, "{}", primary)?;
    }

    Ok(())
}

/// Write contigs with abundance information computed from a k-mer counter.
/// This matches the full C++ PrintRslt output format.
pub fn write_contigs_with_abundance<W: Write>(
    out: &mut W,
    contigs: &[ContigSequence],
    kmers: &KmerCount,
    kmer_len: usize,
    min_contig: usize,
) -> std::io::Result<()> {
    let mut num = 0;

    let mut sorted_contigs: Vec<&ContigSequence> = contigs.iter().collect();
    sorted_contigs.sort();

    for contig in sorted_contigs {
        if contig.len_min() < min_contig {
            continue;
        }

        num += 1;

        // Compute the primary sequence and its abundance
        let mut primary = contig.primary_sequence();
        if contig.circular {
            // For circular contigs, append kmer_len-1 bases for abundance calculation
            let overlap: String = primary.chars().take(kmer_len - 1).collect();
            primary.push_str(&overlap);
        }

        // Calculate average k-mer abundance along the primary sequence
        let abundance = compute_sequence_abundance(&primary, kmers, kmer_len);

        if contig.circular {
            // Remove the overlap from displayed sequence
            let display_len = primary.len() - (kmer_len - 1);
            let display_seq = &primary[..display_len];
            writeln!(out, ">Contig_{}_{}_Circ [topology=circular]", num, abundance)?;
            writeln!(out, "{}", display_seq)?;
        } else {
            writeln!(out, ">Contig_{}_{}", num, abundance)?;
            writeln!(out, "{}", primary)?;
        }

        // TODO: Output variant annotations
    }

    Ok(())
}

/// Compute average k-mer abundance for a sequence
fn compute_sequence_abundance(seq: &str, kmers: &KmerCount, kmer_len: usize) -> f64 {
    if seq.len() < kmer_len {
        return 0.0;
    }

    let mut rh = ReadHolder::new(false);
    rh.push_back_str(seq);

    let mut total_abundance = 0.0;
    let mut count = 0;

    let mut ki = rh.kmer_iter(kmer_len);
    while !ki.at_end() {
        let kmer = ki.get();
        let rkmer = kmer.revcomp(kmer_len);
        let canonical = if kmer < rkmer { kmer } else { rkmer };

        let idx = kmers.find(&canonical);
        if idx < kmers.size() {
            let kmer_count = kmers.get_count(idx);
            total_abundance += (kmer_count & 0xFFFFFFFF) as f64;
        }
        count += 1;
        ki.advance();
    }

    if count > 0 {
        total_abundance / count as f64
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_simple_contigs() {
        let mut contig = ContigSequence::new();
        contig.insert_new_chunk_with("ACGTACGTACGTACGTACGTACGTACGT".chars().collect());

        let mut output = Vec::new();
        write_contigs(&mut output, &[contig], 10).unwrap();
        let result = String::from_utf8(output).unwrap();

        assert!(result.contains(">Contig_1_"));
        assert!(result.contains("ACGTACGTACGTACGTACGTACGTACGT"));
    }

    #[test]
    fn test_write_contigs_respects_min_length() {
        let mut c1 = ContigSequence::new();
        c1.insert_new_chunk_with("ACGTACGT".chars().collect()); // 8bp

        let mut c2 = ContigSequence::new();
        c2.insert_new_chunk_with("A".repeat(200).chars().collect()); // 200bp

        let mut output = Vec::new();
        write_contigs(&mut output, &[c1, c2], 100).unwrap();
        let result = String::from_utf8(output).unwrap();

        // Only the 200bp contig should appear
        assert_eq!(result.lines().filter(|l| l.starts_with('>')).count(), 1);
    }

    #[test]
    fn test_contig_abundance_computation() {
        use crate::reads_getter::ReadsGetter;
        use crate::sorted_counter;

        let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
        let fasta = data_dir.join("small_test.fasta");
        let rg = ReadsGetter::new(&[fasta.to_str().unwrap().to_string()], false).unwrap();
        let reads = rg.reads().to_vec();

        let kmers = sorted_counter::count_kmers_sorted(&reads, 21, 2, true, 32);

        // Pick a sequence that exists in the data
        let (first_kmer, _first_count) = kmers.get_kmer_count(0);
        let seq = first_kmer.to_kmer_string(21);
        let abundance = compute_sequence_abundance(&seq, &kmers, 21);

        // Abundance should be > 0 since this kmer is in the counter
        assert!(abundance > 0.0, "Expected positive abundance, got {}", abundance);
    }
}
