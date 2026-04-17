/// Contig output formatting for SKESA-compatible `PrintRslt` fixtures.
///
/// Produces FASTA output with contig headers including abundance information,
/// variant annotations, and circular topology markers.
use std::io::Write;

use crate::contig::{ContigSequence, Variation};
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
            writeln!(
                out,
                ">Contig_{}_{}_Circ [topology=circular]",
                num, abundance
            )?;
        } else {
            writeln!(out, ">Contig_{}_{}", num, abundance)?;
        }
        writeln!(out, "{}", primary)?;
    }

    Ok(())
}

/// Write contigs with abundance information computed from a k-mer counter.
/// This matches the tested C++ PrintRslt output format; variant annotations still need broader parity fixtures.
pub fn write_contigs_with_abundance<W: Write>(
    out: &mut W,
    contigs: &[ContigSequence],
    kmers: &KmerCount,
    kmer_len: usize,
    min_contig: usize,
) -> std::io::Result<()> {
    let mut num = 0;

    let mut sorted_contigs: Vec<&ContigSequence> = contigs.iter().collect();
    sorted_contigs.sort_by_key(|contig| contig.primary_sequence());

    for contig in sorted_contigs {
        if contig.len_min() < min_contig {
            continue;
        }

        num += 1;

        let scored_chunks = score_contig_chunks(contig, kmers, kmer_len);
        let primary = primary_from_scored_chunks(&scored_chunks);
        let mut abundance_sequence = primary.clone();
        if contig.circular {
            let overlap: String = abundance_sequence.chars().take(kmer_len - 1).collect();
            abundance_sequence.push_str(&overlap);
        }

        let abundance = compute_sequence_abundance(&abundance_sequence, kmers, kmer_len);

        if contig.circular {
            writeln!(
                out,
                ">Contig_{}_{}_Circ [topology=circular]",
                num,
                format_cpp_float(abundance)
            )?;
        } else {
            writeln!(out, ">Contig_{}_{}", num, format_cpp_float(abundance))?;
        }
        writeln!(out, "{}", primary)?;
        write_variant_annotations(out, num, &scored_chunks)?;
    }

    Ok(())
}

#[derive(Clone, Debug)]
struct ScoredVariant {
    score: f64,
    sequence: String,
}

type ScoredChunk = Vec<ScoredVariant>;

fn score_contig_chunks(
    contig: &ContigSequence,
    kmers: &KmerCount,
    kmer_len: usize,
) -> Vec<ScoredChunk> {
    let mut scored = Vec::with_capacity(contig.chunks.len());
    for (chunk_idx, chunk) in contig.chunks.iter().enumerate() {
        if !contig.variable_chunk(chunk_idx) {
            scored.push(vec![ScoredVariant {
                score: 1.0,
                sequence: variation_to_string(&chunk[0]),
            }]);
            continue;
        }

        let left_context = if chunk_idx > 0 {
            let previous = variation_to_string(&contig.chunks[chunk_idx - 1][0]);
            suffix_chars(&previous, kmer_len.saturating_sub(1))
        } else {
            String::new()
        };
        let right_context = if chunk_idx + 1 < contig.chunks.len() {
            let next = variation_to_string(&contig.chunks[chunk_idx + 1][0]);
            prefix_chars(&next, kmer_len.saturating_sub(1))
        } else {
            String::new()
        };

        let mut variants = Vec::with_capacity(chunk.len());
        let mut total = 0.0;
        for variant in chunk {
            let sequence = variation_to_string(variant);
            let mut scored_sequence = String::new();
            scored_sequence.push_str(&left_context);
            scored_sequence.push_str(&sequence);
            scored_sequence.push_str(&right_context);
            let abundance = sum_sequence_abundance(&scored_sequence, kmers, kmer_len);
            total += abundance;
            variants.push(ScoredVariant {
                score: abundance,
                sequence,
            });
        }
        if total > 0.0 {
            for variant in &mut variants {
                variant.score /= total;
            }
        }
        variants.sort_by(|left, right| {
            right
                .score
                .total_cmp(&left.score)
                .then_with(|| right.sequence.cmp(&left.sequence))
        });
        scored.push(variants);
    }
    scored
}

fn primary_from_scored_chunks(chunks: &[ScoredChunk]) -> String {
    let mut primary = String::new();
    for chunk in chunks {
        if let Some(variant) = chunk.first() {
            primary.push_str(&variant.sequence);
        }
    }
    primary
}

fn write_variant_annotations<W: Write>(
    out: &mut W,
    contig_num: usize,
    chunks: &[ScoredChunk],
) -> std::io::Result<()> {
    let mut pos: isize = 0;
    for (chunk_idx, chunk) in chunks.iter().enumerate() {
        let chunk_len = chunk.first().map_or(0, |variant| variant.sequence.len());
        if chunk.len() > 1 {
            let left = if chunk_idx > 0 {
                chunks[chunk_idx - 1]
                    .first()
                    .map_or(0, |variant| variant.sequence.len().min(100))
            } else {
                0
            };
            let right = if chunk_idx + 1 < chunks.len() {
                chunks[chunk_idx + 1]
                    .first()
                    .map_or(0, |variant| variant.sequence.len().min(100))
            } else {
                0
            };
            for (variant_idx, variant) in chunk.iter().enumerate().skip(1) {
                writeln!(
                    out,
                    ">Variant_{}_for_Contig_{}:{}_{}:{}",
                    variant_idx,
                    contig_num,
                    pos - left as isize + 1,
                    pos + chunk_len as isize + right as isize,
                    format_cpp_float(variant.score)
                )?;
                if left > 0 {
                    let previous = &chunks[chunk_idx - 1].first().unwrap().sequence;
                    write!(out, "{}", suffix_chars(previous, left))?;
                }
                write!(out, "{}", variant.sequence)?;
                if right > 0 {
                    let next = &chunks[chunk_idx + 1].first().unwrap().sequence;
                    write!(out, "{}", prefix_chars(next, right))?;
                }
                writeln!(out)?;
            }
        }
        pos += chunk_len as isize;
    }
    Ok(())
}

fn variation_to_string(variant: &Variation) -> String {
    variant.iter().collect()
}

fn prefix_chars(value: &str, count: usize) -> String {
    value.chars().take(count).collect()
}

fn suffix_chars(value: &str, count: usize) -> String {
    let chars: Vec<char> = value.chars().collect();
    chars[chars.len().saturating_sub(count)..].iter().collect()
}

fn sum_sequence_abundance(seq: &str, kmers: &KmerCount, kmer_len: usize) -> f64 {
    if seq.len() < kmer_len {
        return 0.0;
    }

    let mut rh = ReadHolder::new(false);
    rh.push_back_str(seq);
    let mut total_abundance = 0.0;
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
        ki.advance();
    }
    total_abundance
}

pub fn write_all_iterations<W: Write>(
    writer: &mut W,
    result: &crate::assembler::AssemblyResult,
    allow_snps: bool,
    has_seeds: bool,
) -> std::io::Result<()> {
    let seed_sequences: std::collections::HashSet<String> = if has_seeds {
        result
            .all_iterations
            .first()
            .into_iter()
            .flat_map(|contigs| contigs.iter().map(|contig| contig.primary_sequence()))
            .collect()
    } else {
        std::collections::HashSet::new()
    };

    let graph_offset = if has_seeds {
        if let Some(contigs) = result.all_iterations.first() {
            write_iteration_contigs(writer, "Seed", contigs, None, None)?;
        }
        1
    } else {
        0
    };

    for (i, (kmer_len, _)) in result.graphs.iter().enumerate() {
        if let Some(contigs) = result.all_iterations.get(graph_offset + i) {
            write_iteration_contigs(
                writer,
                &format!("kmer{}", kmer_len),
                contigs,
                Some((*kmer_len - 1) as i32),
                Some(&seed_sequences),
            )?;
        }
    }

    if allow_snps {
        let offset = graph_offset + result.graphs.len();
        for (i, (kmer_len, _)) in result.graphs.iter().rev().enumerate() {
            if let Some(contigs) = result.all_iterations.get(offset + i) {
                write_iteration_contigs(
                    writer,
                    &format!("SNP_recovery_kmer{}", kmer_len),
                    contigs,
                    Some((*kmer_len - 1) as i32),
                    Some(&seed_sequences),
                )?;
            }
        }
    }

    Ok(())
}

fn write_iteration_contigs<W: Write>(
    writer: &mut W,
    prefix: &str,
    contigs: &crate::contig::ContigSequenceList,
    repeat_fallback: Option<i32>,
    fallback_exclusions: Option<&std::collections::HashSet<String>>,
) -> std::io::Result<()> {
    let mut contigs = contigs.clone();
    contigs.sort_by_key(|contig| contig.primary_sequence());
    for (i, contig) in contigs.iter().enumerate() {
        let seq = contig.primary_sequence();
        let use_fallback = !fallback_exclusions.is_some_and(|excluded| excluded.contains(&seq));
        let left_repeat = if use_fallback && contig.left_repeat == 0 {
            repeat_fallback.unwrap_or(0)
        } else {
            contig.left_repeat
        };
        let right_repeat = if use_fallback && contig.right_repeat == 0 {
            repeat_fallback.unwrap_or(0)
        } else {
            contig.right_repeat
        };
        writeln!(
            writer,
            ">{}_{} {} {}",
            prefix,
            i + 1,
            left_repeat,
            right_repeat
        )?;
        writeln!(writer, "{}", seq)?;
    }
    Ok(())
}

fn format_cpp_float(value: f64) -> String {
    if value == 0.0 {
        return "0".to_string();
    }

    let abs = value.abs();
    let digits_before = if abs >= 1.0 {
        abs.log10().floor() as i32 + 1
    } else {
        0
    };
    let decimals = (6 - digits_before).max(0) as usize;
    let mut s = format!("{:.*}", decimals, value);
    if s.contains('.') {
        while s.ends_with('0') {
            s.pop();
        }
        if s.ends_with('.') {
            s.pop();
        }
    }
    s
}

/// Compute average k-mer abundance for a sequence
fn compute_sequence_abundance(seq: &str, kmers: &KmerCount, kmer_len: usize) -> f64 {
    if seq.len() < kmer_len {
        return 0.0;
    }

    let count = seq.len() - kmer_len + 1;
    if count > 0 {
        sum_sequence_abundance(seq, kmers, kmer_len) / count as f64
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
    fn test_write_contigs_uses_full_sequence_lines() {
        let sequence = "ACGT".repeat(40);
        let mut contig = ContigSequence::new();
        contig.insert_new_chunk_with(sequence.chars().collect());

        let mut output = Vec::new();
        write_contigs(&mut output, &[contig], 1).unwrap();
        let result = String::from_utf8(output).unwrap();
        let mut lines = result.lines();
        assert!(lines.next().unwrap().starts_with(">Contig_1_"));
        assert_eq!(lines.next().unwrap(), sequence);
        assert!(lines.next().is_none());
    }

    #[test]
    fn test_write_contigs_marks_circular_topology_without_repeating_sequence() {
        let mut contig = ContigSequence::new();
        contig.insert_new_chunk_with("ACGTAC".chars().collect());
        contig.circular = true;

        let mut output = Vec::new();
        write_contigs(&mut output, &[contig], 1).unwrap();
        let result = String::from_utf8(output).unwrap();
        let mut lines = result.lines();

        assert_eq!(
            lines.next().unwrap(),
            ">Contig_1_0_Circ [topology=circular]"
        );
        assert_eq!(lines.next().unwrap(), "ACGTAC");
        assert!(lines.next().is_none());
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

    fn push_canonical_count(kmers: &mut KmerCount, sequence: &str, kmer_len: usize, count: u64) {
        let kmer = crate::kmer::Kmer::from_kmer_str(sequence);
        let rkmer = kmer.revcomp(kmer_len);
        let canonical = if kmer < rkmer { kmer } else { rkmer };
        kmers.push_back(&canonical, count);
    }

    #[test]
    fn test_write_contigs_with_abundance_outputs_variant_annotations() {
        let mut contig = ContigSequence::new();
        contig.insert_new_chunk_with(vec!['A', 'A']);
        contig.chunks.push(vec![vec!['C'], vec!['G']]);
        contig.insert_new_chunk_with(vec!['T', 'T']);

        let mut kmers = KmerCount::new(3);
        for sequence in ["AAG", "AGT", "GTT"] {
            push_canonical_count(&mut kmers, sequence, 3, 10);
        }
        for sequence in ["AAC", "ACT", "CTT"] {
            push_canonical_count(&mut kmers, sequence, 3, 1);
        }
        kmers.sort_and_uniq(0);

        let mut output = Vec::new();
        write_contigs_with_abundance(&mut output, &[contig], &kmers, 3, 1).unwrap();
        let text = String::from_utf8(output).unwrap();
        assert!(
            text.starts_with(
                ">Contig_1_11
AAGTT
"
            ),
            "{text}"
        );
        assert!(
            text.contains(
                ">Variant_1_for_Contig_1:1_5:0.5
AACTT
"
            ),
            "{text}"
        );
    }

    #[test]
    fn test_variant_annotations_follow_score_order_after_primary() {
        let mut contig = ContigSequence::new();
        contig.insert_new_chunk_with(vec!['A', 'A']);
        contig.chunks.push(vec![vec!['C'], vec!['G'], vec!['T']]);
        contig.insert_new_chunk_with(vec!['T', 'T']);

        let mut kmers = KmerCount::new(3);
        for sequence in ["AAG", "AGT", "GTT"] {
            push_canonical_count(&mut kmers, sequence, 3, 10);
        }
        for sequence in ["AAT", "ATT", "TTT"] {
            push_canonical_count(&mut kmers, sequence, 3, 5);
        }
        for sequence in ["AAC", "ACT", "CTT"] {
            push_canonical_count(&mut kmers, sequence, 3, 1);
        }
        kmers.sort_and_uniq(0);

        let mut output = Vec::new();
        write_contigs_with_abundance(&mut output, &[contig], &kmers, 3, 1).unwrap();
        let text = String::from_utf8(output).unwrap();
        let lines: Vec<_> = text.lines().collect();

        assert_eq!(lines[0], ">Contig_1_11");
        assert_eq!(lines[1], "AAGTT");
        assert!(lines[2].starts_with(">Variant_1_for_Contig_1:1_5:"));
        assert_eq!(lines[3], "AACTT");
        assert!(lines[4].starts_with(">Variant_2_for_Contig_1:1_5:"));
        assert_eq!(lines[5], "AATTT");

        let first_score: f64 = lines[2].rsplit(':').next().unwrap().parse().unwrap();
        let second_score: f64 = lines[4].rsplit(':').next().unwrap().parse().unwrap();
        assert!(first_score > second_score, "{text}");
    }

    fn single_contig(sequence: &str) -> ContigSequence {
        let mut contig = ContigSequence::new();
        contig.insert_new_chunk_with(sequence.chars().collect());
        contig
    }

    #[test]
    fn test_write_all_iterations_orders_snp_recovery_by_reverse_graph_order() {
        let result = crate::assembler::AssemblyResult {
            contigs: vec![single_contig("FINAL")],
            all_iterations: vec![
                vec![single_contig("K21")],
                vec![single_contig("K33")],
                vec![single_contig("SNP33")],
                vec![single_contig("SNP21")],
            ],
            graphs: vec![(21, KmerCount::new(21)), (33, KmerCount::new(33))],
            connected_reads: ReadHolder::new(false),
        };

        let mut output = Vec::new();
        write_all_iterations(&mut output, &result, true, false).unwrap();
        let text = String::from_utf8(output).unwrap();

        assert_eq!(
            text,
            ">kmer21_1 20 20\nK21\n\
             >kmer33_1 32 32\nK33\n\
             >SNP_recovery_kmer33_1 32 32\nSNP33\n\
             >SNP_recovery_kmer21_1 20 20\nSNP21\n"
        );
    }

    #[test]
    fn test_write_all_iterations_preserves_explicit_repeat_metadata() {
        let mut contig = ContigSequence::new();
        contig.insert_new_chunk_with("ACGTACGT".chars().collect());
        contig.left_repeat = 3;
        contig.right_repeat = 7;

        let result = crate::assembler::AssemblyResult {
            contigs: vec![contig.clone()],
            all_iterations: vec![vec![contig]],
            graphs: vec![(21, KmerCount::new(21))],
            connected_reads: ReadHolder::new(false),
        };

        let mut output = Vec::new();
        write_all_iterations(&mut output, &result, false, false).unwrap();
        let text = String::from_utf8(output).unwrap();

        assert_eq!(text, ">kmer21_1 3 7\nACGTACGT\n");
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
        assert!(
            abundance > 0.0,
            "Expected positive abundance, got {}",
            abundance
        );
    }
}
