/// Nucleotide-protein alignment utilities.
///
/// Port of SKESA's nuc_prot_align.hpp.
/// Used by SAUTE_PROT for protein-guided assembly.
use crate::genetic_code::GeneticCode;
use crate::glb_align::{self, Cigar, SMatrix};

/// Align a nucleotide sequence against a protein sequence.
/// Translates the nucleotide in all 3 reading frames and aligns each against the protein.
/// Returns the best alignment across all frames.
pub fn nuc_prot_align(
    nuc: &[u8],
    protein: &[u8],
    genetic_code: &GeneticCode,
    gopen: i32,
    gapextend: i32,
) -> Option<(Cigar, usize)> {
    if nuc.len() < 3 || protein.is_empty() {
        return None;
    }

    let blosum = SMatrix::new_blosum62();
    let mut best_score = i32::MIN;
    let mut best_result = None;

    // Try all 3 reading frames
    for frame in 0..3usize {
        if nuc.len() < frame + 3 {
            continue;
        }
        let nuc_slice = &nuc[frame..];
        let translated = genetic_code.translate(
            std::str::from_utf8(nuc_slice).unwrap_or(""),
            false,
        );
        let prot_bytes = translated.as_bytes();

        if prot_bytes.is_empty() {
            continue;
        }

        let cigar = glb_align::lcl_align(prot_bytes, protein, gopen, gapextend, &blosum.matrix);
        let score = cigar.matches(prot_bytes, protein) as i32;

        if score > best_score {
            best_score = score;
            best_result = Some((cigar, frame));
        }
    }

    best_result
}

/// Translate a nucleotide sequence to protein in all 6 frames (3 forward + 3 reverse).
pub fn six_frame_translation(
    nuc: &str,
    genetic_code: &GeneticCode,
) -> Vec<(String, usize, bool)> {
    // (protein, frame, is_reverse)
    let mut results = Vec::new();

    // Forward frames
    for frame in 0..3 {
        if nuc.len() >= frame + 3 {
            let translated = genetic_code.translate(&nuc[frame..], false);
            if !translated.is_empty() {
                results.push((translated, frame, false));
            }
        }
    }

    // Reverse complement frames
    let rc: String = nuc.chars().rev().map(crate::model::complement).collect();
    for frame in 0..3 {
        if rc.len() >= frame + 3 {
            let translated = genetic_code.translate(&rc[frame..], false);
            if !translated.is_empty() {
                results.push((translated, frame, true));
            }
        }
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_six_frame_translation() {
        let gc = GeneticCode::new(1).unwrap();
        let nuc = "ATGGCTAAATGA"; // M A K *
        let frames = six_frame_translation(nuc, &gc);
        assert!(!frames.is_empty());
        // Frame 0 forward should give "MAK*"
        let frame0 = frames.iter().find(|(_, f, r)| *f == 0 && !*r);
        assert!(frame0.is_some());
        assert!(frame0.unwrap().0.starts_with("MAK"));
    }

    #[test]
    fn test_nuc_prot_align() {
        let gc = GeneticCode::new(1).unwrap();
        // Nucleotide encoding "MAK" protein
        let nuc = b"ATGGCTAAA";
        let protein = b"MAK";
        let result = nuc_prot_align(nuc, protein, &gc, 5, 2);
        assert!(result.is_some());
        let (_cigar, frame) = result.unwrap();
        assert_eq!(frame, 0);
    }
}
