//! NCBI genetic code tables for codon-to-amino-acid translation.
//! Port of SKESA's genetic_code.hpp.

/// A genetic code table entry: (amino_acids, starts_stops, code_number, name)
pub struct TableInfo {
    pub amino_acids: &'static str,  // 64 characters: amino acid for each codon
    pub starts_stops: &'static str, // 64 characters: M for start, * for stop
    pub code_number: u32,
    pub name: &'static str,
}

/// Standard codon order: T=0, C=1, A=2, G=3
/// Codon index = first*16 + second*4 + third
const CODON_CHARS: [char; 4] = ['T', 'C', 'A', 'G'];

/// SKESA codon order (reverse): A=0, C=1, T=2, G=3 (matching bin2NT)
const BIN2NT: [char; 4] = ['A', 'C', 'T', 'G'];

/// NCBI genetic code table for codon translation.
///
/// # Example
/// ```
/// use skesa_rs::genetic_code::GeneticCode;
/// let gc = GeneticCode::new(1).unwrap(); // Standard code
/// assert_eq!(gc.translate_codon("ATG"), 'M');
/// assert_eq!(gc.translate("ATGGCTAAA", false), "MAK");
/// ```
pub struct GeneticCode {
    table_idx: usize,
    skesa_translation: Vec<char>,
    skesa_starts: Vec<u8>,
    starts: Vec<String>,
}

impl GeneticCode {
    /// Create a GeneticCode for the given NCBI code number (default: 1 = Standard).
    pub fn new(code_number: u32) -> Result<Self, String> {
        let tables = all_tables();
        let table_idx = tables
            .iter()
            .position(|t| t.code_number == code_number)
            .ok_or_else(|| format!("Unknown genetic code: {}", code_number))?;

        let table = &tables[table_idx];

        // Build start codons
        let starts_bytes = table.starts_stops.as_bytes();
        let mut starts = Vec::new();
        for i in 0..64 {
            if starts_bytes[i] == b'M' {
                starts.push(int_to_codon(i));
            }
        }

        // Build SKESA translation (SKESA uses bin2NT order: A=0,C=1,T=2,G=3)
        let aa_bytes = table.amino_acids.as_bytes();
        let mut skesa_translation = vec!['X'; 64];
        let mut skesa_starts = Vec::new();
        for i in 0..64u8 {
            let mut codon = String::with_capacity(3);
            codon.push(BIN2NT[(i >> 4) as usize]);
            codon.push(BIN2NT[((i >> 2) & 3) as usize]);
            codon.push(BIN2NT[(i & 3) as usize]);
            skesa_translation[i as usize] = aa(aa_bytes, &codon);
            if is_start_codon(starts_bytes, &codon) {
                skesa_starts.push(i);
            }
        }

        Ok(GeneticCode {
            table_idx,
            skesa_translation,
            skesa_starts,
            starts,
        })
    }

    /// Translate a codon string to an amino acid character
    pub fn translate_codon(&self, codon: &str) -> char {
        let tables = all_tables();
        aa(tables[self.table_idx].amino_acids.as_bytes(), codon)
    }

    /// Translate a SKESA-encoded codon (6-bit integer) to amino acid
    pub fn translate_skesa_codon(&self, codon: u8) -> char {
        self.skesa_translation[codon as usize]
    }

    /// Translate a nucleotide sequence to protein
    pub fn translate(&self, nuc: &str, use_alts: bool) -> String {
        let bytes = nuc.as_bytes();
        let mut prot = String::new();
        let mut p = 0;

        if use_alts && bytes.len() >= 3 {
            let codon = &nuc[0..3];
            if self.is_start(codon) {
                prot.push('M');
                p = 3;
            }
        }

        while p + 3 <= bytes.len() {
            prot.push(self.translate_codon(&nuc[p..p + 3]));
            p += 3;
        }
        prot
    }

    /// Check if a codon is a start codon
    pub fn is_start(&self, codon: &str) -> bool {
        self.starts.iter().any(|s| s == codon)
    }

    /// Check if a SKESA-encoded codon is a start
    pub fn is_start_skesa(&self, codon: u8) -> bool {
        self.skesa_starts.contains(&codon)
    }

    /// Get all codons encoding a given amino acid
    pub fn codons_for(&self, amino_acid: char) -> Vec<String> {
        let tables = all_tables();
        let aa_bytes = tables[self.table_idx].amino_acids.as_bytes();
        let aa_upper = amino_acid.to_ascii_uppercase() as u8;
        (0..64)
            .filter(|&i| aa_bytes[i] == aa_upper)
            .map(int_to_codon)
            .collect()
    }

    /// Get all stop codons
    pub fn stops(&self) -> Vec<String> {
        let tables = all_tables();
        let ss = tables[self.table_idx].starts_stops.as_bytes();
        (0..64)
            .filter(|&i| ss[i] == b'*')
            .map(int_to_codon)
            .collect()
    }

    /// Get the code name
    pub fn name(&self) -> &'static str {
        all_tables()[self.table_idx].name
    }

    /// All recognised start codons for this genetic code. Port of
    /// `GeneticCode::Starts` (genetic_code.hpp:115). C++ returns a
    /// `list<string>`; Rust returns a slice into the pre-computed vector.
    pub fn starts(&self) -> &[String] {
        &self.starts
    }

    /// The full table of all 25 NCBI genetic code definitions. Port of
    /// `GeneticCode::AllInfo` (genetic_code.hpp:128). C++ exposes the
    /// internal `array<TableInfo, 25>`; Rust returns a slice of the
    /// equivalent `TableInfo` entries.
    pub fn all_info(&self) -> &'static [TableInfo] {
        all_tables()
    }
}

fn codon_to_int(codon: &str) -> usize {
    let b = codon.as_bytes();
    nt_idx(b[0]) * 16 + nt_idx(b[1]) * 4 + nt_idx(b[2])
}

fn nt_idx(c: u8) -> usize {
    match c.to_ascii_uppercase() {
        b'T' => 0,
        b'C' => 1,
        b'A' => 2,
        b'G' => 3,
        _ => 0,
    }
}

fn int_to_codon(i: usize) -> String {
    let mut s = String::with_capacity(3);
    s.push(CODON_CHARS[i >> 4]);
    s.push(CODON_CHARS[(i >> 2) & 3]);
    s.push(CODON_CHARS[i & 3]);
    s
}

fn aa(aa_table: &[u8], codon: &str) -> char {
    let upper: String = codon.chars().map(|c| c.to_ascii_uppercase()).collect();
    if upper.chars().all(|c| "ACGT".contains(c)) {
        aa_table[codon_to_int(&upper)] as char
    } else {
        'X'
    }
}

fn is_start_codon(starts: &[u8], codon: &str) -> bool {
    let upper: String = codon.chars().map(|c| c.to_ascii_uppercase()).collect();
    if upper.chars().all(|c| "ACGT".contains(c)) {
        starts[codon_to_int(&upper)] == b'M'
    } else {
        false
    }
}

fn all_tables() -> &'static [TableInfo] {
    &GENETIC_CODE_TABLES
}

static GENETIC_CODE_TABLES: [TableInfo; 25] = [
    TableInfo {
        amino_acids: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts_stops: "---M------**--*----M---------------M----------------------------",
        code_number: 1,
        name: "Standard Code",
    },
    TableInfo {
        amino_acids: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        starts_stops: "----------**--------------------MMMM----------**---M------------",
        code_number: 2,
        name: "Vertebrate Mitochondrial Code",
    },
    TableInfo {
        amino_acids: "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts_stops: "----------**----------------------MM---------------M------------",
        code_number: 3,
        name: "Yeast Mitochondrial Code",
    },
    TableInfo {
        amino_acids: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts_stops: "--MM------**-------M------------MMMM---------------M------------",
        code_number: 4,
        name: "Mold/Protozoan/Coelenterate Mitochondrial",
    },
    TableInfo {
        amino_acids: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
        starts_stops: "---M------**--------------------MMMM---------------M------------",
        code_number: 5,
        name: "Invertebrate Mitochondrial Code",
    },
    TableInfo {
        amino_acids: "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts_stops: "--------------*--------------------M----------------------------",
        code_number: 6,
        name: "Ciliate/Dasycladacean/Hexamita Nuclear",
    },
    TableInfo {
        amino_acids: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        starts_stops: "----------**-----------------------M---------------M------------",
        code_number: 9,
        name: "Echinoderm/Flatworm Mitochondrial",
    },
    TableInfo {
        amino_acids: "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts_stops: "----------**-----------------------M----------------------------",
        code_number: 10,
        name: "Euplotid Nuclear Code",
    },
    TableInfo {
        amino_acids: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts_stops: "---M------**--*----M------------MMMM---------------M------------",
        code_number: 11,
        name: "Bacterial/Archaeal/Plant Plastid",
    },
    TableInfo {
        amino_acids: "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts_stops: "----------**--*----M---------------M----------------------------",
        code_number: 12,
        name: "Alternative Yeast Nuclear Code",
    },
    TableInfo {
        amino_acids: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
        starts_stops: "---M------**----------------------MM---------------M------------",
        code_number: 13,
        name: "Ascidian Mitochondrial Code",
    },
    TableInfo {
        amino_acids: "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        starts_stops: "-----------*-----------------------M----------------------------",
        code_number: 14,
        name: "Alternative Flatworm Mitochondrial",
    },
    TableInfo {
        amino_acids: "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts_stops: "----------*---*--------------------M----------------------------",
        code_number: 16,
        name: "Chlorophycean Mitochondrial Code",
    },
    TableInfo {
        amino_acids: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        starts_stops: "----------**-----------------------M---------------M------------",
        code_number: 21,
        name: "Trematode Mitochondrial Code",
    },
    TableInfo {
        amino_acids: "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts_stops: "------*---*---*--------------------M----------------------------",
        code_number: 22,
        name: "Scenedesmus obliquus Mitochondrial",
    },
    TableInfo {
        amino_acids: "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts_stops: "--*-------**--*-----------------M--M---------------M------------",
        code_number: 23,
        name: "Thraustochytrium Mitochondrial Code",
    },
    TableInfo {
        amino_acids: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        starts_stops: "---M------**-------M---------------M---------------M------------",
        code_number: 24,
        name: "Pterobranchia Mitochondrial Code",
    },
    TableInfo {
        amino_acids: "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts_stops: "---M------**-----------------------M---------------M------------",
        code_number: 25,
        name: "Candidate Division SR1/Gracilibacteria",
    },
    TableInfo {
        amino_acids: "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts_stops: "----------**--*----M---------------M----------------------------",
        code_number: 26,
        name: "Pachysolen tannophilus Nuclear",
    },
    TableInfo {
        amino_acids: "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts_stops: "--------------*--------------------M----------------------------",
        code_number: 27,
        name: "Karyorelict Nuclear Code",
    },
    TableInfo {
        amino_acids: "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts_stops: "----------**--*--------------------M----------------------------",
        code_number: 28,
        name: "Condylostoma Nuclear Code",
    },
    TableInfo {
        amino_acids: "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts_stops: "--------------*--------------------M----------------------------",
        code_number: 29,
        name: "Mesodinium Nuclear Code",
    },
    TableInfo {
        amino_acids: "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts_stops: "--------------*--------------------M----------------------------",
        code_number: 30,
        name: "Peritrich Nuclear Code",
    },
    TableInfo {
        amino_acids: "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts_stops: "----------**-----------------------M----------------------------",
        code_number: 31,
        name: "Blastocrithidia Nuclear Code",
    },
    TableInfo {
        amino_acids: "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        starts_stops: "---M-------*-------M---------------M---------------M------------",
        code_number: 33,
        name: "Cephalodiscidae Mitochondrial",
    },
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_code() {
        let gc = GeneticCode::new(1).unwrap();
        assert_eq!(gc.name(), "Standard Code");
        assert_eq!(gc.translate_codon("ATG"), 'M');
        assert_eq!(gc.translate_codon("TAA"), '*');
        assert_eq!(gc.translate_codon("GCT"), 'A'); // Alanine
    }

    #[test]
    fn test_translate() {
        let gc = GeneticCode::new(1).unwrap();
        assert_eq!(gc.translate("ATGGCTAAA", false), "MAK");
        assert_eq!(gc.translate("ATGGCTAAA", true), "MAK");
    }

    #[test]
    fn test_stops() {
        let gc = GeneticCode::new(1).unwrap();
        let stops = gc.stops();
        assert!(stops.contains(&"TAA".to_string()));
        assert!(stops.contains(&"TAG".to_string()));
        assert!(stops.contains(&"TGA".to_string()));
    }

    #[test]
    fn test_start_codons() {
        let gc = GeneticCode::new(1).unwrap();
        assert!(gc.is_start("ATG"));
        assert!(!gc.is_start("AAA"));
    }

    #[test]
    fn test_mitochondrial_code() {
        let gc = GeneticCode::new(2).unwrap();
        assert_eq!(gc.name(), "Vertebrate Mitochondrial Code");
        // In mitochondrial code, TGA codes for W instead of stop
        assert_eq!(gc.translate_codon("TGA"), 'W');
    }

    #[test]
    fn test_invalid_code() {
        assert!(GeneticCode::new(999).is_err());
    }

    #[test]
    fn test_codons_for_amino_acid() {
        let gc = GeneticCode::new(1).unwrap();
        let ala_codons = gc.codons_for('A');
        assert_eq!(ala_codons.len(), 4); // GCT, GCC, GCA, GCG
        for codon in &ala_codons {
            assert_eq!(gc.translate_codon(codon), 'A');
        }
    }
}
