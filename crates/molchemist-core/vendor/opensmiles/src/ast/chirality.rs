use std::fmt;

/// Chirality specification per the OpenSMILES specification.
///
/// Uses `#[repr(u8)]` starting at 1 so that `Option<Chirality>` benefits
/// from niche optimization and occupies exactly **1 byte** (0 = `None`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
#[rustfmt::skip]
pub enum Chirality {
    // Tetrahedral
    TH1 = 1, TH2,
    // Allenal
    AL1, AL2,
    // Square Planar
    SP1, SP2, SP3,
    // Trigonal Bipyramidal (1-20)
    TB1, TB2, TB3, TB4, TB5, TB6, TB7, TB8, TB9, TB10,
    TB11, TB12, TB13, TB14, TB15, TB16, TB17, TB18, TB19, TB20,
    // Octahedral (1-30)
    OH1, OH2, OH3, OH4, OH5, OH6, OH7, OH8, OH9, OH10,
    OH11, OH12, OH13, OH14, OH15, OH16, OH17, OH18, OH19, OH20,
    OH21, OH22, OH23, OH24, OH25, OH26, OH27, OH28, OH29, OH30,
}

impl Chirality {
    /// Construct a trigonal-bipyramidal variant from an index (1–20).
    pub fn tb(n: u8) -> Option<Self> {
        const TB: [Chirality; 20] = [
            Chirality::TB1,
            Chirality::TB2,
            Chirality::TB3,
            Chirality::TB4,
            Chirality::TB5,
            Chirality::TB6,
            Chirality::TB7,
            Chirality::TB8,
            Chirality::TB9,
            Chirality::TB10,
            Chirality::TB11,
            Chirality::TB12,
            Chirality::TB13,
            Chirality::TB14,
            Chirality::TB15,
            Chirality::TB16,
            Chirality::TB17,
            Chirality::TB18,
            Chirality::TB19,
            Chirality::TB20,
        ];
        TB.get((n as usize).wrapping_sub(1)).copied()
    }

    /// Construct an octahedral variant from an index (1–30).
    pub fn oh(n: u8) -> Option<Self> {
        const OH: [Chirality; 30] = [
            Chirality::OH1,
            Chirality::OH2,
            Chirality::OH3,
            Chirality::OH4,
            Chirality::OH5,
            Chirality::OH6,
            Chirality::OH7,
            Chirality::OH8,
            Chirality::OH9,
            Chirality::OH10,
            Chirality::OH11,
            Chirality::OH12,
            Chirality::OH13,
            Chirality::OH14,
            Chirality::OH15,
            Chirality::OH16,
            Chirality::OH17,
            Chirality::OH18,
            Chirality::OH19,
            Chirality::OH20,
            Chirality::OH21,
            Chirality::OH22,
            Chirality::OH23,
            Chirality::OH24,
            Chirality::OH25,
            Chirality::OH26,
            Chirality::OH27,
            Chirality::OH28,
            Chirality::OH29,
            Chirality::OH30,
        ];
        OH.get((n as usize).wrapping_sub(1)).copied()
    }
}

impl fmt::Display for Chirality {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Chirality::TH1 => write!(f, "@"),
            Chirality::TH2 => write!(f, "@@"),
            Chirality::SP1 => write!(f, "@SP1"),
            Chirality::SP2 => write!(f, "@SP2"),
            Chirality::SP3 => write!(f, "@SP3"),
            Chirality::AL1 => write!(f, "@AL1"),
            Chirality::AL2 => write!(f, "@AL2"),
            n => {
                let val = *n as u8;
                if val <= Chirality::TB20 as u8 {
                    let num = val - Chirality::TB1 as u8 + 1;
                    write!(f, "@TB{}", num)
                } else {
                    let num = val - Chirality::OH1 as u8 + 1;
                    write!(f, "@OH{}", num)
                }
            }
        }
    }
}
