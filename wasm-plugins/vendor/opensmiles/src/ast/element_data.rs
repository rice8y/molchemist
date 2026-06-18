use super::atom::{AtomSymbol, OrganicAtom};

/// Static chemical data for an element.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ElementData {
    pub atomic_number: u8,
    pub standard_mass: f64,
    pub valence_electrons: u8,
}

impl AtomSymbol {
    /// Returns the static chemical data for this element.
    /// Evaluated at compile time when possible.
    pub const fn element_data(&self) -> ElementData {
        match self {
            AtomSymbol::H => ElementData {
                atomic_number: 1,
                standard_mass: 1.008,
                valence_electrons: 1,
            },
            AtomSymbol::He => ElementData {
                atomic_number: 2,
                standard_mass: 4.003,
                valence_electrons: 2,
            },
            AtomSymbol::Li => ElementData {
                atomic_number: 3,
                standard_mass: 6.941,
                valence_electrons: 1,
            },
            AtomSymbol::Be => ElementData {
                atomic_number: 4,
                standard_mass: 9.012,
                valence_electrons: 2,
            },
            AtomSymbol::Ne => ElementData {
                atomic_number: 10,
                standard_mass: 20.180,
                valence_electrons: 8,
            },
            AtomSymbol::Na => ElementData {
                atomic_number: 11,
                standard_mass: 22.990,
                valence_electrons: 1,
            },
            AtomSymbol::Mg => ElementData {
                atomic_number: 12,
                standard_mass: 24.305,
                valence_electrons: 2,
            },
            AtomSymbol::Al => ElementData {
                atomic_number: 13,
                standard_mass: 26.982,
                valence_electrons: 3,
            },
            AtomSymbol::Si => ElementData {
                atomic_number: 14,
                standard_mass: 28.086,
                valence_electrons: 4,
            },
            AtomSymbol::Ar => ElementData {
                atomic_number: 18,
                standard_mass: 39.948,
                valence_electrons: 8,
            },
            AtomSymbol::K => ElementData {
                atomic_number: 19,
                standard_mass: 39.098,
                valence_electrons: 1,
            },
            AtomSymbol::Ca => ElementData {
                atomic_number: 20,
                standard_mass: 40.078,
                valence_electrons: 2,
            },
            AtomSymbol::Sc => ElementData {
                atomic_number: 21,
                standard_mass: 44.956,
                valence_electrons: 3,
            },
            AtomSymbol::Ti => ElementData {
                atomic_number: 22,
                standard_mass: 47.867,
                valence_electrons: 4,
            },
            AtomSymbol::V => ElementData {
                atomic_number: 23,
                standard_mass: 50.942,
                valence_electrons: 5,
            },
            AtomSymbol::Cr => ElementData {
                atomic_number: 24,
                standard_mass: 51.996,
                valence_electrons: 6,
            },
            AtomSymbol::Mn => ElementData {
                atomic_number: 25,
                standard_mass: 54.938,
                valence_electrons: 7,
            },
            AtomSymbol::Fe => ElementData {
                atomic_number: 26,
                standard_mass: 55.845,
                valence_electrons: 8,
            },
            AtomSymbol::Co => ElementData {
                atomic_number: 27,
                standard_mass: 58.933,
                valence_electrons: 9,
            },
            AtomSymbol::Ni => ElementData {
                atomic_number: 28,
                standard_mass: 58.693,
                valence_electrons: 10,
            },
            AtomSymbol::Cu => ElementData {
                atomic_number: 29,
                standard_mass: 63.546,
                valence_electrons: 11,
            },
            AtomSymbol::Zn => ElementData {
                atomic_number: 30,
                standard_mass: 65.38,
                valence_electrons: 2,
            },
            AtomSymbol::Ga => ElementData {
                atomic_number: 31,
                standard_mass: 69.723,
                valence_electrons: 3,
            },
            AtomSymbol::Ge => ElementData {
                atomic_number: 32,
                standard_mass: 72.630,
                valence_electrons: 4,
            },
            AtomSymbol::As => ElementData {
                atomic_number: 33,
                standard_mass: 74.922,
                valence_electrons: 5,
            },
            AtomSymbol::Se => ElementData {
                atomic_number: 34,
                standard_mass: 78.971,
                valence_electrons: 6,
            },
            AtomSymbol::Kr => ElementData {
                atomic_number: 36,
                standard_mass: 83.798,
                valence_electrons: 8,
            },
            AtomSymbol::Rb => ElementData {
                atomic_number: 37,
                standard_mass: 85.468,
                valence_electrons: 1,
            },
            AtomSymbol::Sr => ElementData {
                atomic_number: 38,
                standard_mass: 87.62,
                valence_electrons: 2,
            },
            AtomSymbol::Y => ElementData {
                atomic_number: 39,
                standard_mass: 88.906,
                valence_electrons: 3,
            },
            AtomSymbol::Zr => ElementData {
                atomic_number: 40,
                standard_mass: 91.224,
                valence_electrons: 4,
            },
            AtomSymbol::Nb => ElementData {
                atomic_number: 41,
                standard_mass: 92.906,
                valence_electrons: 5,
            },
            AtomSymbol::Mo => ElementData {
                atomic_number: 42,
                standard_mass: 95.95,
                valence_electrons: 6,
            },
            AtomSymbol::Tc => ElementData {
                atomic_number: 43,
                standard_mass: 97.0,
                valence_electrons: 7,
            },
            AtomSymbol::Ru => ElementData {
                atomic_number: 44,
                standard_mass: 101.07,
                valence_electrons: 8,
            },
            AtomSymbol::Rh => ElementData {
                atomic_number: 45,
                standard_mass: 102.91,
                valence_electrons: 9,
            },
            AtomSymbol::Pd => ElementData {
                atomic_number: 46,
                standard_mass: 106.42,
                valence_electrons: 10,
            },
            AtomSymbol::Ag => ElementData {
                atomic_number: 47,
                standard_mass: 107.87,
                valence_electrons: 11,
            },
            AtomSymbol::Cd => ElementData {
                atomic_number: 48,
                standard_mass: 112.41,
                valence_electrons: 2,
            },
            AtomSymbol::In => ElementData {
                atomic_number: 49,
                standard_mass: 114.82,
                valence_electrons: 3,
            },
            AtomSymbol::Sn => ElementData {
                atomic_number: 50,
                standard_mass: 118.71,
                valence_electrons: 4,
            },
            AtomSymbol::Sb => ElementData {
                atomic_number: 51,
                standard_mass: 121.76,
                valence_electrons: 5,
            },
            AtomSymbol::Te => ElementData {
                atomic_number: 52,
                standard_mass: 127.60,
                valence_electrons: 6,
            },
            AtomSymbol::Xe => ElementData {
                atomic_number: 54,
                standard_mass: 131.29,
                valence_electrons: 8,
            },
            AtomSymbol::Cs => ElementData {
                atomic_number: 55,
                standard_mass: 132.91,
                valence_electrons: 1,
            },
            AtomSymbol::Ba => ElementData {
                atomic_number: 56,
                standard_mass: 137.33,
                valence_electrons: 2,
            },
            AtomSymbol::La => ElementData {
                atomic_number: 57,
                standard_mass: 138.91,
                valence_electrons: 3,
            },
            AtomSymbol::Ce => ElementData {
                atomic_number: 58,
                standard_mass: 140.12,
                valence_electrons: 4,
            },
            AtomSymbol::Pr => ElementData {
                atomic_number: 59,
                standard_mass: 140.91,
                valence_electrons: 5,
            },
            AtomSymbol::Nd => ElementData {
                atomic_number: 60,
                standard_mass: 144.24,
                valence_electrons: 6,
            },
            AtomSymbol::Pm => ElementData {
                atomic_number: 61,
                standard_mass: 145.0,
                valence_electrons: 7,
            },
            AtomSymbol::Sm => ElementData {
                atomic_number: 62,
                standard_mass: 150.36,
                valence_electrons: 8,
            },
            AtomSymbol::Eu => ElementData {
                atomic_number: 63,
                standard_mass: 151.96,
                valence_electrons: 9,
            },
            AtomSymbol::Gd => ElementData {
                atomic_number: 64,
                standard_mass: 157.25,
                valence_electrons: 10,
            },
            AtomSymbol::Tb => ElementData {
                atomic_number: 65,
                standard_mass: 158.93,
                valence_electrons: 11,
            },
            AtomSymbol::Dy => ElementData {
                atomic_number: 66,
                standard_mass: 162.50,
                valence_electrons: 12,
            },
            AtomSymbol::Ho => ElementData {
                atomic_number: 67,
                standard_mass: 164.93,
                valence_electrons: 13,
            },
            AtomSymbol::Er => ElementData {
                atomic_number: 68,
                standard_mass: 167.26,
                valence_electrons: 14,
            },
            AtomSymbol::Tm => ElementData {
                atomic_number: 69,
                standard_mass: 168.93,
                valence_electrons: 15,
            },
            AtomSymbol::Yb => ElementData {
                atomic_number: 70,
                standard_mass: 173.05,
                valence_electrons: 16,
            },
            AtomSymbol::Lu => ElementData {
                atomic_number: 71,
                standard_mass: 174.97,
                valence_electrons: 3,
            },
            AtomSymbol::Hf => ElementData {
                atomic_number: 72,
                standard_mass: 178.49,
                valence_electrons: 4,
            },
            AtomSymbol::Ta => ElementData {
                atomic_number: 73,
                standard_mass: 180.95,
                valence_electrons: 5,
            },
            AtomSymbol::W => ElementData {
                atomic_number: 74,
                standard_mass: 183.84,
                valence_electrons: 6,
            },
            AtomSymbol::Re => ElementData {
                atomic_number: 75,
                standard_mass: 186.21,
                valence_electrons: 7,
            },
            AtomSymbol::Os => ElementData {
                atomic_number: 76,
                standard_mass: 190.23,
                valence_electrons: 8,
            },
            AtomSymbol::Ir => ElementData {
                atomic_number: 77,
                standard_mass: 192.22,
                valence_electrons: 9,
            },
            AtomSymbol::Pt => ElementData {
                atomic_number: 78,
                standard_mass: 195.08,
                valence_electrons: 10,
            },
            AtomSymbol::Au => ElementData {
                atomic_number: 79,
                standard_mass: 196.97,
                valence_electrons: 11,
            },
            AtomSymbol::Hg => ElementData {
                atomic_number: 80,
                standard_mass: 200.59,
                valence_electrons: 2,
            },
            AtomSymbol::Tl => ElementData {
                atomic_number: 81,
                standard_mass: 204.38,
                valence_electrons: 3,
            },
            AtomSymbol::Pb => ElementData {
                atomic_number: 82,
                standard_mass: 207.2,
                valence_electrons: 4,
            },
            AtomSymbol::Bi => ElementData {
                atomic_number: 83,
                standard_mass: 208.98,
                valence_electrons: 5,
            },
            AtomSymbol::Po => ElementData {
                atomic_number: 84,
                standard_mass: 209.0,
                valence_electrons: 6,
            },
            AtomSymbol::At => ElementData {
                atomic_number: 85,
                standard_mass: 210.0,
                valence_electrons: 7,
            },
            AtomSymbol::Rn => ElementData {
                atomic_number: 86,
                standard_mass: 222.0,
                valence_electrons: 8,
            },
            AtomSymbol::Fr => ElementData {
                atomic_number: 87,
                standard_mass: 223.0,
                valence_electrons: 1,
            },
            AtomSymbol::Ra => ElementData {
                atomic_number: 88,
                standard_mass: 226.0,
                valence_electrons: 2,
            },
            AtomSymbol::Ac => ElementData {
                atomic_number: 89,
                standard_mass: 227.0,
                valence_electrons: 3,
            },
            AtomSymbol::Th => ElementData {
                atomic_number: 90,
                standard_mass: 232.04,
                valence_electrons: 4,
            },
            AtomSymbol::Pa => ElementData {
                atomic_number: 91,
                standard_mass: 231.04,
                valence_electrons: 5,
            },
            AtomSymbol::U => ElementData {
                atomic_number: 92,
                standard_mass: 238.03,
                valence_electrons: 6,
            },
            AtomSymbol::Np => ElementData {
                atomic_number: 93,
                standard_mass: 237.0,
                valence_electrons: 7,
            },
            AtomSymbol::Pu => ElementData {
                atomic_number: 94,
                standard_mass: 244.0,
                valence_electrons: 8,
            },
            AtomSymbol::Am => ElementData {
                atomic_number: 95,
                standard_mass: 243.0,
                valence_electrons: 9,
            },
            AtomSymbol::Cm => ElementData {
                atomic_number: 96,
                standard_mass: 247.0,
                valence_electrons: 10,
            },
            AtomSymbol::Bk => ElementData {
                atomic_number: 97,
                standard_mass: 247.0,
                valence_electrons: 11,
            },
            AtomSymbol::Cf => ElementData {
                atomic_number: 98,
                standard_mass: 251.0,
                valence_electrons: 12,
            },
            AtomSymbol::Es => ElementData {
                atomic_number: 99,
                standard_mass: 252.0,
                valence_electrons: 13,
            },
            AtomSymbol::Fm => ElementData {
                atomic_number: 100,
                standard_mass: 257.0,
                valence_electrons: 14,
            },
            AtomSymbol::Md => ElementData {
                atomic_number: 101,
                standard_mass: 258.0,
                valence_electrons: 15,
            },
            AtomSymbol::No => ElementData {
                atomic_number: 102,
                standard_mass: 259.0,
                valence_electrons: 16,
            },
            AtomSymbol::Lr => ElementData {
                atomic_number: 103,
                standard_mass: 266.0,
                valence_electrons: 3,
            },
            AtomSymbol::Rf => ElementData {
                atomic_number: 104,
                standard_mass: 267.0,
                valence_electrons: 4,
            },
            AtomSymbol::Db => ElementData {
                atomic_number: 105,
                standard_mass: 268.0,
                valence_electrons: 5,
            },
            AtomSymbol::Sg => ElementData {
                atomic_number: 106,
                standard_mass: 269.0,
                valence_electrons: 6,
            },
            AtomSymbol::Bh => ElementData {
                atomic_number: 107,
                standard_mass: 270.0,
                valence_electrons: 7,
            },
            AtomSymbol::Hs => ElementData {
                atomic_number: 108,
                standard_mass: 277.0,
                valence_electrons: 8,
            },
            AtomSymbol::Mt => ElementData {
                atomic_number: 109,
                standard_mass: 278.0,
                valence_electrons: 9,
            },
            AtomSymbol::Ds => ElementData {
                atomic_number: 110,
                standard_mass: 281.0,
                valence_electrons: 10,
            },
            AtomSymbol::Rg => ElementData {
                atomic_number: 111,
                standard_mass: 282.0,
                valence_electrons: 11,
            },
            AtomSymbol::Cn => ElementData {
                atomic_number: 112,
                standard_mass: 285.0,
                valence_electrons: 2,
            },
            AtomSymbol::Nh => ElementData {
                atomic_number: 113,
                standard_mass: 286.0,
                valence_electrons: 3,
            },
            AtomSymbol::Fl => ElementData {
                atomic_number: 114,
                standard_mass: 289.0,
                valence_electrons: 4,
            },
            AtomSymbol::Mc => ElementData {
                atomic_number: 115,
                standard_mass: 290.0,
                valence_electrons: 5,
            },
            AtomSymbol::Lv => ElementData {
                atomic_number: 116,
                standard_mass: 293.0,
                valence_electrons: 6,
            },
            AtomSymbol::Ts => ElementData {
                atomic_number: 117,
                standard_mass: 294.0,
                valence_electrons: 7,
            },
            AtomSymbol::Og => ElementData {
                atomic_number: 118,
                standard_mass: 294.0,
                valence_electrons: 8,
            },
            AtomSymbol::Wildcard => ElementData {
                atomic_number: 0,
                standard_mass: 0.0,
                valence_electrons: 0,
            },
            AtomSymbol::Organic(o) => o.element_data(),
        }
    }

    /// Returns the atomic number (Z) of this element.
    pub const fn atomic_number(&self) -> u8 {
        self.element_data().atomic_number
    }

    /// Returns the standard atomic mass in Daltons.
    pub const fn standard_mass(&self) -> f64 {
        self.element_data().standard_mass
    }

    /// Returns the number of valence electrons for main-group elements.
    pub const fn valence_electrons(&self) -> u8 {
        self.element_data().valence_electrons
    }
}

impl OrganicAtom {
    /// Returns the static chemical data for this organic atom.
    /// Evaluated at compile time when possible.
    pub const fn element_data(&self) -> ElementData {
        match self {
            OrganicAtom::B => ElementData {
                atomic_number: 5,
                standard_mass: 10.81,
                valence_electrons: 3,
            },
            OrganicAtom::C => ElementData {
                atomic_number: 6,
                standard_mass: 12.011,
                valence_electrons: 4,
            },
            OrganicAtom::N => ElementData {
                atomic_number: 7,
                standard_mass: 14.007,
                valence_electrons: 5,
            },
            OrganicAtom::O => ElementData {
                atomic_number: 8,
                standard_mass: 15.999,
                valence_electrons: 6,
            },
            OrganicAtom::F => ElementData {
                atomic_number: 9,
                standard_mass: 18.998,
                valence_electrons: 7,
            },
            OrganicAtom::P => ElementData {
                atomic_number: 15,
                standard_mass: 30.974,
                valence_electrons: 5,
            },
            OrganicAtom::S => ElementData {
                atomic_number: 16,
                standard_mass: 32.06,
                valence_electrons: 6,
            },
            OrganicAtom::Cl => ElementData {
                atomic_number: 17,
                standard_mass: 35.45,
                valence_electrons: 7,
            },
            OrganicAtom::Br => ElementData {
                atomic_number: 35,
                standard_mass: 79.904,
                valence_electrons: 7,
            },
            OrganicAtom::I => ElementData {
                atomic_number: 53,
                standard_mass: 126.90,
                valence_electrons: 7,
            },
        }
    }
}

/// Returns the mass of a specific isotope.
/// For common isotopes, returns the known mass; otherwise, uses the mass number as approximation.
pub const fn isotope_mass(element: &AtomSymbol, mass_number: u16) -> f64 {
    match (element.atomic_number(), mass_number) {
        // Hydrogen isotopes
        (1, 1) => 1.00783, // protium
        (1, 2) => 2.01410, // deuterium
        (1, 3) => 3.01605, // tritium
        // Carbon isotopes
        (6, 12) => 12.0,     // carbon-12 (exact by definition)
        (6, 13) => 13.00335, // carbon-13
        (6, 14) => 14.00324, // carbon-14
        // Nitrogen isotopes
        (7, 14) => 14.00307,
        (7, 15) => 15.00011,
        // Oxygen isotopes
        (8, 16) => 15.99491,
        (8, 17) => 16.99913,
        (8, 18) => 17.99916,
        // Fallback: mass number as approximation
        (_, m) => m as f64,
    }
}

impl super::atom::Atom {
    /// Returns the mass: isotopic if specified, standard otherwise.
    pub fn mass(&self) -> f64 {
        match self.isotope() {
            Some(mass_number) => isotope_mass(self.element(), mass_number),
            None => self.element().standard_mass(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn organic_atoms_have_correct_data() {
        assert_eq!(AtomSymbol::Organic(OrganicAtom::B).atomic_number(), 5);
        assert_eq!(AtomSymbol::Organic(OrganicAtom::C).atomic_number(), 6);
        assert_eq!(AtomSymbol::Organic(OrganicAtom::N).atomic_number(), 7);
        assert_eq!(AtomSymbol::Organic(OrganicAtom::O).atomic_number(), 8);
        assert_eq!(AtomSymbol::Organic(OrganicAtom::F).atomic_number(), 9);
        assert_eq!(AtomSymbol::Organic(OrganicAtom::P).atomic_number(), 15);
        assert_eq!(AtomSymbol::Organic(OrganicAtom::S).atomic_number(), 16);
        assert_eq!(AtomSymbol::Organic(OrganicAtom::Cl).atomic_number(), 17);
        assert_eq!(AtomSymbol::Organic(OrganicAtom::Br).atomic_number(), 35);
        assert_eq!(AtomSymbol::Organic(OrganicAtom::I).atomic_number(), 53);
    }

    #[test]
    fn aromatic_elements_have_correct_valence_electrons() {
        assert_eq!(AtomSymbol::Organic(OrganicAtom::B).valence_electrons(), 3);
        assert_eq!(AtomSymbol::Organic(OrganicAtom::C).valence_electrons(), 4);
        assert_eq!(AtomSymbol::Organic(OrganicAtom::N).valence_electrons(), 5);
        assert_eq!(AtomSymbol::Organic(OrganicAtom::O).valence_electrons(), 6);
        assert_eq!(AtomSymbol::Organic(OrganicAtom::P).valence_electrons(), 5);
        assert_eq!(AtomSymbol::Organic(OrganicAtom::S).valence_electrons(), 6);
        assert_eq!(AtomSymbol::As.valence_electrons(), 5);
        assert_eq!(AtomSymbol::Se.valence_electrons(), 6);
        assert_eq!(AtomSymbol::Te.valence_electrons(), 6);
    }

    #[test]
    fn wildcard_returns_zero() {
        let data = AtomSymbol::Wildcard.element_data();
        assert_eq!(data.atomic_number, 0);
        assert_eq!(data.standard_mass, 0.0);
        assert_eq!(data.valence_electrons, 0);
    }

    #[test]
    fn atomic_numbers_are_unique() {
        let all_elements: Vec<AtomSymbol> = vec![
            AtomSymbol::H,
            AtomSymbol::He,
            AtomSymbol::Li,
            AtomSymbol::Be,
            AtomSymbol::Ne,
            AtomSymbol::Na,
            AtomSymbol::Mg,
            AtomSymbol::Al,
            AtomSymbol::Si,
            AtomSymbol::Ar,
            AtomSymbol::K,
            AtomSymbol::Ca,
            AtomSymbol::Sc,
            AtomSymbol::Ti,
            AtomSymbol::V,
            AtomSymbol::Cr,
            AtomSymbol::Mn,
            AtomSymbol::Fe,
            AtomSymbol::Co,
            AtomSymbol::Ni,
            AtomSymbol::Cu,
            AtomSymbol::Zn,
            AtomSymbol::Ga,
            AtomSymbol::Ge,
            AtomSymbol::As,
            AtomSymbol::Se,
            AtomSymbol::Kr,
            AtomSymbol::Organic(OrganicAtom::B),
            AtomSymbol::Organic(OrganicAtom::C),
            AtomSymbol::Organic(OrganicAtom::N),
            AtomSymbol::Organic(OrganicAtom::O),
            AtomSymbol::Organic(OrganicAtom::F),
            AtomSymbol::Organic(OrganicAtom::P),
            AtomSymbol::Organic(OrganicAtom::S),
            AtomSymbol::Organic(OrganicAtom::Cl),
            AtomSymbol::Organic(OrganicAtom::Br),
            AtomSymbol::Organic(OrganicAtom::I),
        ];

        let mut seen = HashSet::new();
        for elem in all_elements {
            let z = elem.atomic_number();
            assert!(
                seen.insert(z),
                "Duplicate atomic number {} for {:?}",
                z,
                elem
            );
        }
    }

    #[test]
    fn standard_mass_is_positive() {
        let elements = vec![
            AtomSymbol::Organic(OrganicAtom::C),
            AtomSymbol::Organic(OrganicAtom::N),
            AtomSymbol::Organic(OrganicAtom::O),
            AtomSymbol::H,
            AtomSymbol::Fe,
        ];
        for elem in elements {
            assert!(
                elem.standard_mass() > 0.0,
                "{:?} should have positive mass",
                elem
            );
        }
    }

    #[test]
    fn isotope_mass_known_isotopes() {
        let c = AtomSymbol::Organic(OrganicAtom::C);
        assert!((isotope_mass(&c, 12) - 12.0).abs() < 1e-5);
        assert!((isotope_mass(&c, 13) - 13.00335).abs() < 1e-4);

        let h = AtomSymbol::H;
        assert!((isotope_mass(&h, 2) - 2.01410).abs() < 1e-4);
    }

    #[test]
    fn isotope_mass_unknown_fallback() {
        let c = AtomSymbol::Organic(OrganicAtom::C);
        // Unknown isotope: should return mass number as float
        assert!((isotope_mass(&c, 99) - 99.0).abs() < 1e-5);
    }

    #[test]
    fn atom_mass_with_and_without_isotope() {
        use super::super::atom::Atom;
        let c = Atom::new(AtomSymbol::Organic(OrganicAtom::C), 0, None).unwrap();
        assert!((c.mass() - 12.011).abs() < 1e-3);

        let c13 = Atom::new(AtomSymbol::Organic(OrganicAtom::C), 0, Some(13)).unwrap();
        assert!((c13.mass() - 13.00335).abs() < 1e-4);
    }
}
