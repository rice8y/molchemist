use std::fmt;
use std::str::FromStr;

use crate::AtomError;

/// An element symbol, covering the full periodic table plus the organic subset and wildcard.
///
/// Atoms parsed outside brackets (e.g. `C`, `N`, `Cl`) are represented as
/// `Organic(OrganicAtom)` variants, which carry implicit-hydrogen rules.
/// Atoms inside brackets (e.g. `[Fe]`, `[13C]`) use the element variants directly.
///
/// The wildcard atom `*` is represented by `Wildcard`.
#[derive(Debug, Clone, PartialEq, Copy)]
#[rustfmt::skip]
pub enum AtomSymbol {
    H, He,
    Li, Be, Ne,
    Na, Mg, Al, Si, Ar,
    K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Kr,
    Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, Xe,
    Cs, Ba, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn,
    Fr, Ra, Lr, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Nh, Fl, Mc, Lv, Ts, Og,
    La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb,
    Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No,
    /// Wildcard atom (`*`).
    Wildcard,
    /// An atom from the organic subset (B, C, N, O, P, S, F, Cl, Br, I).
    Organic(OrganicAtom)
}

impl AtomSymbol {
    /// Returns `true` if this element can participate in an aromatic ring.
    pub fn can_be_aromatic(self) -> bool {
        match self {
            AtomSymbol::Organic(OrganicAtom::C) => true,
            AtomSymbol::Organic(OrganicAtom::N) => true,
            AtomSymbol::Organic(OrganicAtom::O) => true,
            AtomSymbol::Organic(OrganicAtom::S) => true,
            AtomSymbol::Organic(OrganicAtom::P) => true,
            AtomSymbol::Se => true,
            AtomSymbol::As => true,
            AtomSymbol::Organic(OrganicAtom::B) => true,
            AtomSymbol::Te => true,
            // Wildcard can be part of an aromatic ring per OpenSMILES spec
            AtomSymbol::Wildcard => true,
            _ => false,
        }
    }
}

impl fmt::Display for AtomSymbol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AtomSymbol::H => write!(f, "H"),
            AtomSymbol::He => write!(f, "He"),
            AtomSymbol::Li => write!(f, "Li"),
            AtomSymbol::Be => write!(f, "Be"),
            AtomSymbol::Ne => write!(f, "Ne"),
            AtomSymbol::Na => write!(f, "Na"),
            AtomSymbol::Mg => write!(f, "Mg"),
            AtomSymbol::Al => write!(f, "Al"),
            AtomSymbol::Si => write!(f, "Si"),
            AtomSymbol::Ar => write!(f, "Ar"),
            AtomSymbol::K => write!(f, "K"),
            AtomSymbol::Ca => write!(f, "Ca"),
            AtomSymbol::Sc => write!(f, "Sc"),
            AtomSymbol::Ti => write!(f, "Ti"),
            AtomSymbol::V => write!(f, "V"),
            AtomSymbol::Cr => write!(f, "Cr"),
            AtomSymbol::Mn => write!(f, "Mn"),
            AtomSymbol::Fe => write!(f, "Fe"),
            AtomSymbol::Co => write!(f, "Co"),
            AtomSymbol::Ni => write!(f, "Ni"),
            AtomSymbol::Cu => write!(f, "Cu"),
            AtomSymbol::Zn => write!(f, "Zn"),
            AtomSymbol::Ga => write!(f, "Ga"),
            AtomSymbol::Ge => write!(f, "Ge"),
            AtomSymbol::As => write!(f, "As"),
            AtomSymbol::Se => write!(f, "Se"),
            AtomSymbol::Kr => write!(f, "Kr"),
            AtomSymbol::Rb => write!(f, "Rb"),
            AtomSymbol::Sr => write!(f, "Sr"),
            AtomSymbol::Y => write!(f, "Y"),
            AtomSymbol::Zr => write!(f, "Zr"),
            AtomSymbol::Nb => write!(f, "Nb"),
            AtomSymbol::Mo => write!(f, "Mo"),
            AtomSymbol::Tc => write!(f, "Tc"),
            AtomSymbol::Ru => write!(f, "Ru"),
            AtomSymbol::Rh => write!(f, "Rh"),
            AtomSymbol::Pd => write!(f, "Pd"),
            AtomSymbol::Ag => write!(f, "Ag"),
            AtomSymbol::Cd => write!(f, "Cd"),
            AtomSymbol::In => write!(f, "In"),
            AtomSymbol::Sn => write!(f, "Sn"),
            AtomSymbol::Sb => write!(f, "Sb"),
            AtomSymbol::Te => write!(f, "Te"),
            AtomSymbol::Xe => write!(f, "Xe"),
            AtomSymbol::Cs => write!(f, "Cs"),
            AtomSymbol::Ba => write!(f, "Ba"),
            AtomSymbol::La => write!(f, "La"),
            AtomSymbol::Ce => write!(f, "Ce"),
            AtomSymbol::Pr => write!(f, "Pr"),
            AtomSymbol::Nd => write!(f, "Nd"),
            AtomSymbol::Pm => write!(f, "Pm"),
            AtomSymbol::Sm => write!(f, "Sm"),
            AtomSymbol::Eu => write!(f, "Eu"),
            AtomSymbol::Gd => write!(f, "Gd"),
            AtomSymbol::Tb => write!(f, "Tb"),
            AtomSymbol::Dy => write!(f, "Dy"),
            AtomSymbol::Ho => write!(f, "Ho"),
            AtomSymbol::Er => write!(f, "Er"),
            AtomSymbol::Tm => write!(f, "Tm"),
            AtomSymbol::Yb => write!(f, "Yb"),
            AtomSymbol::Lu => write!(f, "Lu"),
            AtomSymbol::Hf => write!(f, "Hf"),
            AtomSymbol::Ta => write!(f, "Ta"),
            AtomSymbol::W => write!(f, "W"),
            AtomSymbol::Re => write!(f, "Re"),
            AtomSymbol::Os => write!(f, "Os"),
            AtomSymbol::Ir => write!(f, "Ir"),
            AtomSymbol::Pt => write!(f, "Pt"),
            AtomSymbol::Au => write!(f, "Au"),
            AtomSymbol::Hg => write!(f, "Hg"),
            AtomSymbol::Tl => write!(f, "Tl"),
            AtomSymbol::Pb => write!(f, "Pb"),
            AtomSymbol::Bi => write!(f, "Bi"),
            AtomSymbol::Po => write!(f, "Po"),
            AtomSymbol::At => write!(f, "At"),
            AtomSymbol::Rn => write!(f, "Rn"),
            AtomSymbol::Fr => write!(f, "Fr"),
            AtomSymbol::Ra => write!(f, "Ra"),
            AtomSymbol::Ac => write!(f, "Ac"),
            AtomSymbol::Th => write!(f, "Th"),
            AtomSymbol::Pa => write!(f, "Pa"),
            AtomSymbol::U => write!(f, "U"),
            AtomSymbol::Np => write!(f, "Np"),
            AtomSymbol::Pu => write!(f, "Pu"),
            AtomSymbol::Am => write!(f, "Am"),
            AtomSymbol::Cm => write!(f, "Cm"),
            AtomSymbol::Bk => write!(f, "Bk"),
            AtomSymbol::Cf => write!(f, "Cf"),
            AtomSymbol::Es => write!(f, "Es"),
            AtomSymbol::Fm => write!(f, "Fm"),
            AtomSymbol::Md => write!(f, "Md"),
            AtomSymbol::No => write!(f, "No"),
            AtomSymbol::Lr => write!(f, "Lr"),
            AtomSymbol::Rf => write!(f, "Rf"),
            AtomSymbol::Db => write!(f, "Db"),
            AtomSymbol::Sg => write!(f, "Sg"),
            AtomSymbol::Bh => write!(f, "Bh"),
            AtomSymbol::Hs => write!(f, "Hs"),
            AtomSymbol::Mt => write!(f, "Mt"),
            AtomSymbol::Ds => write!(f, "Ds"),
            AtomSymbol::Rg => write!(f, "Rg"),
            AtomSymbol::Cn => write!(f, "Cn"),
            AtomSymbol::Nh => write!(f, "Nh"),
            AtomSymbol::Fl => write!(f, "Fl"),
            AtomSymbol::Mc => write!(f, "Mc"),
            AtomSymbol::Lv => write!(f, "Lv"),
            AtomSymbol::Ts => write!(f, "Ts"),
            AtomSymbol::Og => write!(f, "Og"),
            AtomSymbol::Wildcard => write!(f, "*"),
            AtomSymbol::Organic(organic) => write!(f, "{}", organic),
        }
    }
}

impl FromStr for AtomSymbol {
    type Err = AtomError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let normalized = s.to_uppercase();
        match normalized.as_str() {
            "H" => Ok(AtomSymbol::H),
            "HE" => Ok(AtomSymbol::He),
            "LI" => Ok(AtomSymbol::Li),
            "BE" => Ok(AtomSymbol::Be),
            "B" => Ok(AtomSymbol::Organic(OrganicAtom::B)),
            "C" => Ok(AtomSymbol::Organic(OrganicAtom::C)),
            "N" => Ok(AtomSymbol::Organic(OrganicAtom::N)),
            "O" => Ok(AtomSymbol::Organic(OrganicAtom::O)),
            "F" => Ok(AtomSymbol::Organic(OrganicAtom::F)),
            "NE" => Ok(AtomSymbol::Ne),
            "NA" => Ok(AtomSymbol::Na),
            "MG" => Ok(AtomSymbol::Mg),
            "AL" => Ok(AtomSymbol::Al),
            "SI" => Ok(AtomSymbol::Si),
            "P" => Ok(AtomSymbol::Organic(OrganicAtom::P)),
            "S" => Ok(AtomSymbol::Organic(OrganicAtom::S)),
            "CL" => Ok(AtomSymbol::Organic(OrganicAtom::Cl)),
            "AR" => Ok(AtomSymbol::Ar),
            "K" => Ok(AtomSymbol::K),
            "CA" => Ok(AtomSymbol::Ca),
            "SC" => Ok(AtomSymbol::Sc),
            "TI" => Ok(AtomSymbol::Ti),
            "V" => Ok(AtomSymbol::V),
            "CR" => Ok(AtomSymbol::Cr),
            "MN" => Ok(AtomSymbol::Mn),
            "FE" => Ok(AtomSymbol::Fe),
            "CO" => Ok(AtomSymbol::Co),
            "NI" => Ok(AtomSymbol::Ni),
            "CU" => Ok(AtomSymbol::Cu),
            "ZN" => Ok(AtomSymbol::Zn),
            "GA" => Ok(AtomSymbol::Ga),
            "GE" => Ok(AtomSymbol::Ge),
            "AS" => Ok(AtomSymbol::As),
            "SE" => Ok(AtomSymbol::Se),
            "BR" => Ok(AtomSymbol::Organic(OrganicAtom::Br)),
            "KR" => Ok(AtomSymbol::Kr),
            "RB" => Ok(AtomSymbol::Rb),
            "SR" => Ok(AtomSymbol::Sr),
            "Y" => Ok(AtomSymbol::Y),
            "ZR" => Ok(AtomSymbol::Zr),
            "NB" => Ok(AtomSymbol::Nb),
            "MO" => Ok(AtomSymbol::Mo),
            "TC" => Ok(AtomSymbol::Tc),
            "RU" => Ok(AtomSymbol::Ru),
            "RH" => Ok(AtomSymbol::Rh),
            "PD" => Ok(AtomSymbol::Pd),
            "AG" => Ok(AtomSymbol::Ag),
            "CD" => Ok(AtomSymbol::Cd),
            "IN" => Ok(AtomSymbol::In),
            "SN" => Ok(AtomSymbol::Sn),
            "SB" => Ok(AtomSymbol::Sb),
            "TE" => Ok(AtomSymbol::Te),
            "I" => Ok(AtomSymbol::Organic(OrganicAtom::I)),
            "XE" => Ok(AtomSymbol::Xe),
            "CS" => Ok(AtomSymbol::Cs),
            "BA" => Ok(AtomSymbol::Ba),
            "LA" => Ok(AtomSymbol::La),
            "CE" => Ok(AtomSymbol::Ce),
            "PR" => Ok(AtomSymbol::Pr),
            "ND" => Ok(AtomSymbol::Nd),
            "PM" => Ok(AtomSymbol::Pm),
            "SM" => Ok(AtomSymbol::Sm),
            "EU" => Ok(AtomSymbol::Eu),
            "GD" => Ok(AtomSymbol::Gd),
            "TB" => Ok(AtomSymbol::Tb),
            "DY" => Ok(AtomSymbol::Dy),
            "HO" => Ok(AtomSymbol::Ho),
            "ER" => Ok(AtomSymbol::Er),
            "TM" => Ok(AtomSymbol::Tm),
            "YB" => Ok(AtomSymbol::Yb),
            "LU" => Ok(AtomSymbol::Lu),
            "HF" => Ok(AtomSymbol::Hf),
            "TA" => Ok(AtomSymbol::Ta),
            "W" => Ok(AtomSymbol::W),
            "RE" => Ok(AtomSymbol::Re),
            "OS" => Ok(AtomSymbol::Os),
            "IR" => Ok(AtomSymbol::Ir),
            "PT" => Ok(AtomSymbol::Pt),
            "AU" => Ok(AtomSymbol::Au),
            "HG" => Ok(AtomSymbol::Hg),
            "TL" => Ok(AtomSymbol::Tl),
            "PB" => Ok(AtomSymbol::Pb),
            "BI" => Ok(AtomSymbol::Bi),
            "PO" => Ok(AtomSymbol::Po),
            "AT" => Ok(AtomSymbol::At),
            "RN" => Ok(AtomSymbol::Rn),
            "FR" => Ok(AtomSymbol::Fr),
            "RA" => Ok(AtomSymbol::Ra),
            "AC" => Ok(AtomSymbol::Ac),
            "TH" => Ok(AtomSymbol::Th),
            "PA" => Ok(AtomSymbol::Pa),
            "U" => Ok(AtomSymbol::U),
            "NP" => Ok(AtomSymbol::Np),
            "PU" => Ok(AtomSymbol::Pu),
            "AM" => Ok(AtomSymbol::Am),
            "CM" => Ok(AtomSymbol::Cm),
            "BK" => Ok(AtomSymbol::Bk),
            "CF" => Ok(AtomSymbol::Cf),
            "ES" => Ok(AtomSymbol::Es),
            "FM" => Ok(AtomSymbol::Fm),
            "MD" => Ok(AtomSymbol::Md),
            "NO" => Ok(AtomSymbol::No),
            "LR" => Ok(AtomSymbol::Lr),
            "RF" => Ok(AtomSymbol::Rf),
            "DB" => Ok(AtomSymbol::Db),
            "SG" => Ok(AtomSymbol::Sg),
            "BH" => Ok(AtomSymbol::Bh),
            "HS" => Ok(AtomSymbol::Hs),
            "MT" => Ok(AtomSymbol::Mt),
            "DS" => Ok(AtomSymbol::Ds),
            "RG" => Ok(AtomSymbol::Rg),
            "CN" => Ok(AtomSymbol::Cn),
            "NH" => Ok(AtomSymbol::Nh),
            "FL" => Ok(AtomSymbol::Fl),
            "MC" => Ok(AtomSymbol::Mc),
            "LV" => Ok(AtomSymbol::Lv),
            "TS" => Ok(AtomSymbol::Ts),
            "OG" => Ok(AtomSymbol::Og),
            "*" => Ok(AtomSymbol::Wildcard),
            _ => Err(AtomError::UnknownElement(s.to_string())),
        }
    }
}

/// An atom from the SMILES organic subset.
///
/// These atoms can appear outside brackets and have well-defined implicit
/// hydrogen rules based on their standard valences.
#[derive(Debug, Clone, PartialEq, Copy)]
pub enum OrganicAtom {
    B,
    C,
    N,
    O,
    P,
    S,
    F,
    Cl,
    Br,
    I,
}

impl OrganicAtom {
    /// Returns the standard valences for this organic atom.
    ///
    /// The lowest valence ≥ current bond order sum determines implicit hydrogen count.
    pub fn valence(&self) -> &'static [u8] {
        match self {
            OrganicAtom::B => &[3],
            OrganicAtom::C => &[4],
            OrganicAtom::N => &[3, 5],
            OrganicAtom::O => &[2],
            OrganicAtom::P => &[3, 5],
            OrganicAtom::S => &[2, 4, 6],
            OrganicAtom::F => &[1],
            OrganicAtom::Cl => &[1],
            OrganicAtom::Br => &[1],
            OrganicAtom::I => &[1],
        }
    }

    /// Calculates the number of implicit hydrogens for this atom.
    ///
    /// `bond_order_sum` is the sum of bond orders of all bonds to this atom.
    /// `aromatic` applies the OpenSMILES aromatic subvalence rule:
    /// if subvalence > 1, implicit H = subvalence − 1; otherwise 0.
    pub fn implicit_hydrogens(&self, bond_order_sum: u8, aromatic: bool) -> u8 {
        for v in self.valence() {
            if *v >= bond_order_sum {
                let subvalence = *v - bond_order_sum;
                if aromatic {
                    if subvalence > 1 {
                        return subvalence - 1;
                    } else {
                        return 0;
                    }
                }
                return subvalence;
            }
        }

        0
    }
}

impl fmt::Display for OrganicAtom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OrganicAtom::B => write!(f, "B"),
            OrganicAtom::C => write!(f, "C"),
            OrganicAtom::N => write!(f, "N"),
            OrganicAtom::O => write!(f, "O"),
            OrganicAtom::P => write!(f, "P"),
            OrganicAtom::S => write!(f, "S"),
            OrganicAtom::F => write!(f, "F"),
            OrganicAtom::Cl => write!(f, "Cl"),
            OrganicAtom::Br => write!(f, "Br"),
            OrganicAtom::I => write!(f, "I"),
        }
    }
}

/// A chemical atom with element, formal charge, and isotope.
#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    element: AtomSymbol,
    charge: i8,           // -15 to +15
    isotope: Option<u16>, // 0 to 999 or unspecified
}

impl Atom {
    pub fn new(element: AtomSymbol, charge: i8, isotope: Option<u16>) -> Result<Atom, AtomError> {
        if !(-15..=15).contains(&charge) {
            return Err(AtomError::InvalidCharge(charge));
        }

        if let Some(value) = isotope {
            if value > 999 {
                return Err(AtomError::InvalidIsotope(value));
            }
        }

        Ok(Atom {
            element,
            charge,
            isotope,
        })
    }

    /// Calculates implicit hydrogens for organic-subset atoms.
    ///
    /// Returns `Ok(0)` for bracket atoms (non-organic), where hydrogen count
    /// is always explicit. Returns an error if `bond_order_sum` is missing
    /// for an organic atom.
    pub fn implicit_hydrogens(
        &self,
        bond_order_sum: Option<u8>,
        aromatic: bool,
    ) -> Result<u8, AtomError> {
        if let AtomSymbol::Organic(organic) = self.element() {
            Ok(organic
                .implicit_hydrogens(bond_order_sum.ok_or(AtomError::MissingBondOrder)?, aromatic))
        } else {
            Ok(0)
        }
    }

    /// Returns `true` if this atom is from the organic subset (B, C, N, O, P, S, F, Cl, Br, I).
    pub fn is_organic(&self) -> bool {
        matches!(self.element, AtomSymbol::Organic(_))
    }

    /// Returns the element symbol.
    pub fn element(&self) -> &AtomSymbol {
        &self.element
    }

    /// Returns the formal charge (range: −15 to +15).
    pub fn charge(&self) -> i8 {
        self.charge
    }

    /// Returns the isotope mass number, if specified (e.g. `[13C]` → `Some(13)`).
    pub fn isotope(&self) -> Option<u16> {
        self.isotope
    }
}
