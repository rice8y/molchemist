use crate::{BondError, NodeIndex};

/// The type of a covalent bond between two atoms.
///
/// Corresponds directly to the bond symbols in SMILES notation:
///
/// | Variant       | Symbol | Description                            |
/// |---------------|--------|----------------------------------------|
/// | `Simple`      | `-`    | Single bond (default when omitted)     |
/// | `Double`      | `=`    | Double bond                            |
/// | `Triple`      | `#`    | Triple bond                            |
/// | `Quadruple`   | `$`    | Quadruple bond                         |
/// | `Aromatic`    | `:`    | Aromatic bond                          |
/// | `Up`          | `/`    | Directional bond (E/Z stereochemistry) |
/// | `Down`        | `\`    | Directional bond (E/Z stereochemistry) |
/// | `Disconnected`| `.`    | No bond (disconnected fragments)       |
#[derive(Debug, Clone, PartialEq, Copy)]
pub enum BondType {
    Simple,
    Double,
    Triple,
    Quadruple,
    Aromatic,
    Disconnected,
    Up,
    Down,
}

impl TryFrom<&char> for BondType {
    type Error = BondError;

    fn try_from(c: &char) -> Result<Self, Self::Error> {
        match c {
            '-' => Ok(BondType::Simple),
            '=' => Ok(BondType::Double),
            '#' => Ok(BondType::Triple),
            '$' => Ok(BondType::Quadruple),
            '.' => Ok(BondType::Disconnected),
            ':' => Ok(BondType::Aromatic),
            '/' => Ok(BondType::Up),
            '\\' => Ok(BondType::Down),
            _ => Err(BondError::UnknownBond(*c)),
        }
    }
}

impl BondType {
    /// Returns the number of electrons involved in this bond.
    pub fn electrons_involved(&self) -> u8 {
        match self {
            BondType::Simple => 2,
            BondType::Double => 4,
            BondType::Triple => 6,
            BondType::Quadruple => 8,
            BondType::Aromatic => 3,
            BondType::Disconnected => 0,
            BondType::Up => 2,
            BondType::Down => 2,
        }
    }

    /// Returns a priority used when sorting neighbors in the spanning tree.
    ///
    /// Higher-order bonds get higher priority so they become tree (chain) edges
    /// rather than ring-closure back edges, keeping ring closures simple in output.
    pub(crate) fn bond_order_priority(&self) -> u8 {
        match self {
            BondType::Quadruple => 4,
            BondType::Triple => 3,
            BondType::Double => 2,
            BondType::Aromatic => 1,
            _ => 0,
        }
    }

    /// Returns the bond order contribution × 2 for implicit hydrogen calculation.
    ///
    /// Per the OpenSMILES spec, aromatic bonds count as 1.0 (not 1.5)
    /// for the purpose of calculating implicit hydrogens. The value is
    /// multiplied by 2 to avoid floating-point arithmetic.
    pub(crate) fn bond_order_x2_for_implicit_h(&self) -> u8 {
        match self {
            BondType::Simple => 2,
            BondType::Double => 4,
            BondType::Triple => 6,
            BondType::Quadruple => 8,
            BondType::Aromatic => 2,
            BondType::Disconnected => 0,
            BondType::Up => 2,
            BondType::Down => 2,
        }
    }
}

/// A bond connecting two atom nodes in a [`Molecule`](crate::Molecule).
///
/// Bonds are undirected but stored with a `source` and `target` index
/// (both indices into [`Molecule::nodes()`](crate::Molecule::nodes)).
#[derive(Debug, Clone, PartialEq)]
pub struct Bond {
    kind: BondType,
    source: NodeIndex,
    target: NodeIndex,
}

impl Bond {
    pub(crate) fn new(kind: BondType, source: NodeIndex, target: NodeIndex) -> Bond {
        Bond {
            kind,
            source,
            target,
        }
    }

    /// Returns the bond type.
    pub fn kind(&self) -> BondType {
        self.kind
    }

    /// Returns the index of the source node.
    pub fn source(&self) -> NodeIndex {
        self.source
    }

    /// Returns the index of the target node.
    pub fn target(&self) -> NodeIndex {
        self.target
    }
}
