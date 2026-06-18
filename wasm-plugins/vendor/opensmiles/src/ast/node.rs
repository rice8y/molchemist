use crate::{
    ast::atom::{Atom, AtomSymbol},
    ast::chirality::Chirality,
    NodeError,
};

/// A fully resolved atom node in the molecule graph.
///
/// Each `Node` holds an [`Atom`] (element, charge, isotope) along with
/// aromaticity, explicit hydrogen count, atom class, and chirality — all
/// resolved at parse time.
///
/// Nodes are connected to each other via [`Bond`](crate::Bond)s in a
/// [`Molecule`](crate::Molecule).
#[derive(Debug, Clone, PartialEq)]
pub struct Node {
    atom: Atom,
    aromatic: bool,
    hydrogens: u8,
    class: Option<u16>,
    chirality: Option<Chirality>,
}

impl Node {
    pub fn new(
        atom: Atom,
        aromatic: bool,
        hydrogens: u8,
        class: Option<u16>,
        chirality: Option<Chirality>,
    ) -> Result<Node, NodeError> {
        if hydrogens > 9 {
            return Err(NodeError::InvalidHydrogen(hydrogens));
        }

        match class {
            None => (),
            Some(value) => {
                if value > 999 {
                    return Err(NodeError::InvalidClass(value));
                }
            }
        }

        if aromatic && !atom.element().can_be_aromatic() {
            return Err(NodeError::InvalidAromaticElement(*atom.element()));
        }

        Ok(Node {
            atom,
            aromatic,
            hydrogens,
            class,
            chirality,
        })
    }

    /// Returns the atom (element, charge, isotope).
    pub fn atom(&self) -> &Atom {
        &self.atom
    }

    /// Returns `true` if this atom is aromatic (written in lowercase in SMILES).
    pub fn aromatic(&self) -> bool {
        self.aromatic
    }

    /// Returns the number of hydrogens attached to this atom.
    ///
    /// For bracket atoms (e.g. `[CH3]`) this is the explicit hcount.
    /// For organic-subset atoms (e.g. `C`) this is the implicit hydrogen
    /// count calculated from valence and bond order sum.
    pub fn hydrogens(&self) -> u8 {
        self.hydrogens
    }

    /// Returns the atom class, if specified (e.g. `[C:1]` → `Some(1)`).
    pub fn class(&self) -> Option<u16> {
        self.class
    }

    /// Returns the chirality specification, if any.
    pub fn chirality(&self) -> Option<Chirality> {
        self.chirality
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct NodeBuilder {
    atom: Atom,
    aromatic: Option<bool>,
    hydrogens: Option<u8>,
    class: Option<u16>,
    chirality: Option<Chirality>,
}

impl NodeBuilder {
    pub(crate) fn new(
        element: AtomSymbol,
        charge: i8,
        isotope: Option<u16>,
        aromatic: Option<bool>,
        hydrogens: Option<u8>,
        class: Option<u16>,
        chirality: Option<Chirality>,
    ) -> Result<NodeBuilder, NodeError> {
        let atom = Atom::new(element, charge, isotope)?;

        Ok(NodeBuilder {
            atom,
            aromatic,
            hydrogens,
            class,
            chirality,
        })
    }

    pub(crate) fn aromatic(&self) -> Option<bool> {
        self.aromatic
    }

    pub(crate) fn set_hydrogens(&mut self, h: u8) -> &mut Self {
        self.hydrogens = Some(h);
        self
    }

    pub(crate) fn build(mut self, bond_order_sum: Option<u8>) -> Result<Node, NodeError> {
        if self.hydrogens.is_none() {
            let aromatic = self.aromatic.unwrap_or(false);
            self.set_hydrogens(self.atom.implicit_hydrogens(bond_order_sum, aromatic)?);
        }

        Node::new(
            self.atom,
            self.aromatic.ok_or(NodeError::UndefinedAromatic)?,
            self.hydrogens.ok_or(NodeError::UndefinedHydrogen)?,
            self.class,
            self.chirality,
        )
    }
}
