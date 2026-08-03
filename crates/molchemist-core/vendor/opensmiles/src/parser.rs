use std::collections::HashMap;
use std::iter::Peekable;
use std::str::{Chars, FromStr};

use crate::error::ParserError;
use crate::{
    AtomSymbol, BondType, Chirality, Molecule, MoleculeBuilder, NodeError, NodeIndex, OrganicAtom,
};

struct Parser<'a> {
    chars: Peekable<Chars<'a>>,
    position: usize,
    builder: MoleculeBuilder,
    next_bond_type: Option<BondType>,
    next_bond_source: Option<NodeIndex>,
    branch_bond_type: Option<BondType>,
    cycles_target: HashMap<u8, (NodeIndex, Option<BondType>)>, // (node_index, bond_type_at_open)
    node_offset: NodeIndex, // Offset for global node indexing (used in branches)
    deferred_ring_bonds: Vec<(NodeIndex, NodeIndex, BondType)>, // (main_target, local_source, bond_type)
    pending_dot: bool,
    allow_terminator: bool,
}

impl<'a> Parser<'a> {
    fn new(input: &'a str) -> Self {
        Parser {
            chars: input.chars().peekable(),
            position: 0,
            builder: MoleculeBuilder::new(),
            next_bond_type: None,
            next_bond_source: None,
            branch_bond_type: None,
            cycles_target: HashMap::new(),
            node_offset: 0,
            deferred_ring_bonds: Vec::new(),
            pending_dot: false,
            allow_terminator: true,
        }
    }

    fn new_with_offset(
        input: &'a str,
        position_offset: usize,
        node_offset: NodeIndex,
        cycles_target: HashMap<u8, (NodeIndex, Option<BondType>)>,
    ) -> Self {
        Parser {
            chars: input.chars().peekable(),
            position: position_offset,
            builder: MoleculeBuilder::new(),
            next_bond_type: None,
            next_bond_source: None,
            branch_bond_type: None,
            cycles_target,
            node_offset,
            deferred_ring_bonds: Vec::new(),
            pending_dot: false,
            allow_terminator: false,
        }
    }

    fn next(&mut self) -> Option<char> {
        self.position += 1;
        self.chars.next()
    }

    fn peek(&mut self) -> Option<&char> {
        self.chars.peek()
    }

    #[allow(clippy::type_complexity)]
    fn parse(
        mut self,
    ) -> Result<
        (
            MoleculeBuilder,
            Option<BondType>,                           // branch_bond_type
            Option<BondType>, // next_bond_type (for dangling bond detection)
            HashMap<u8, (NodeIndex, Option<BondType>)>, // cycles_target
            Vec<(NodeIndex, NodeIndex, BondType)>, // deferred_ring_bonds
        ),
        ParserError,
    > {
        while let Some(c) = self.next() {
            // Atom
            if c.is_ascii_alphabetic() || c == '*' {
                let elem = self.parse_element_symbol(c, false)?;
                // Aromaticity is indicated by lowercase letters (c, n, o, etc.)
                // Wildcard '*' outside brackets is non-aromatic by default
                let aromatic = Some(c.is_ascii_lowercase());
                self.builder
                    .add_atom(elem, 0, None, aromatic, None, None, None)?;
                self.connect_current_atom()?;
                self.pending_dot = false;
            // Brackets Atom
            } else if c == '[' {
                let (elem, charge, isotope, aromatic, hydrogen, class, chirality) =
                    self.parse_bracket_atom()?;
                self.builder
                    .add_atom(elem, charge, isotope, aromatic, hydrogen, class, chirality)?;
                self.connect_current_atom()?;
                self.pending_dot = false;

            // Dot separator (disconnected fragments) — resets the chain
            } else if c == '.' {
                if self.pending_dot || self.next_bond_type.is_some() {
                    return Err(ParserError::MisplacedDot);
                }
                if self.next_bond_source.is_none() {
                    if !self.builder.nodes().is_empty() || self.branch_bond_type.is_some() {
                        return Err(ParserError::MisplacedDot);
                    }
                    self.branch_bond_type = Some(BondType::Disconnected);
                }
                self.next_bond_source = None;
                self.pending_dot = true;

            // Explicit bond
            } else if c == '-'
                || c == '='
                || c == '#'
                || c == '$'
                || c == ':'
                || c == '/'
                || c == '\\'
            {
                if self.pending_dot {
                    return Err(ParserError::MisplacedDot);
                }
                if self.next_bond_type.is_some() {
                    return Err(ParserError::ConsecutiveBondSymbols);
                }
                self.next_bond_type = Some(BondType::try_from(&c)?);
                if self.builder.nodes().is_empty() {
                    self.branch_bond_type = self.next_bond_type;
                    self.next_bond_type = None;
                }

            // Branches
            } else if c == '(' {
                if self.next_bond_source.is_none() {
                    return Err(ParserError::BranchWithoutPrecedingAtom);
                }
                self.parse_branch()?;
            // cycles
            } else if c == '%' || c.is_ascii_digit() {
                let cycle_number: u8 = if c == '%' {
                    let first = self.next().ok_or(ParserError::UnexpectedEndOfInput(
                        "cycle number".to_string(),
                    ))?;
                    let second = self.next().ok_or(ParserError::UnexpectedEndOfInput(
                        "cycle number".to_string(),
                    ))?;
                    let first_u8: u8 = first
                        .to_digit(10)
                        .ok_or(ParserError::UnexpectedCharacter(first, self.position))?
                        as u8;
                    let second_u8: u8 = second
                        .to_digit(10)
                        .ok_or(ParserError::UnexpectedCharacter(second, self.position))?
                        as u8;
                    first_u8 * 10 + second_u8
                } else {
                    c.to_digit(10).expect("Unreachable error") as u8
                };

                // If the key already exists, close the ring
                if let Some((target, bond_type_at_open)) =
                    self.cycles_target.get(&cycle_number).copied()
                {
                    let local_index = self
                        .next_bond_source
                        .ok_or(ParserError::UnexpectedCharacter(c, self.position))?;
                    let global_source = self.node_offset + local_index;
                    let bond_type_at_close = self.next_bond_type.take();

                    // An atom cannot be bonded to itself (e.g., C11)
                    if global_source == target {
                        return Err(ParserError::SelfBond(cycle_number));
                    }

                    // Two atoms cannot be joined by more than one bond (e.g., C12CCCCC12)
                    if self.has_bond_between(target, global_source) {
                        return Err(ParserError::DuplicateBond(target, global_source));
                    }

                    // Validate ring bond types match if both are explicitly specified
                    if let (Some(open), Some(close)) = (bond_type_at_open, bond_type_at_close) {
                        if open != close {
                            return Err(ParserError::MismatchedRingBond(cycle_number));
                        }
                    }

                    // Check if the target is from a parent parser (before our node_offset)
                    if target < self.node_offset {
                        // Defer this bond - it connects to a node in the parent
                        // Use explicit bond type if specified, otherwise Simple
                        // (aromaticity will be determined by parent when creating the bond)
                        let ring_bond_type = bond_type_at_open
                            .or(bond_type_at_close)
                            .unwrap_or(BondType::Simple);
                        self.deferred_ring_bonds
                            .push((target, local_index, ring_bond_type));
                    } else {
                        // Target is within this parser's nodes
                        // Determine bond type from explicit specification or aromaticity
                        let ring_bond_type =
                            bond_type_at_open.or(bond_type_at_close).unwrap_or_else(|| {
                                // Adjust target index for local builder space
                                let local_target = target - self.node_offset;
                                let source_aromatic =
                                    self.builder.nodes()[local_index as usize].aromatic();
                                let target_aromatic =
                                    self.builder.nodes()[local_target as usize].aromatic();
                                if source_aromatic == Some(true) && target_aromatic == Some(true) {
                                    BondType::Aromatic
                                } else {
                                    BondType::Simple
                                }
                            });
                        self.connect_ring_closure(target, ring_bond_type)?;
                    }
                    // Remove the key so it can be reused
                    self.cycles_target.remove(&cycle_number);
                // Otherwise, this is the start of a new ring
                } else {
                    let local_index = self
                        .next_bond_source
                        .ok_or(ParserError::UnexpectedCharacter(c, self.position))?;
                    let global_index = self.node_offset + local_index;
                    let bond_type_at_open = self.next_bond_type.take();
                    self.cycles_target
                        .insert(cycle_number, (global_index, bond_type_at_open));
                }
            // Whitespace terminates the SMILES string (OpenSMILES spec)
            } else if c == ' ' || c == '\t' || c == '\n' || c == '\r' {
                if self.allow_terminator {
                    break;
                }
                return Err(ParserError::UnexpectedCharacter(c, self.position));
            } else {
                return Err(ParserError::UnexpectedCharacter(c, self.position));
            }
        }

        if self.pending_dot {
            return Err(ParserError::MisplacedDot);
        }

        Ok((
            self.builder,
            self.branch_bond_type,
            self.next_bond_type,
            self.cycles_target,
            self.deferred_ring_bonds,
        ))
    }

    fn parse_branch(&mut self) -> Result<(), ParserError> {
        let mut s = String::new();
        let mut parenthesis_count: i32 = 1;
        let position = self.position;
        while let Some(c) = self.next() {
            if c == '(' {
                parenthesis_count += 1;
            }

            if c == ')' {
                parenthesis_count -= 1;

                if parenthesis_count == 0 {
                    break;
                }
            }
            s.push(c);
        }
        if parenthesis_count > 0 {
            return Err(ParserError::UnclosedParenthesis);
        }
        if parenthesis_count < 0 {
            return Err(ParserError::UnopenedParenthesis);
        }
        if s.is_empty() {
            return Err(ParserError::EmptyBranch);
        }

        // Calculate the global node offset for the branch
        let branch_node_offset = self
            .node_offset
            .checked_add(self.builder.nodes().len() as NodeIndex)
            .ok_or(ParserError::TooManyNodes)?;

        // Pass cycles_target to branch so rings can span branch boundaries
        let branch_parser = Parser::new_with_offset(
            &s,
            position,
            branch_node_offset,
            std::mem::take(&mut self.cycles_target),
        );
        let (
            branch_builder,
            branch_bond_type,
            branch_next_bond_type,
            updated_cycles,
            deferred_bonds,
        ) = branch_parser.parse()?;

        if branch_builder.nodes().is_empty() || branch_next_bond_type.is_some() {
            return Err(ParserError::BondWithoutFollowingAtom);
        }

        // Restore the updated cycles_target
        self.cycles_target = updated_cycles;

        // Add the branch to the main builder
        // A dot at the start of a branch means disconnected — no bond to parent
        let bond_type = branch_bond_type.unwrap_or(BondType::Simple);
        let connect_source = if bond_type == BondType::Disconnected {
            None
        } else {
            self.next_bond_source
        };
        self.builder
            .add_branch(branch_builder, bond_type, connect_source);

        // Create deferred ring bonds (rings opened in parent, closed in branch)
        for (main_target, branch_local_source, mut ring_bond_type) in deferred_bonds {
            // branch_local_source needs to be adjusted to main molecule space
            let main_source = branch_node_offset + branch_local_source;
            let local_source = main_source - self.node_offset;

            // If the deferred bond still points outside this parser's builder,
            // bubble it upward again so an ancestor can connect it.
            if main_target < self.node_offset {
                self.deferred_ring_bonds
                    .push((main_target, local_source, ring_bond_type));
                continue;
            }
            let local_target = main_target - self.node_offset;

            // Fix: if no explicit bond type was specified and both atoms are aromatic,
            // the implicit bond should be Aromatic, not Simple
            if ring_bond_type == BondType::Simple {
                let target_aromatic = self
                    .builder
                    .nodes()
                    .get(local_target as usize)
                    .and_then(|node| node.aromatic());
                let source_aromatic = self
                    .builder
                    .nodes()
                    .get(local_source as usize)
                    .and_then(|node| node.aromatic());
                if target_aromatic == Some(true) && source_aromatic == Some(true) {
                    ring_bond_type = BondType::Aromatic;
                }
            }

            self.builder
                .add_bond(local_target, local_source, ring_bond_type);
        }

        if self.next_bond_source.is_none() {
            self.next_bond_source = Some(0);
        }
        Ok(())
    }

    /// Parse element symbol. `in_bracket` controls whether all two-letter elements
    /// are allowed (true) or only organic subset Cl/Br (false).
    ///
    /// Returns the parsed `AtomSymbol` directly, avoiding intermediate String allocations.
    fn parse_element_symbol(
        &mut self,
        c: char,
        in_bracket: bool,
    ) -> Result<AtomSymbol, ParserError> {
        if c == '*' {
            return Ok(AtomSymbol::Wildcard);
        }

        if c.is_ascii_uppercase() {
            if let Some(&next_c) = self.peek() {
                if next_c.is_ascii_lowercase() {
                    if in_bracket {
                        // In brackets, all valid two-letter elements are allowed
                        let buf = [c as u8, next_c as u8];
                        let two_letter = std::str::from_utf8(&buf).unwrap();
                        if let Ok(sym) = AtomSymbol::from_str(two_letter) {
                            self.next();
                            return Ok(sym);
                        }
                    } else {
                        // Outside brackets, only Cl and Br are valid two-letter elements
                        match (c, next_c) {
                            ('C', 'l') => {
                                self.next();
                                return Ok(AtomSymbol::Organic(OrganicAtom::Cl));
                            }
                            ('B', 'r') => {
                                self.next();
                                return Ok(AtomSymbol::Organic(OrganicAtom::Br));
                            }
                            _ => {}
                        }
                    }
                }
            }
        } else if in_bracket && c.is_ascii_lowercase() {
            // Aromatic two-letter symbols: se, as, te (OpenSMILES spec)
            if let Some(&next_c) = self.peek() {
                if next_c.is_ascii_lowercase() {
                    let buf = [c as u8, next_c as u8];
                    let two_letter = std::str::from_utf8(&buf).unwrap();
                    if let Ok(sym) = AtomSymbol::from_str(two_letter) {
                        self.next();
                        return Ok(sym);
                    }
                }
            }
        }

        // Single character element
        let buf = [c.to_ascii_uppercase() as u8];
        let s = std::str::from_utf8(&buf).unwrap();
        AtomSymbol::from_str(s).map_err(|e| ParserError::NodeError(NodeError::AtomError(e)))
    }

    #[allow(clippy::type_complexity)]
    fn parse_bracket_atom(
        &mut self,
    ) -> Result<
        (
            AtomSymbol,
            i8,
            Option<u16>,
            Option<bool>,
            Option<u8>,
            Option<u16>,
            Option<Chirality>,
        ),
        ParserError,
    > {
        let isotope = self.parse_isotope()?;

        let first_char = self.next().ok_or(ParserError::UnexpectedEndOfInput(
            "Element identifier".to_string(),
        ))?;
        if !first_char.is_ascii_alphabetic() && first_char != '*' {
            return Err(ParserError::MissingElementInBracketAtom);
        }
        let elem = self.parse_element_symbol(first_char, true)?;

        let mut chirality = None;
        let mut hydrogen: Option<u8> = None;
        let mut charge: i8 = 0;
        let mut class = None;
        let mut property_rank = 0u8;

        loop {
            match self.peek() {
                Some(&']') | None => break,
                Some(&'@') => {
                    if property_rank >= 1 {
                        return Err(ParserError::InvalidBracketPropertyOrder(self.position + 1));
                    }
                    property_rank = 1;
                    chirality = self.parse_chirality()?;
                }
                Some(&'H') => {
                    if property_rank >= 2 {
                        return Err(ParserError::InvalidBracketPropertyOrder(self.position + 1));
                    }
                    property_rank = 2;
                    hydrogen = self.parse_hydrogen()?;
                }
                Some(&'+') | Some(&'-') => {
                    if property_rank >= 3 {
                        return Err(ParserError::InvalidBracketPropertyOrder(self.position + 1));
                    }
                    property_rank = 3;
                    charge = self.parse_charge()?;
                }
                Some(&':') => {
                    if property_rank >= 4 {
                        return Err(ParserError::InvalidBracketPropertyOrder(self.position + 1));
                    }
                    class = Some(self.parse_class()?);
                    break;
                }
                Some(&c) => {
                    self.next();
                    return Err(ParserError::UnexpectedCharacter(c, self.position));
                }
            }
        }

        let hydrogen = hydrogen.or(Some(0));

        match self.next() {
            Some(']') => (),
            None => Err(ParserError::UnexpectedEndOfInput("]".to_string()))?,
            Some(c) => Err(ParserError::UnexpectedCharacter(c, self.position))?,
        }

        // A hydrogen atom cannot have a hydrogen count (e.g., [HH1] is illegal)
        if elem == AtomSymbol::H {
            if let Some(h) = hydrogen {
                if h > 0 {
                    return Err(ParserError::HydrogenWithHydrogenCount);
                }
            }
        }

        // Aromaticity is determined by whether the element was written in lowercase
        let aromatic = Some(first_char.is_ascii_lowercase());

        Ok((elem, charge, isotope, aromatic, hydrogen, class, chirality))
    }

    fn parse_isotope(&mut self) -> Result<Option<u16>, ParserError> {
        let mut builder = String::new();
        while self.peek().is_some_and(|c| c.is_ascii_digit()) {
            builder.push(self.next().unwrap());
        }
        if builder.is_empty() {
            Ok(None)
        } else {
            builder
                .parse::<u16>()
                .map(Some)
                .map_err(|_| ParserError::IsotopeOutOfRange(builder))
        }
    }

    fn parse_class(&mut self) -> Result<u16, ParserError> {
        self.next(); // consume ':'
        let mut builder = String::new();
        while self.peek().is_some_and(|c| c.is_ascii_digit()) {
            builder.push(self.next().unwrap());
        }
        if builder.is_empty() {
            return Err(ParserError::MissingAtomClass);
        }
        let value = builder
            .parse::<u16>()
            .map_err(|_| ParserError::AtomClassOutOfRange(builder.clone()))?;
        if value > 9999 {
            return Err(ParserError::AtomClassOutOfRange(builder));
        }
        Ok(value)
    }

    fn parse_chirality(&mut self) -> Result<Option<Chirality>, ParserError> {
        if self.peek() != Some(&'@') {
            return Ok(None);
        }
        self.next(); // consume first '@'

        match self.peek() {
            Some(&'@') => {
                self.next();
                Ok(Some(Chirality::TH2))
            }
            Some(&'T') => {
                self.next();
                match self.next() {
                    Some('H') => self.parse_chirality_index(1, 2, |n| match n {
                        1 => Some(Chirality::TH1),
                        2 => Some(Chirality::TH2),
                        _ => None,
                    }),
                    Some('B') => self.parse_chirality_index(1, 20, |n| Chirality::tb(n as u8)),
                    Some(c) => Err(ParserError::InvalidChiralitySpec(
                        format!("@T{}", c),
                        self.position,
                    )),
                    None => Err(ParserError::UnexpectedEndOfInput(
                        "chirality class".to_string(),
                    )),
                }
            }
            Some(&'A') => {
                self.next();
                match self.next() {
                    Some('L') => self.parse_chirality_index(1, 2, |n| match n {
                        1 => Some(Chirality::AL1),
                        2 => Some(Chirality::AL2),
                        _ => None,
                    }),
                    Some(c) => Err(ParserError::InvalidChiralitySpec(
                        format!("@A{}", c),
                        self.position,
                    )),
                    None => Err(ParserError::UnexpectedEndOfInput(
                        "chirality class".to_string(),
                    )),
                }
            }
            Some(&'S') => {
                self.next();
                match self.next() {
                    Some('P') => self.parse_chirality_index(1, 3, |n| match n {
                        1 => Some(Chirality::SP1),
                        2 => Some(Chirality::SP2),
                        3 => Some(Chirality::SP3),
                        _ => None,
                    }),
                    Some(c) => Err(ParserError::InvalidChiralitySpec(
                        format!("@S{}", c),
                        self.position,
                    )),
                    None => Err(ParserError::UnexpectedEndOfInput(
                        "chirality class".to_string(),
                    )),
                }
            }
            Some(&'O') => {
                self.next();
                match self.next() {
                    Some('H') => self.parse_chirality_index(1, 30, |n| Chirality::oh(n as u8)),
                    Some(c) => Err(ParserError::InvalidChiralitySpec(
                        format!("@O{}", c),
                        self.position,
                    )),
                    None => Err(ParserError::UnexpectedEndOfInput(
                        "chirality class".to_string(),
                    )),
                }
            }
            _ => Ok(Some(Chirality::TH1)),
        }
    }

    /// Parse a chirality index (1 or 2 digit number) and map it via `f`.
    /// Returns an error if the number is outside `[min, max]` or if `f` returns None.
    fn parse_chirality_index(
        &mut self,
        min: u32,
        max: u32,
        f: impl FnOnce(u32) -> Option<Chirality>,
    ) -> Result<Option<Chirality>, ParserError> {
        let first = self.next().ok_or(ParserError::UnexpectedEndOfInput(
            "chirality index".to_string(),
        ))?;
        let pos = self.position;
        let first_digit = first
            .to_digit(10)
            .ok_or(ParserError::InvalidChiralityClass(first.to_string(), pos))?;

        let n = if let Some(&next_c) = self.peek() {
            if let Some(second_digit) = next_c.to_digit(10) {
                self.next();
                first_digit * 10 + second_digit
            } else {
                first_digit
            }
        } else {
            first_digit
        };

        if n < min || n > max {
            return Err(ParserError::InvalidChiralityClass(
                n.to_string(),
                self.position,
            ));
        }

        f(n).map(Some)
            .ok_or_else(|| ParserError::InvalidChiralityClass(n.to_string(), self.position))
    }

    fn parse_hydrogen(&mut self) -> Result<Option<u8>, ParserError> {
        match self.peek() {
            None => Err(ParserError::UnexpectedEndOfInput("]".to_string())),
            Some(&'H') => {
                self.next();
                let Some(&digit) = self.peek() else {
                    return Ok(Some(1));
                };
                if !digit.is_ascii_digit() {
                    return Ok(Some(1));
                }
                self.next();
                if self.peek().is_some_and(|next| next.is_ascii_digit()) {
                    let mut invalid = digit.to_string();
                    while self.peek().is_some_and(|next| next.is_ascii_digit()) {
                        invalid.push(self.next().unwrap());
                    }
                    return Err(ParserError::HydrogenOutOfRange(invalid));
                }
                Ok(Some(digit.to_digit(10).unwrap() as u8))
            }
            _ => Ok(Some(0)),
        }
    }

    fn parse_charge(&mut self) -> Result<i8, ParserError> {
        let sign = self
            .next()
            .ok_or_else(|| ParserError::InvalidChargeSyntax(String::new()))?;
        let multiplier = match sign {
            '+' => 1i8,
            '-' => -1i8,
            _ => return Err(ParserError::ChargeWithoutSign),
        };

        if self.peek().is_some_and(|next| *next == sign) {
            self.next();
            return Ok(multiplier * 2);
        }

        let mut digits = String::new();
        while self.peek().is_some_and(|next| next.is_ascii_digit()) {
            digits.push(self.next().unwrap());
        }
        if digits.is_empty() {
            return Ok(multiplier);
        }
        if digits.len() > 2 {
            return Err(ParserError::InvalidChargeSyntax(format!("{sign}{digits}")));
        }
        let magnitude = digits
            .parse::<i8>()
            .map_err(|_| ParserError::ChargeOutOfRange(digits.clone()))?;
        Ok(multiplier * magnitude)
    }

    fn connect_current_atom(&mut self) -> Result<(), ParserError> {
        if self.builder.nodes().is_empty() {
            return Err(ParserError::NoAtomToBond);
        }
        let current_atom = self.get_current_atom_index()?;

        if let Some(src) = self.next_bond_source {
            self.add_bond_between(src, current_atom);
        }
        self.next_bond_source = Some(current_atom);
        Ok(())
    }

    fn connect_ring_closure(
        &mut self,
        target: NodeIndex,
        bond_type: BondType,
    ) -> Result<(), ParserError> {
        if self.builder.nodes().is_empty() {
            return Err(ParserError::NoAtomToBond);
        }
        let current_atom = self.next_bond_source.ok_or(ParserError::NoAtomToBond)?;
        // `target` is a global atom index; convert to local branch-builder index.
        // This is always valid because `connect_ring_closure` is only called when
        // `target >= self.node_offset` (see ring-closure dispatch above).
        let local_target = target - self.node_offset;
        self.builder.add_bond(current_atom, local_target, bond_type);
        Ok(())
    }

    fn get_current_atom_index(&self) -> Result<NodeIndex, ParserError> {
        let current_atom: NodeIndex = (self.builder.nodes().len() - 1)
            .try_into()
            .map_err(|_| ParserError::TooManyNodes)?;
        Ok(current_atom)
    }

    fn add_bond_between(&mut self, source: NodeIndex, target: NodeIndex) {
        // Explicit bonds take priority, otherwise determine implicit bond type
        let bond_type = self.next_bond_type.take().unwrap_or(
            if self.builder.nodes()[source as usize].aromatic() == Some(true)
                && self.builder.nodes()[target as usize].aromatic() == Some(true)
            {
                BondType::Aromatic
            } else {
                BondType::Simple
            },
        );

        self.builder.add_bond(source, target, bond_type);
    }

    /// Check if a bond already exists between two atoms (in either direction).
    fn has_bond_between(&self, a: NodeIndex, b: NodeIndex) -> bool {
        // Check local bonds
        for bond in self.builder.bonds() {
            let (s, t) = (
                self.node_offset + bond.source(),
                self.node_offset + bond.target(),
            );
            if (s == a && t == b) || (s == b && t == a) {
                return true;
            }
        }
        // Check deferred ring bonds
        for &(target, source, _) in &self.deferred_ring_bonds {
            let s = self.node_offset + source;
            if (s == a && target == b) || (s == b && target == a) {
                return true;
            }
        }
        false
    }
}

/// Parses a SMILES string into a [`Molecule`].
///
/// Follows the [OpenSMILES specification](http://opensmiles.org/opensmiles.html).
/// Whitespace (space, tab, newline) terminates the SMILES string as per the spec,
/// so `"CCO ethanol"` parses successfully as ethanol.
///
/// # Errors
///
/// Returns a [`ParserError`] if the input is not valid SMILES.
///
/// # Examples
///
/// ```
/// use molchemist_core::parse;
///
/// let ethanol = parse("CCO").unwrap();
/// let benzene = parse("c1ccccc1").unwrap();
/// let chiral  = parse("[C@H](F)(Cl)Br").unwrap();
///
/// assert!(parse("C(C").is_err()); // unclosed parenthesis
/// ```
pub fn parse(input: &str) -> Result<Molecule, ParserError> {
    let parser = Parser::new(input);
    let (builder, branch_bond_type, next_bond_type, cycles_target, _) = parser.parse()?;

    if builder.nodes().is_empty() {
        return Err(ParserError::EmptyInput);
    }

    // Check for unclosed rings at the top level
    if !cycles_target.is_empty() {
        return Err(ParserError::UnclosedRing(
            cycles_target.into_keys().collect(),
        ));
    }

    // Check for bond at start (e.g., "-C") - only invalid at top level
    if branch_bond_type.is_some() {
        return Err(ParserError::BondWithoutPrecedingAtom);
    }

    // Check for dangling bond at end (e.g., "C=")
    // Note: if branch_bond_type was Some, we already returned above
    if next_bond_type.is_some() {
        return Err(ParserError::BondWithoutFollowingAtom);
    }

    let molecule = builder.build()?;
    require_valid_double_bond_stereo(&molecule)?;

    crate::ast::aromaticity::require_valid_aromaticity(&molecule)?;

    Ok(molecule)
}

fn require_valid_double_bond_stereo(molecule: &Molecule) -> Result<(), ParserError> {
    let mut used_directional_bonds = std::collections::HashSet::new();

    for double_bond in molecule.bonds() {
        if double_bond.kind() != BondType::Double {
            continue;
        }
        let left = double_bond.source();
        let right = double_bond.target();
        let left_markers = directional_markers_at(molecule, left, right);
        let right_markers = directional_markers_at(molecule, right, left);

        require_consistent_directional_markers(left, &left_markers)?;
        require_consistent_directional_markers(right, &right_markers)?;

        if left_markers.is_empty() && right_markers.is_empty() {
            continue;
        }

        if left_markers.is_empty() || right_markers.is_empty() {
            continue;
        }

        used_directional_bonds.extend(left_markers.iter().map(|(bond_idx, _)| *bond_idx));
        used_directional_bonds.extend(right_markers.iter().map(|(bond_idx, _)| *bond_idx));
    }

    for (bond_idx, bond) in molecule.bonds().iter().enumerate() {
        if matches!(bond.kind(), BondType::Up | BondType::Down)
            && !used_directional_bonds.contains(&bond_idx)
        {
            return Err(ParserError::OrphanDirectionalBond(bond_idx));
        }
    }
    Ok(())
}

fn directional_markers_at(
    molecule: &Molecule,
    center: NodeIndex,
    other_center: NodeIndex,
) -> Vec<(usize, i8)> {
    molecule
        .bonds()
        .iter()
        .enumerate()
        .filter_map(|(bond_idx, bond)| {
            let base_sign = match bond.kind() {
                BondType::Up => 1,
                BondType::Down => -1,
                _ => return None,
            };
            if bond.source() == center && bond.target() != other_center {
                Some((bond_idx, base_sign))
            } else if bond.target() == center && bond.source() != other_center {
                Some((bond_idx, -base_sign))
            } else {
                None
            }
        })
        .collect()
}

fn require_consistent_directional_markers(
    center: NodeIndex,
    markers: &[(usize, i8)],
) -> Result<(), ParserError> {
    if markers.len() > 2 || (markers.len() == 2 && markers[0].1 == markers[1].1) {
        return Err(ParserError::ConflictingDoubleBondStereo(center));
    }
    Ok(())
}

#[cfg(test)]
mod strict_tests {
    use super::*;

    #[test]
    fn accepts_the_required_atom_class_range() {
        assert_eq!(parse("[C:1000]").unwrap().nodes()[0].class(), Some(1000));
        assert_eq!(parse("[C:9999]").unwrap().nodes()[0].class(), Some(9999));
    }

    #[test]
    fn rejects_missing_or_overflowing_bracket_numbers() {
        for smiles in ["[C:]", "[C:65536]", "[65536C]"] {
            assert!(parse(smiles).is_err(), "{smiles}");
        }
    }

    #[test]
    fn rejects_empty_input() {
        for smiles in ["", " ", "\tcomment"] {
            assert!(
                matches!(parse(smiles), Err(ParserError::EmptyInput)),
                "{smiles:?}"
            );
        }
    }

    #[test]
    fn rejects_repeated_or_out_of_order_bracket_properties() {
        for smiles in [
            "[C+H]", "[CHH]", "[C@@@]", "[C+-]", "[C++2]", "[C+2+]", "[CH00]",
        ] {
            assert!(parse(smiles).is_err(), "{smiles}");
        }
    }

    #[test]
    fn rejects_relaxed_branch_dot_and_bond_forms_by_default() {
        for smiles in [
            "C((C))O", "(CO)N", "C.", ".C", "C..C", "C(=)O", "C(C=)O", "C-=C", "CC/=C/C",
        ] {
            assert!(parse(smiles).is_err(), "{smiles}");
        }
    }

    #[test]
    fn keeps_strict_disconnected_and_ring_number_forms() {
        assert!(parse("C.C").is_ok());
        assert!(parse("C(.O)N").is_ok());
        assert!(parse("C1.C1").is_ok());
    }

    #[test]
    fn rejects_incomplete_or_conflicting_double_bond_stereo() {
        for smiles in ["C/C=C", "C/C=CC", r"C/C(\F)=C/C", "F/C"] {
            assert!(parse(smiles).is_err(), "{smiles}");
        }
        assert!(parse(r"F/C=C\F").is_ok());
        assert!(parse("F/C=C/CC=CC").is_ok());
    }

    #[test]
    fn rejects_semantically_invalid_aromatic_notation() {
        for smiles in ["c", "CccccC", "c1cccc1", "C:C"] {
            assert!(parse(smiles).is_err(), "{smiles}");
        }
        for smiles in ["c1ccccc1", "o1cccc1", "[nH]1cccc1", "c1ccc2ccccc2c1"] {
            assert!(parse(smiles).is_ok(), "{smiles}");
        }
    }
}
