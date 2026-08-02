use std::ffi::{c_char, c_int};
use std::sync::{Mutex, OnceLock};

use crate::{smiles_layout_input, smiles_to_commands_with_coords, Command, RenderMode};

const LAYOUT_MAGIC: &[u8; 4] = b"MCG2";

unsafe extern "C" {
    fn molchemist_coordgen_layout(
        input: *const u8,
        input_size: usize,
        output: *mut u8,
        output_capacity: usize,
        output_size: *mut usize,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
}

pub fn layout_payload(input: &[u8]) -> Result<Vec<u8>, String> {
    let (atom_count, bond_count) = payload_counts(input)?;
    let output_capacity = 12usize
        .checked_add(
            atom_count
                .checked_mul(8)
                .ok_or_else(|| "Coordinate output size overflowed".to_string())?,
        )
        .and_then(|size| size.checked_add(bond_count))
        .ok_or_else(|| "Coordinate output size overflowed".to_string())?;

    static COORDGEN_LOCK: OnceLock<Mutex<()>> = OnceLock::new();
    let _guard = COORDGEN_LOCK
        .get_or_init(|| Mutex::new(()))
        .lock()
        .map_err(|_| "CoordgenLibs lock was poisoned".to_string())?;

    let mut output = vec![0u8; output_capacity];
    let mut output_size = 0usize;
    let mut error = [0u8; 1024];
    let status = unsafe {
        molchemist_coordgen_layout(
            input.as_ptr(),
            input.len(),
            output.as_mut_ptr(),
            output.len(),
            &mut output_size,
            error.as_mut_ptr().cast(),
            error.len(),
        )
    };

    if status != 0 {
        let length = error
            .iter()
            .position(|byte| *byte == 0)
            .unwrap_or(error.len());
        let message = String::from_utf8_lossy(&error[..length]);
        return Err(if message.is_empty() {
            format!("CoordgenLibs failed with status {status}")
        } else {
            message.into_owned()
        });
    }

    if output_size > output.len() {
        return Err("CoordgenLibs returned an invalid output length".to_string());
    }
    output.truncate(output_size);
    Ok(output)
}

pub fn smiles_to_commands(smiles: &str, mode: RenderMode) -> Result<Vec<Command>, String> {
    let layout_input = smiles_layout_input(smiles, mode)?;
    let coordinates = layout_payload(&layout_input)?;
    smiles_to_commands_with_coords(smiles, &coordinates, mode)
}

fn payload_counts(input: &[u8]) -> Result<(usize, usize), String> {
    if input.len() < 12 || &input[..4] != LAYOUT_MAGIC {
        return Err("Invalid layout payload".to_string());
    }

    let atom_count = u32::from_le_bytes(input[4..8].try_into().unwrap()) as usize;
    let bond_count = u32::from_le_bytes(input[8..12].try_into().unwrap()) as usize;
    Ok((atom_count, bond_count))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        sdf_record_to_commands_with_coords, sdf_record_to_layout_input, smiles_to_layout_input,
    };

    const COLLAPSED_SDF: &str = concat!(
        "collapsed\n",
        "  molchemist\n",
        "\n",
        "  3  2  0  0  0  0  0  0  0  0999 V2000\n",
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
        "    0.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n",
        "    0.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n",
        "  1  2  1  1  0  0  0\n",
        "  1  3  1  0  0  0  0\n",
        "M  END\n",
    );

    fn stereo_bonds(commands: &[Command], output: &mut Vec<String>) {
        for command in commands {
            match command {
                Command::Fragment { links, .. } => output.extend(
                    links
                        .iter()
                        .filter(|link| link.bond_type.starts_with("cram-"))
                        .map(|link| link.bond_type.clone()),
                ),
                Command::Bond { bond_type, .. } if bond_type.starts_with("cram-") => {
                    output.push(bond_type.clone());
                }
                Command::Branch { body } => stereo_bonds(body, output),
                _ => {}
            }
        }
    }

    fn bond_types(commands: &[Command], output: &mut Vec<String>) {
        for command in commands {
            match command {
                Command::Fragment { links, .. } => {
                    output.extend(links.iter().map(|link| link.bond_type.clone()));
                }
                Command::Bond { bond_type, .. } => output.push(bond_type.clone()),
                Command::Branch { body } => bond_types(body, output),
                Command::ComponentBreak => {}
            }
        }
    }

    fn annotations(commands: &[Command], output: &mut Vec<String>) {
        for command in commands {
            match command {
                Command::Fragment {
                    annotation: Some(annotation),
                    ..
                } => output.push(annotation.clone()),
                Command::Branch { body } => annotations(body, output),
                _ => {}
            }
        }
    }

    #[test]
    fn native_coordgen_returns_expected_payload_shape() {
        let input = smiles_to_layout_input(b"N[C@@H](C)C(=O)O").unwrap();
        let output = layout_payload(&input).unwrap();

        assert_eq!(&output[..4], b"MCC2");
        assert_eq!(u32::from_le_bytes(output[4..8].try_into().unwrap()), 6);
        assert_eq!(u32::from_le_bytes(output[8..12].try_into().unwrap()), 5);
        assert_eq!(output.len(), 12 + 6 * 8 + 5);
    }

    #[test]
    fn native_coordgen_accepts_smiles_quadruple_bonds_in_every_mode() {
        for mode in [
            RenderMode::Full,
            RenderMode::Abbreviate,
            RenderMode::Skeletal,
        ] {
            let commands = smiles_to_commands("[Cr]$[Cr]", mode).unwrap();
            let mut bonds = Vec::new();
            bond_types(&commands, &mut bonds);
            assert_eq!(bonds, vec!["quadruple"], "mode {mode:?}");
        }
    }

    #[test]
    fn leading_tetrahedral_enantiomers_produce_opposite_bond_styles() {
        let clockwise = smiles_to_commands("[C@H](F)(Cl)Br", RenderMode::Skeletal).unwrap();
        let counterclockwise = smiles_to_commands("[C@@H](F)(Cl)Br", RenderMode::Skeletal).unwrap();
        let mut clockwise_stereo = Vec::new();
        let mut counterclockwise_stereo = Vec::new();
        stereo_bonds(&clockwise, &mut clockwise_stereo);
        stereo_bonds(&counterclockwise, &mut counterclockwise_stereo);

        assert_eq!(clockwise_stereo, vec!["cram-dashed-left"]);
        assert_eq!(counterclockwise_stereo, vec!["cram-filled-left"]);
    }

    #[test]
    fn alanine_enantiomers_produce_expected_absolute_bond_styles() {
        let l_alanine = smiles_to_commands("N[C@@H](C)C(=O)O", RenderMode::Skeletal).unwrap();
        let d_alanine = smiles_to_commands("N[C@H](C)C(=O)O", RenderMode::Skeletal).unwrap();
        let mut l_stereo = Vec::new();
        let mut d_stereo = Vec::new();
        stereo_bonds(&l_alanine, &mut l_stereo);
        stereo_bonds(&d_alanine, &mut d_stereo);

        assert_eq!(l_stereo, vec!["cram-filled-right"]);
        assert_eq!(d_stereo, vec!["cram-dashed-right"]);
    }

    #[test]
    fn implicit_and_explicit_hydrogen_alanine_keep_the_same_absolute_configuration() {
        let implicit = smiles_to_commands("N[C@@H](C)C(=O)O", RenderMode::Skeletal).unwrap();
        let explicit = smiles_to_commands("N[C@@]([H])(C)C(=O)O", RenderMode::Skeletal).unwrap();
        let mut implicit_stereo = Vec::new();
        let mut explicit_stereo = Vec::new();
        stereo_bonds(&implicit, &mut implicit_stereo);
        stereo_bonds(&explicit, &mut explicit_stereo);

        assert_eq!(implicit_stereo.first(), explicit_stereo.first());
        assert_eq!(
            implicit_stereo.first().map(String::as_str),
            Some("cram-filled-right")
        );
    }

    #[test]
    fn full_mode_keeps_the_expanded_stereochemical_hydrogen_visible() {
        let commands = smiles_to_commands("[C@H](F)(Cl)Br", RenderMode::Full).unwrap();
        let mut stereo = Vec::new();
        stereo_bonds(&commands, &mut stereo);

        assert!(stereo.len() >= 1);
        assert!(commands.iter().any(|command| matches!(
            command,
            Command::Fragment { element, name, .. } if element == "H" && name == "a4"
        )));
    }

    #[test]
    fn extended_chirality_uses_native_geometry_when_topology_is_supported() {
        let cases = [
            ("[Pt@SP2](F)(Cl)(Br)I", 0),
            ("[As@TB5](F)(Cl)(Br)(N)S", 2),
            ("[Co@OH5](F)(Cl)(Br)(I)(N)S", 2),
            ("NC(Br)=[C@AL1]=C(O)C", 2),
        ];

        for (smiles, expected_stereo_bonds) in cases {
            let commands = smiles_to_commands(smiles, RenderMode::Skeletal).unwrap();
            let mut stereo = Vec::new();
            let mut fallback_annotations = Vec::new();
            stereo_bonds(&commands, &mut stereo);
            annotations(&commands, &mut fallback_annotations);

            assert_eq!(stereo.len(), expected_stereo_bonds, "{smiles}");
            assert!(fallback_annotations.is_empty(), "{smiles}");
        }
    }

    #[test]
    fn native_coordgen_recovers_a_collapsed_sdf_layout() {
        let layout_input = sdf_record_to_layout_input(COLLAPSED_SDF.as_bytes(), 1).unwrap();
        let coordinates = layout_payload(&layout_input).unwrap();
        let commands = sdf_record_to_commands_with_coords(
            COLLAPSED_SDF,
            &coordinates,
            RenderMode::Skeletal,
            1,
        )
        .unwrap();
        let mut stereo = Vec::new();
        stereo_bonds(&commands, &mut stereo);

        assert_eq!(stereo.len(), 1);
        assert!(commands.iter().any(|command| matches!(
            command,
            Command::Bond { length_scale, .. } if *length_scale > 0.5
        )));
    }
}
