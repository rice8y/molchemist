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
    use crate::smiles_to_layout_input;

    #[test]
    fn native_coordgen_returns_expected_payload_shape() {
        let input = smiles_to_layout_input(b"N[C@@H](C)C(=O)O").unwrap();
        let output = layout_payload(&input).unwrap();

        assert_eq!(&output[..4], b"MCC2");
        assert_eq!(u32::from_le_bytes(output[4..8].try_into().unwrap()), 6);
        assert_eq!(u32::from_le_bytes(output[8..12].try_into().unwrap()), 5);
        assert_eq!(output.len(), 12 + 6 * 8 + 5);
    }
}
