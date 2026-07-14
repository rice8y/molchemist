use wasm_minimal_protocol::*;

initiate_protocol!();

#[wasm_func]
pub fn sdf_to_ast(sdf_data: &[u8], options: &[u8]) -> Result<Vec<u8>, String> {
    molchemist_core::sdf_to_ast(sdf_data, options)
}

#[wasm_func]
pub fn smiles_to_layout_input(smiles_data: &[u8]) -> Result<Vec<u8>, String> {
    molchemist_core::smiles_to_layout_input(smiles_data)
}

#[wasm_func]
pub fn smiles_to_full_layout_input(smiles_data: &[u8]) -> Result<Vec<u8>, String> {
    molchemist_core::smiles_to_full_layout_input(smiles_data)
}

#[wasm_func]
pub fn smiles_to_ast(
    smiles_data: &[u8],
    coords_data: &[u8],
    options: &[u8],
) -> Result<Vec<u8>, String> {
    molchemist_core::smiles_to_ast(smiles_data, coords_data, options)
}

#[wasm_func]
pub fn sdf_to_code(
    sdf_data: &[u8],
    options: &[u8],
    base_sep: &[u8],
    indent: &[u8],
) -> Result<Vec<u8>, String> {
    let sdf = std::str::from_utf8(sdf_data).map_err(|error| error.to_string())?;
    let base_sep = std::str::from_utf8(base_sep).map_err(|error| error.to_string())?;
    let indent = parse_indent(indent)?;
    let commands = molchemist_core::sdf_to_commands(sdf, render_mode(options))?;
    Ok(molchemist_core::format_alchemist(&commands, base_sep, indent).into_bytes())
}

#[wasm_func]
pub fn smiles_to_code(
    smiles_data: &[u8],
    coords_data: &[u8],
    options: &[u8],
    base_sep: &[u8],
    indent: &[u8],
) -> Result<Vec<u8>, String> {
    let smiles = std::str::from_utf8(smiles_data).map_err(|error| error.to_string())?;
    let base_sep = std::str::from_utf8(base_sep).map_err(|error| error.to_string())?;
    let indent = parse_indent(indent)?;
    let commands =
        molchemist_core::smiles_to_commands_with_coords(smiles, coords_data, render_mode(options))?;
    Ok(molchemist_core::format_alchemist(&commands, base_sep, indent).into_bytes())
}

fn render_mode(options: &[u8]) -> molchemist_core::RenderMode {
    std::str::from_utf8(options)
        .ok()
        .and_then(|mode| molchemist_core::RenderMode::parse(mode).ok())
        .unwrap_or_default()
}

fn parse_indent(indent: &[u8]) -> Result<usize, String> {
    let indent = std::str::from_utf8(indent).map_err(|error| error.to_string())?;
    let indent = indent
        .parse::<usize>()
        .map_err(|_| "Indent width must be an integer".to_string())?;
    if (1..=8).contains(&indent) {
        Ok(indent)
    } else {
        Err("Indent width must be from 1 through 8".to_string())
    }
}
