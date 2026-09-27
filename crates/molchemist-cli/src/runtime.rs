// SPDX-License-Identifier: Apache-2.0
// Adapted from Typst 0.15.0's WebAssembly plugin host and modified for molchemist-cli.

use wasmi::{Caller, Config, Engine, ExternType, Instance, Linker, Module, Store, Val, ValType};

use crate::RenderMode;

const CORE_WASM: &[u8] = include_bytes!("../wasm/molchemist_plugin.wasm");
const LAYOUT_WASM: &[u8] = include_bytes!("../wasm/molchemist_smiles_plugin.wasm");

#[cfg(test)]
#[path = "wasm_optimization_tests.rs"]
mod optimization_tests;

/// Generates the same Alchemist source as molchemist's Typst plugins.
pub struct Generator {
    core: Plugin,
    layout: Option<Plugin>,
}

impl Generator {
    pub fn new() -> Result<Self, String> {
        Ok(Self {
            core: Plugin::new(CORE_WASM)?,
            layout: None,
        })
    }

    pub fn mol2_to_sdf(&mut self, data: &str, record: usize) -> Result<String, String> {
        let output = self.core.call(
            "mol2_to_sdf",
            &[data.as_bytes(), record.to_string().as_bytes()],
        )?;
        String::from_utf8(output).map_err(|e| e.to_string())
    }

    pub fn inspect_reaction(
        &mut self,
        data: &str,
        format: &str,
        infer: bool,
        limit: usize,
    ) -> Result<serde_json::Value, String> {
        let output = self.core.call(
            "reaction_inspection",
            &[
                data.as_bytes(),
                format.as_bytes(),
                if infer { b"true" } else { b"false" },
                limit.to_string().as_bytes(),
            ],
        )?;
        ciborium::from_reader(output.as_slice()).map_err(|e| e.to_string())
    }

    pub fn coordinate_to_code(
        &mut self,
        data: &str,
        format: &str,
        mode: RenderMode,
        record: usize,
        atom_sep: &str,
        config: &serde_json::Value,
    ) -> Result<String, String> {
        let mut config = config.clone();
        let record = record.to_string();
        let ast = if format == "smiles" {
            let function = if mode == RenderMode::Full {
                "smiles_to_full_layout_input"
            } else {
                "smiles_to_layout_input"
            };
            let input = self.core.call(function, &[data.as_bytes()])?;
            let coords = self.layout_coordinates(&input)?;
            self.core.call(
                "smiles_to_ast",
                &[data.as_bytes(), &coords, mode.as_str().as_bytes()],
            )?
        } else {
            let function = if config["infer-stereo"] == true {
                "sdf_stereo3d_layout_input"
            } else if config["layout"] == "reflow" {
                "sdf_force_layout_input"
            } else {
                "sdf_record_to_layout_input"
            };
            let input = self
                .core
                .call(function, &[data.as_bytes(), record.as_bytes()])?;
            let coords = if input.is_empty() {
                Vec::new()
            } else {
                self.layout_coordinates(&input)?
            };
            let positions = self.core.call(
                "sdf_depiction_positions",
                &[data.as_bytes(), &coords, record.as_bytes()],
            )?;
            config["positions"] =
                ciborium::from_reader::<serde_json::Value, _>(positions.as_slice())
                    .map_err(|e| e.to_string())?;
            if config["infer-stereo"] == true {
                self.core.call(
                    "sdf_stereo3d_ast",
                    &[
                        data.as_bytes(),
                        &coords,
                        mode.as_str().as_bytes(),
                        record.as_bytes(),
                    ],
                )?
            } else if coords.is_empty() {
                self.core.call(
                    "sdf_record_to_ast",
                    &[data.as_bytes(), mode.as_str().as_bytes(), record.as_bytes()],
                )?
            } else {
                self.core.call(
                    if config["layout"] == "reflow" {
                        "sdf_reoriented_ast"
                    } else {
                        "sdf_record_to_ast_with_coords"
                    },
                    &[
                        data.as_bytes(),
                        &coords,
                        mode.as_str().as_bytes(),
                        record.as_bytes(),
                    ],
                )?
            }
        };
        let mut options = Vec::new();
        ciborium::into_writer(&config, &mut options).map_err(|e| e.to_string())?;
        let output = self
            .core
            .call("coordinate_code", &[&ast, &options, atom_sep.as_bytes()])?;
        String::from_utf8(output).map_err(|e| e.to_string())
    }

    pub fn expand_superatoms(&mut self, data: &str, record: usize) -> Result<String, String> {
        let output = self.core.call(
            "expand_superatoms",
            &[data.as_bytes(), record.to_string().as_bytes()],
        )?;
        String::from_utf8(output).map_err(|e| e.to_string())
    }

    pub fn inspect_rgroup(&mut self, data: &str) -> Result<serde_json::Value, String> {
        let output = self.core.call("rgroup_inspection", &[data.as_bytes()])?;
        ciborium::from_reader(output.as_slice()).map_err(|e| e.to_string())
    }

    pub fn rgroup_to_code(
        &mut self,
        data: &str,
        mode: RenderMode,
        atom_sep: &str,
        config: &serde_json::Value,
    ) -> Result<String, String> {
        let doc = self.inspect_rgroup(data)?;
        let root = self.coordinate_to_code(
            doc["root"].as_str().ok_or("Missing R-group root")?,
            "mol",
            mode,
            1,
            atom_sep,
            config,
        )?;
        let mut output = String::from_utf8(self.core.call("composition_code", &[])?)
            .map_err(|e| e.to_string())?;
        output.push_str(&format!("\n#let _rgroup-root = [\n{root}\n]\n"));
        let mut names = Vec::new();
        for (i, member) in doc["members"]
            .as_array()
            .ok_or("Missing R-group members")?
            .iter()
            .enumerate()
        {
            let mut style = config.clone();
            let mut points = serde_json::Map::new();
            for entry in member["attachments"]
                .as_array()
                .ok_or("Missing attachment sites")?
            {
                points.insert(format!("a{}", entry[0]), entry[1].clone());
            }
            style["attachment-points"] = points.into();
            let code = self.coordinate_to_code(
                member["data"].as_str().ok_or("Missing R-group member")?,
                "mol",
                mode,
                1,
                atom_sep,
                &style,
            )?;
            output.push_str(&format!("#let _rgroup-{i} = [\n{code}\n]\n"));
            names.push(format!(
                "(group: {}, drawing: _rgroup-{i})",
                member["group"]
            ));
        }
        let mut conditions = Vec::new();
        for condition in doc["conditions"]
            .as_array()
            .ok_or("Missing R-group conditions")?
        {
            conditions.push(format!(
                "(group: {}, thenGroup: {}, restH: {}, occurrence: \"{}\")",
                condition["group"],
                condition["thenGroup"],
                condition["restH"],
                super::escape_typst_string(condition["occurrence"].as_str().unwrap_or_default())
            ));
        }
        let tuple = |values: Vec<String>| {
            if values.is_empty() {
                "()".into()
            } else {
                format!("({},)", values.join(","))
            }
        };
        output.push_str(&format!(
            "#_rgroup-panel(_rgroup-root, {}, {})",
            tuple(names),
            tuple(conditions)
        ));
        Ok(output)
    }

    fn layout_coordinates(&mut self, input: &[u8]) -> Result<Vec<u8>, String> {
        let plugin = match &mut self.layout {
            Some(plugin) => plugin,
            None => self.layout.insert(Plugin::new(LAYOUT_WASM)?),
        };
        plugin.call("layout_coordinates", &[input])
    }

    pub fn reaction_to_code(
        &mut self,
        analysis: &serde_json::Value,
        mode: RenderMode,
        atom_sep: &str,
        config: &serde_json::Value,
        conditions: &str,
    ) -> Result<String, String> {
        let changed = analysis["changes"]
            .as_array()
            .ok_or("Missing reaction changes")?
            .iter()
            .flat_map(|c| c["maps"].as_array().into_iter().flatten())
            .filter_map(serde_json::Value::as_u64)
            .collect::<std::collections::BTreeSet<_>>();
        let mut output = String::from_utf8(self.core.call("composition_code", &[])?)
            .map_err(|e| e.to_string())?;
        output.push('\n');
        let mut groups = Vec::new();
        for (key, atom_side) in [
            ("reactants", "reactant"),
            ("agents", "agent"),
            ("products", "product"),
        ] {
            let molecules = analysis["reaction"][key]
                .as_array()
                .ok_or("Missing reaction side")?;
            let mut names = Vec::new();
            for (index, mol) in molecules.iter().enumerate() {
                let mut style = config.clone();
                let empty = Vec::new();
                let highlights = analysis[format!("{atom_side}Atoms")]
                    .as_array()
                    .unwrap_or(&empty)
                    .iter()
                    .filter(|m| {
                        m["molecule"].as_u64() == Some(index as u64)
                            && m["map"].as_u64().is_some_and(|n| changed.contains(&n))
                    })
                    .map(|m| serde_json::Value::String(format!("a{}", m["atom"])))
                    .collect::<Vec<_>>();
                style["highlight-atoms"] = highlights.into();
                style["highlight-bonds"] = analysis["changes"]
                    .as_array()
                    .unwrap_or(&empty)
                    .iter()
                    .filter(|c| c[atom_side][0].as_u64() == Some(index as u64))
                    .map(|c| {
                        serde_json::json!([
                            format!("a{}", c[atom_side][1]),
                            format!("a{}", c[atom_side][2])
                        ])
                    })
                    .collect::<Vec<_>>()
                    .into();
                let code = self.coordinate_to_code(
                    mol["data"].as_str().ok_or("Missing reaction molecule")?,
                    mol["format"].as_str().ok_or("Missing molecule format")?,
                    mode,
                    1,
                    atom_sep,
                    &style,
                )?;
                let name = format!("_reaction-{key}-{index}");
                output.push_str(&format!("#let {name} = [\n{code}\n]\n"));
                if index > 0 {
                    names.push("box(baseline: 0.35em)[$+$]".to_string());
                }
                names.push(name);
            }
            groups.push(if names.is_empty() {
                "[]".into()
            } else {
                format!("_molecule-row(({},))", names.join(", "))
            });
        }
        let conditions = super::escape_typst_string(conditions);
        let agents = if groups[1] == "[]" {
            "()".into()
        } else {
            format!("({},)", groups[1])
        };
        output.push_str(&format!(
            "#_reaction-row({}, {agents}, {}, \"{conditions}\")",
            groups[0], groups[2]
        ));
        Ok(output)
    }

    /// Export one SDF record using measured-label coordinate layout.
    pub fn sdf_to_code(
        &mut self,
        sdf: &str,
        mode: RenderMode,
        atom_sep: &str,
    ) -> Result<String, String> {
        self.sdf_record_to_code(sdf, mode, 1, atom_sep)
    }

    pub fn inspect_sdf_record_json(
        &mut self,
        sdf: &str,
        record: usize,
        pretty: bool,
    ) -> Result<String, String> {
        let record = record.to_string();
        let output = self.core.call(
            "sdf_record_to_inspection",
            &[sdf.as_bytes(), record.as_bytes()],
        )?;
        let inspection: serde_json::Value = ciborium::from_reader(output.as_slice())
            .map_err(|error| format!("core plugin returned invalid inspection CBOR: {error}"))?;
        if pretty {
            serde_json::to_string_pretty(&inspection).map_err(|error| error.to_string())
        } else {
            serde_json::to_string(&inspection).map_err(|error| error.to_string())
        }
    }

    pub fn sdf_record_to_code(
        &mut self,
        sdf: &str,
        mode: RenderMode,
        record: usize,
        atom_sep: &str,
    ) -> Result<String, String> {
        self.coordinate_to_code(sdf, "sdf", mode, record, atom_sep, &serde_json::json!({}))
    }

    pub fn smiles_to_code(
        &mut self,
        smiles: &str,
        mode: RenderMode,
        atom_sep: &str,
    ) -> Result<String, String> {
        self.coordinate_to_code(smiles, "smiles", mode, 1, atom_sep, &serde_json::json!({}))
    }
}

struct Plugin {
    instance: Instance,
    store: Store<CallData>,
}

impl Plugin {
    fn new(bytes: &[u8]) -> Result<Self, String> {
        let mut config = Config::default();
        config.wasm_relaxed_simd(false);
        let engine = Engine::new(&config);
        let module = Module::new(&engine, bytes)
            .map_err(|error| format!("failed to load embedded WebAssembly module: {error}"))?;

        if !matches!(module.get_export("memory"), Some(ExternType::Memory(_))) {
            return Err("embedded WebAssembly module does not export memory".to_string());
        }

        let mut linker = Linker::new(&engine);
        linker
            .func_wrap(
                "typst_env",
                "wasm_minimal_protocol_send_result_to_host",
                send_result_to_host,
            )
            .map_err(|error| format!("failed to link plugin result callback: {error}"))?;
        linker
            .func_wrap(
                "typst_env",
                "wasm_minimal_protocol_write_args_to_buffer",
                write_args_to_buffer,
            )
            .map_err(|error| format!("failed to link plugin argument callback: {error}"))?;

        let mut store = Store::new(&engine, CallData::default());
        let instance = linker
            .instantiate_and_start(&mut store, &module)
            .map_err(|error| {
                format!("failed to instantiate embedded WebAssembly module: {error}")
            })?;
        Ok(Self { instance, store })
    }

    fn call(&mut self, function: &str, args: &[&[u8]]) -> Result<Vec<u8>, String> {
        let handle = self
            .instance
            .get_export(&self.store, function)
            .and_then(|export| export.into_func())
            .ok_or_else(|| format!("embedded plugin does not export `{function}`"))?;
        let ty = handle.ty(&self.store);
        if ty.params().iter().any(|value| *value != ValType::I32) || ty.results() != [ValType::I32]
        {
            return Err(format!(
                "embedded plugin function `{function}` has an invalid signature"
            ));
        }
        if ty.params().len() != args.len() {
            return Err(format!(
                "embedded plugin function `{function}` expects {} arguments, but {} were provided",
                ty.params().len(),
                args.len(),
            ));
        }

        let mut lengths = Vec::with_capacity(args.len());
        for arg in args {
            let length = i32::try_from(arg.len())
                .map_err(|_| format!("argument for `{function}` exceeds the WebAssembly limit"))?;
            lengths.push(Val::I32(length));
        }

        let data = self.store.data_mut();
        data.args = args.iter().map(|arg| arg.to_vec()).collect();
        data.output.clear();
        data.memory_error = None;

        let mut status = Val::I32(-1);
        handle
            .call(&mut self.store, &lengths, std::slice::from_mut(&mut status))
            .map_err(|error| format!("embedded plugin `{function}` panicked: {error}"))?;

        if let Some(error) = self.store.data_mut().memory_error.take() {
            return Err(format!(
                "embedded plugin tried to {} outside its memory at {:#x} for {} bytes",
                if error.write { "write" } else { "read" },
                error.offset,
                error.length,
            ));
        }

        let output = std::mem::take(&mut self.store.data_mut().output);
        match status {
            Val::I32(0) => Ok(output),
            Val::I32(1) => match String::from_utf8(output) {
                Ok(message) => Err(message),
                Err(_) => Err(format!(
                    "embedded plugin `{function}` returned a non-UTF-8 error"
                )),
            },
            _ => Err(format!(
                "embedded plugin `{function}` did not respect the Typst plugin protocol"
            )),
        }
    }
}

#[derive(Default)]
struct CallData {
    args: Vec<Vec<u8>>,
    output: Vec<u8>,
    memory_error: Option<MemoryAccessError>,
}

struct MemoryAccessError {
    offset: u32,
    length: u32,
    write: bool,
}

fn write_args_to_buffer(mut caller: Caller<'_, CallData>, pointer: u32) {
    let memory = caller
        .get_export("memory")
        .and_then(|export| export.into_memory())
        .expect("validated plugin memory export");
    let arguments = std::mem::take(&mut caller.data_mut().args);
    let mut offset = pointer as usize;
    for argument in arguments {
        if memory.write(&mut caller, offset, &argument).is_err() {
            caller.data_mut().memory_error = Some(MemoryAccessError {
                offset: offset as u32,
                length: argument.len() as u32,
                write: true,
            });
            return;
        }
        offset += argument.len();
    }
}

fn send_result_to_host(mut caller: Caller<'_, CallData>, pointer: u32, length: u32) {
    let memory = caller
        .get_export("memory")
        .and_then(|export| export.into_memory())
        .expect("validated plugin memory export");
    let mut output = std::mem::take(&mut caller.data_mut().output);
    output.resize(length as usize, 0);
    if memory.read(&caller, pointer as usize, &mut output).is_err() {
        caller.data_mut().memory_error = Some(MemoryAccessError {
            offset: pointer,
            length,
            write: false,
        });
        return;
    }
    caller.data_mut().output = output;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn embedded_plugins_generate_alchemist_source() {
        let mut generator = Generator::new().unwrap();
        let output = generator
            .smiles_to_code("c1ccccc1", RenderMode::Skeletal, "3em")
            .unwrap();
        assert!(output.starts_with("#import \"@preview/alchemist:0.2.0\": *\n"));
        assert!(output.contains("\"bondType\": \"double\""));
    }

    #[test]
    fn plugin_errors_are_returned_without_host_noise() {
        let mut generator = Generator::new().unwrap();
        let error = generator
            .smiles_to_code("C(", RenderMode::Full, "3em")
            .unwrap_err();
        assert!(!error.is_empty());
    }
}
