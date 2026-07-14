// SPDX-License-Identifier: Apache-2.0
// Adapted from Typst 0.15.0's WebAssembly plugin host and modified for molchemist-cli.

use wasmi::{Caller, Config, Engine, ExternType, Instance, Linker, Module, Store, Val, ValType};

use crate::RenderMode;

const CORE_WASM: &[u8] = include_bytes!("../wasm/molchemist_plugin.wasm");
const LAYOUT_WASM: &[u8] = include_bytes!("../wasm/molchemist_smiles_plugin.wasm");

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

    pub fn sdf_to_code(
        &mut self,
        sdf: &str,
        mode: RenderMode,
        atom_sep: &str,
        indent: usize,
    ) -> Result<String, String> {
        let indent = indent.to_string();
        let output = self.core.call(
            "sdf_to_code",
            &[
                sdf.as_bytes(),
                mode.as_str().as_bytes(),
                atom_sep.as_bytes(),
                indent.as_bytes(),
            ],
        )?;
        String::from_utf8(output).map_err(|_| "core plugin returned non-UTF-8 source".to_string())
    }

    pub fn smiles_to_code(
        &mut self,
        smiles: &str,
        mode: RenderMode,
        atom_sep: &str,
        indent: usize,
    ) -> Result<String, String> {
        let layout_function = match mode {
            RenderMode::Full => "smiles_to_full_layout_input",
            RenderMode::Abbreviate | RenderMode::Skeletal => "smiles_to_layout_input",
        };
        let layout_input = self.core.call(layout_function, &[smiles.as_bytes()])?;
        let layout = match &mut self.layout {
            Some(layout) => layout,
            None => self.layout.insert(Plugin::new(LAYOUT_WASM)?),
        };
        let coordinates = layout.call("layout_coordinates", &[&layout_input])?;
        let indent = indent.to_string();
        let output = self.core.call(
            "smiles_to_code",
            &[
                smiles.as_bytes(),
                &coordinates,
                mode.as_str().as_bytes(),
                atom_sep.as_bytes(),
                indent.as_bytes(),
            ],
        )?;
        String::from_utf8(output).map_err(|_| "core plugin returned non-UTF-8 source".to_string())
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
            .smiles_to_code("c1ccccc1", RenderMode::Skeletal, "3em", 2)
            .unwrap();
        assert!(output.starts_with("#let base-sep = 3em\n#skeletize({\n"));
        assert!(output.contains("double(absolute:"));
    }

    #[test]
    fn plugin_errors_are_returned_without_host_noise() {
        let mut generator = Generator::new().unwrap();
        let error = generator
            .smiles_to_code("C(", RenderMode::Full, "3em", 2)
            .unwrap_err();
        assert!(!error.is_empty());
    }
}
