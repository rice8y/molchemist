#ifndef MOLCHEMIST_SMILES_PROTOCOL_H
#define MOLCHEMIST_SMILES_PROTOCOL_H

#include <stddef.h>
#include <stdint.h>

#ifndef EMSCRIPTEN_KEEPALIVE
#include <emscripten.h>
#endif

#ifndef PROTOCOL_FUNCTION
#define PROTOCOL_FUNCTION __attribute__((import_module("typst_env"))) extern "C"
#endif

PROTOCOL_FUNCTION void wasm_minimal_protocol_send_result_to_host(const uint8_t* ptr, size_t len);
PROTOCOL_FUNCTION void wasm_minimal_protocol_write_args_to_buffer(uint8_t* ptr);

#endif
