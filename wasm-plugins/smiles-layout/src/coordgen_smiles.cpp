#include "protocol.h"

#include <cstring>
#include <string>
#include <vector>

#include "coordgen_engine.h"

extern "C" void _initialize();

namespace
{
int send_result(const std::vector<uint8_t>& result)
{
    wasm_minimal_protocol_send_result_to_host(result.data(), result.size());
    return 0;
}
int send_error(const std::string& message)
{
    wasm_minimal_protocol_send_result_to_host(
        reinterpret_cast<const uint8_t*>(message.data()),
        message.size());
    return 1;
}
}

extern "C" EMSCRIPTEN_KEEPALIVE int layout_coordinates(size_t buffer_len)
{
    static bool runtime_initialized = false;
    if (!runtime_initialized)
    {
        _initialize();
        runtime_initialized = true;
    }

    std::vector<uint8_t> input(buffer_len);
    if (buffer_len > 0)
    {
        wasm_minimal_protocol_write_args_to_buffer(input.data());
    }

    std::vector<uint8_t> output;
    std::string error;
    if (!molchemist::layout_coordinates(input.data(), input.size(), output, error))
    {
        return send_error(error);
    }

    return send_result(output);
}
