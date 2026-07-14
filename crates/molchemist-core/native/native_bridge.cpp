#include "coordgen_engine.h"

#include <algorithm>
#include <cstring>
#include <exception>

namespace
{
void write_error(const std::string& message, char* error, size_t error_capacity)
{
    if (error == nullptr || error_capacity == 0)
    {
        return;
    }

    const size_t length = std::min(message.size(), error_capacity - 1);
    std::memcpy(error, message.data(), length);
    error[length] = '\0';
}
}
extern "C" int molchemist_coordgen_layout(
    const uint8_t* input,
    size_t input_size,
    uint8_t* output,
    size_t output_capacity,
    size_t* output_size,
    char* error,
    size_t error_capacity)
{
    try
    {
        std::vector<uint8_t> result;
        std::string message;
        if (!molchemist::layout_coordinates(input, input_size, result, message))
        {
            write_error(message, error, error_capacity);
            return 1;
        }

        if (result.size() > output_capacity)
        {
            write_error("Coordinate output buffer is too small", error, error_capacity);
            return 2;
        }

        if (!result.empty())
        {
            std::memcpy(output, result.data(), result.size());
        }
        *output_size = result.size();
        return 0;
    }
    catch (const std::exception& exception)
    {
        write_error(exception.what(), error, error_capacity);
        return 3;
    }
    catch (...)
    {
        write_error("CoordgenLibs failed with an unknown native exception", error, error_capacity);
        return 3;
    }
}
