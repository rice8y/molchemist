#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace molchemist
{
bool layout_coordinates(
    const uint8_t* input_data,
    size_t input_size,
    std::vector<uint8_t>& output,
    std::string& error);
}
