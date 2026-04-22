#include "protocol.h"

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "sketcherMinimizer.h"

extern "C" void _initialize();

namespace
{
    constexpr unsigned char kLayoutMagic[4] = {'M', 'C', 'G', '2'};
    constexpr unsigned char kCoordMagic[4] = {'M', 'C', 'C', '2'};
    constexpr uint32_t kMissingIndex = 0xFFFFFFFFu;

    struct AtomStereoRecord
    {
        uint32_t atom = 0;
        uint32_t looking_from = kMissingIndex;
        uint32_t atom1 = kMissingIndex;
        uint32_t atom2 = kMissingIndex;
        uint8_t direction = 0;
    };

    struct DoubleBondStereoRecord
    {
        uint32_t bond = 0;
        uint32_t atom1 = 0;
        uint32_t atom2 = 0;
        uint8_t is_z = 0;
    };

    bool read_u32_le(const std::vector<uint8_t>& input, size_t& offset, uint32_t& value)
    {
        if (offset + 4 > input.size())
        {
            return false;
        }

        value = static_cast<uint32_t>(input[offset]) |
                (static_cast<uint32_t>(input[offset + 1]) << 8) |
                (static_cast<uint32_t>(input[offset + 2]) << 16) |
                (static_cast<uint32_t>(input[offset + 3]) << 24);
        offset += 4;
        return true;
    }

    void write_u32_le(std::vector<uint8_t>& output, uint32_t value)
    {
        output.push_back(static_cast<uint8_t>(value & 0xFF));
        output.push_back(static_cast<uint8_t>((value >> 8) & 0xFF));
        output.push_back(static_cast<uint8_t>((value >> 16) & 0xFF));
        output.push_back(static_cast<uint8_t>((value >> 24) & 0xFF));
    }

    void write_f32_le(std::vector<uint8_t>& output, float value)
    {
        uint32_t bits = 0;
        static_assert(sizeof(bits) == sizeof(value), "f32 size mismatch");
        std::memcpy(&bits, &value, sizeof(value));
        write_u32_le(output, bits);
    }

    int send_result(const std::vector<uint8_t>& result)
    {
        wasm_minimal_protocol_send_result_to_host(result.data(), result.size());
        return 0;
    }

    int send_error(const char* message)
    {
        wasm_minimal_protocol_send_result_to_host(
            reinterpret_cast<const uint8_t*>(message),
            std::strlen(message));
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

    if (input.size() < 20 || std::memcmp(input.data(), kLayoutMagic, 4) != 0)
    {
        return send_error("Invalid layout payload");
    }

    size_t offset = 4;
    uint32_t atom_count = 0;
    uint32_t bond_count = 0;
    uint32_t atom_stereo_count = 0;
    uint32_t double_bond_stereo_count = 0;
    if (!read_u32_le(input, offset, atom_count) ||
        !read_u32_le(input, offset, bond_count) ||
        !read_u32_le(input, offset, atom_stereo_count) ||
        !read_u32_le(input, offset, double_bond_stereo_count))
    {
        return send_error("Invalid layout payload header");
    }

    const size_t expected_size =
        20 +
        static_cast<size_t>(atom_count) +
        static_cast<size_t>(bond_count) * 9 +
        static_cast<size_t>(atom_stereo_count) * 17 +
        static_cast<size_t>(double_bond_stereo_count) * 13;
    if (input.size() != expected_size)
    {
        return send_error("Layout payload length does not match atom/bond counts");
    }

    auto* molecule = new sketcherMinimizerMolecule();
    std::vector<sketcherMinimizerAtom*> atoms;
    atoms.reserve(atom_count);
    std::vector<sketcherMinimizerBond*> bonds;
    bonds.reserve(bond_count);

    for (uint32_t i = 0; i < atom_count; ++i)
    {
        const auto atomic_number = input[offset++];
        auto* atom = molecule->addNewAtom();
        atom->setAtomicNumber(atomic_number);
        atoms.push_back(atom);
    }

    for (uint32_t i = 0; i < bond_count; ++i)
    {
        uint32_t atom1 = 0;
        uint32_t atom2 = 0;
        if (!read_u32_le(input, offset, atom1) || !read_u32_le(input, offset, atom2))
        {
            delete molecule;
            return send_error("Bond payload ended unexpectedly");
        }

        const auto order = input[offset++];
        if (atom1 >= atom_count || atom2 >= atom_count || order == 0 || order > 3)
        {
            delete molecule;
            return send_error("Layout payload contains an invalid bond");
        }

        auto* bond = molecule->addNewBond(atoms[atom1], atoms[atom2]);
        bond->setBondOrder(static_cast<int>(order));
        bonds.push_back(bond);
    }

    std::vector<AtomStereoRecord> atom_stereo;
    atom_stereo.reserve(atom_stereo_count);
    for (uint32_t i = 0; i < atom_stereo_count; ++i)
    {
        AtomStereoRecord record;
        if (!read_u32_le(input, offset, record.atom) ||
            !read_u32_le(input, offset, record.looking_from) ||
            !read_u32_le(input, offset, record.atom1) ||
            !read_u32_le(input, offset, record.atom2))
        {
            delete molecule;
            return send_error("Atom stereo payload ended unexpectedly");
        }
        record.direction = input[offset++];
        atom_stereo.push_back(record);
    }

    std::vector<DoubleBondStereoRecord> double_bond_stereo;
    double_bond_stereo.reserve(double_bond_stereo_count);
    for (uint32_t i = 0; i < double_bond_stereo_count; ++i)
    {
        DoubleBondStereoRecord record;
        if (!read_u32_le(input, offset, record.bond) ||
            !read_u32_le(input, offset, record.atom1) ||
            !read_u32_le(input, offset, record.atom2))
        {
            delete molecule;
            return send_error("Double-bond stereo payload ended unexpectedly");
        }
        record.is_z = input[offset++];
        double_bond_stereo.push_back(record);
    }

    sketcherMinimizerMolecule::assignBondsAndNeighbors(
        molecule->getAtoms(), molecule->getBonds());

    for (const auto& record : atom_stereo)
    {
        const bool bad_neighbor_index =
            (record.looking_from != kMissingIndex && record.looking_from >= atom_count) ||
            (record.atom1 != kMissingIndex && record.atom1 >= atom_count) ||
            (record.atom2 != kMissingIndex && record.atom2 >= atom_count);

        if (record.atom >= atom_count || bad_neighbor_index || (record.direction != 1 && record.direction != 2))
        {
            delete molecule;
            return send_error("Atom stereo payload contains an invalid atom index");
        }

        sketcherMinimizerAtomChiralityInfo info;
        info.lookingFrom = record.looking_from == kMissingIndex ? nullptr : atoms[record.looking_from];
        info.atom1 = record.atom1 == kMissingIndex ? nullptr : atoms[record.atom1];
        info.atom2 = record.atom2 == kMissingIndex ? nullptr : atoms[record.atom2];
        info.direction = record.direction == 1
            ? sketcherMinimizerAtomChiralityInfo::clockwise
            : sketcherMinimizerAtomChiralityInfo::counterClockwise;

        atoms[record.atom]->setStereoChemistry(info);
        atoms[record.atom]->setAbsoluteStereoFromChiralityInfo();
    }

    for (const auto& record : double_bond_stereo)
    {
        if (record.bond >= bond_count || record.atom1 >= atom_count || record.atom2 >= atom_count)
        {
            delete molecule;
            return send_error("Double-bond stereo payload contains an invalid index");
        }

        sketcherMinimizerBondStereoInfo info;
        info.atom1 = atoms[record.atom1];
        info.atom2 = atoms[record.atom2];
        info.stereo = record.is_z
            ? sketcherMinimizerBondStereoInfo::cis
            : sketcherMinimizerBondStereoInfo::trans;

        auto* bond = bonds[record.bond];
        bond->isZEActive = true;
        bond->setStereoChemistry(info);
        bond->setAbsoluteStereoFromStereoInfo();
    }

    sketcherMinimizer minimizer;
    minimizer.initialize(molecule);
    minimizer.runGenerateCoordinates();

    std::vector<uint8_t> output;
    output.reserve(12 + static_cast<size_t>(atom_count) * 8 + bond_count);
    output.insert(output.end(), std::begin(kCoordMagic), std::end(kCoordMagic));
    write_u32_le(output, atom_count);
    write_u32_le(output, bond_count);

    for (auto* atom : atoms)
    {
        const auto point = atom->getCoordinates();
        write_f32_le(output, point.x());
        write_f32_le(output, point.y());
    }

    for (auto* bond : bonds)
    {
        uint8_t style = 0;
        if (bond->hasStereochemistryDisplay)
        {
            if (bond->isWedge)
            {
                style = bond->isReversed ? 2 : 1;
            }
            else
            {
                style = bond->isReversed ? 4 : 3;
            }
        }
        output.push_back(style);
    }

    return send_result(output);
}
