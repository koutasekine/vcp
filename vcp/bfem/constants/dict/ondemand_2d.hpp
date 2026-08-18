// vcp/bfem/constants/dict/ondemand_2d.hpp
//
// CONST-B2b: 2D ondemand ledger (dictionary layer separation, ruling R21).
// SHIPPED EMPTY: this file is the user-environment growth surface for
// exact-key entries produced by ondemand_l2_projection_constant_sq
// (dict/resolve.hpp).  The LIBRARY NEVER WRITES THIS FILE; entries are
// inserted only by the explicit probe sandbox/probes/
// constb2b_ondemand_append.cpp (or by hand), always above the
// VCP_B2B_ONDEMAND_APPEND_POINT anchor inside the marker block.  Locally
// grown revisions stay in the user environment and are not committed to
// the canonical git history (shipped-layer revisions go through a
// generation track instead).
//
// Entry semantics are IDENTICAL to dict/registry_2d.hpp: exact canonical
// similarity key as integer fractions, Chat^2 = C_d^2 / h^2 as a
// 12-significant-digit upward decimal string; consumption
// C_d(K)^2 <= Chat^2 * h^2(K).  resolve.hpp reads this table on an equal
// footing with the registry (exact key match, source = registry_exact).
//
// Lexical regime B: decimal strings are authorized INSIDE the marker
// block only (cm1_ch_tests G7); shipped in-table count is 0.

#ifndef VCP_BFEM_CONSTANTS_DICT_ONDEMAND_2D_HPP
#define VCP_BFEM_CONSTANTS_DICT_ONDEMAND_2D_HPP

#include <vcp/bfem/constants/dict/dict_entry.hpp>

namespace vcp {
namespace bfem {
namespace constants {
namespace detail {

// VCP_CONSTANTS_TABLE_BEGIN  (ondemand ledger; grows in the user environment)

inline const l2_projection_registry_entry_2d*
l2_projection_ondemand_2d_entries(int& count) {
    static const l2_projection_registry_entry_2d entries[] = {
        // VCP_B2B_ONDEMAND_APPEND_POINT (insert entries ABOVE this line)
        { { 0LL, 0LL, 0LL }, { 1LL, 1LL, 1LL }, -1, "sentinel" }
    };
    count = static_cast<int>(sizeof entries / sizeof entries[0]) - 1;
    return entries;
}

// VCP_CONSTANTS_TABLE_END

} // namespace detail

inline const l2_projection_registry_entry_2d*
l2_projection_ondemand_2d_table(int& count) {
    return detail::l2_projection_ondemand_2d_entries(count);
}

} // namespace constants
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_CONSTANTS_DICT_ONDEMAND_2D_HPP
