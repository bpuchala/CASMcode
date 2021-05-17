#ifndef CASM_clex_io_json_ConfigMapping
#define CASM_clex_io_json_ConfigMapping

#include <memory>

#include "casm/casm_io/dataformatter/DataFormatterDecl.hh"

namespace CASM {
namespace ConfigMapping {
struct Settings;
}
class Structure;
class Supercell;
class jsonParser;

// jsonParser &to_json(ConfigMapping::Settings const &_set, jsonParser &_json);

/// Read ConfigMapping::Settings from JSON
jsonParser const &from_json(
    ConfigMapping::Settings &_set, jsonParser const &_json,
    std::shared_ptr<Structure const> const &shared_prim,
    DataFormatterDictionary<Supercell> const &supercell_query_dict);

}  // namespace CASM

#endif
