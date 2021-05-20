#ifndef CASM_crystallography_StrucMapping_json_io
#define CASM_crystallography_StrucMapping_json_io

namespace CASM {

namespace xtal {
struct MappingNode;
}

class jsonParser;

jsonParser &to_json(xtal::MappingNode const &mapping_node, jsonParser &json);

}  // namespace CASM

#endif
