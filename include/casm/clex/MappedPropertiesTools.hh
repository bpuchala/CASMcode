#ifndef CASM_MappedPropertiesTools
#define CASM_MappedPropertiesTools

#include <set>
#include <string>

namespace CASM {
namespace xtal {
class BasicStructure;
class SimpleStructure;
}  // namespace xtal

struct MappedProperties;
class Permutation;
class PermuteIterator;
class SymOp;

namespace sym {

/// \brief Apply symmetry to MappedProperties, using permutation for mapping
/// site properties
MappedProperties copy_apply(SymOp const &symop, Permutation const &permutation,
                            MappedProperties const &props);

/// \brief Apply symmetry to MappedProperties, using combined permutation for
/// mapping site properties
MappedProperties copy_apply(PermuteIterator const &op,
                            MappedProperties const &props);

}  // namespace sym

/// Construct MappedProperties from a mapped SimpleStructure
MappedProperties make_mapped_properties(
    xtal::SimpleStructure const &mapped_structure,
    xtal::BasicStructure const &prim);

/// Check if all required properties are included in MappedProperties
bool has_all_required_properties(
    MappedProperties const &properties,
    std::vector<std::string> const &required_properties);

}  // namespace CASM

#endif
