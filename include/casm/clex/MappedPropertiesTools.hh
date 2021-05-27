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
class PermuteIterator;

///
MappedProperties copy_apply(PermuteIterator const &op,
                            MappedProperties const &props);

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
