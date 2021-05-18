#include "casm/app/import/methods/ImportSupercells.hh"

#include "casm/clex/Supercell.hh"
#include "casm/global/definitions.hh"

namespace CASM {

std::string ImportSupercells::desc() const {
  return "Import supercells options ...";
}

std::string ImportSupercells::name() const {
  return traits<Supercell>::short_name;
}

void ImportSupercells::run(PrimClex &primclex, jsonParser const &json_options,
                           jsonParser const &cli_options_as_json) const {
  std::cout << "run `import --type scel`" << std::endl;
}

}  // namespace CASM
