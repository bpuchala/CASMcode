#ifndef CASM_import_ImportConfigurations
#define CASM_import_ImportConfigurations

#include "casm/app/import/methods/ImportConfigurations.hh"

namespace CASM {

std::string ImportConfigurations::desc() const {
  return "Import configurations options ...";
}

std::string ImportConfigurations::name() const {
  return traits<Configuration>::short_name;
}

void ImportConfigurations::run(PrimClex &primclex,
                               jsonParser const &json_options,
                               jsonParser const &cli_options_as_json) const {
  std::cout << "run `import --type config`" << std::endl;
}

}  // namespace CASM

#endif
