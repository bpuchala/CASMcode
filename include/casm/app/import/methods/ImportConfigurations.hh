#ifndef CASM_import_ImportConfigurations
#define CASM_import_ImportConfigurations

#include "casm/app/import/ImportInterface.hh"

namespace CASM {

/// Import configurations
class ImportConfigurations : public notstd::Cloneable {
  CLONEABLE(ImportConfigurations)
 public:
  std::string desc() const override;

  std::string name() const override;

  void run(PrimClex &primclex, jsonParser const &json_options,
           jsonParser const &cli_options_as_json) const override;
};

}  // namespace CASM

#endif
