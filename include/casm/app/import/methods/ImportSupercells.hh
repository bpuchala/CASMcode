#ifndef CASM_import_ImportSupercells
#define CASM_import_ImportSupercells

#include "casm/app/import/ImportInterface.hh"

namespace CASM {

/// Import supercells
class ImportSupercells : public notstd::Cloneable {
  CLONEABLE(ImportSupercells)
 public:
  std::string desc() const override;

  std::string name() const override;

  void run(PrimClex &primclex, jsonParser const &json_options,
           jsonParser const &cli_options_as_json) const override;
};

}  // namespace CASM

#endif
