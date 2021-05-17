#ifndef CASM_ImportInterface
#define CASM_ImportInterface

#include <string>

#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

class Log;
class PrimClex;
class ParentInputParser;
class jsonParser;

namespace Completer {
class ImportOption;
}

/// Interface for `casm import` methods
///
/// Notes:
/// - See the interfaces implemented in `casm/app/import/methods` for examples
class ImportInterfaceBase : public notstd::Cloneable {
  ABSTRACT_CLONEABLE(ImportInterfaceBase)
 public:
  /// Import method description. What will be printed by `casm import --type
  /// TypeName`.
  virtual std::string desc() const = 0;

  /// Import type name (i.e. "config", "scel")
  virtual std::string name() const = 0;

  /// Import method implementation
  ///
  /// \param primclex PrimClex gives access to project resources
  /// \param json_options JSON input, as by --input or --settings CLI options
  /// \param cli_options_as_json CLI options converted to JSON. Includes:
  ///        --desc as "desc": array of string (type names to print help for)
  ///        --help as "help": bool (print/return help, including list of
  ///        available methods)
  ///        --type as "method": string (name of enumeration method)
  ///        --settings as "settings": string (path to settings JSON file)
  ///        --input as "input": string (a JSON string)
  ///        --properties as "import_properties": bool
  ///        --data as "import_properties": bool
  ///        --copy-structure-files as "copy_structures_files": bool
  ///        --copy-additional-files as "copy_structures_files" and
  ///        "copy_additional_files"
  ///        --batch as "batch": string (name of batch file)
  ///        --pos as "structures": array of string (location of
  ///        structure files)
  ///        --structures as "structures": array of string (location of
  ///        structure files)
  ///        --selection as "selection": fs::path (selection used to import
  ///        or update structure files from training_data)
  ///
  /// Notes:
  /// - It is up to the individual method to determine how to combine and use
  /// `json_options` and
  ///   `cli_options_as_json`, but by convention prefer using CLI input in case
  ///   of collisions.
  virtual void run(PrimClex &primclex, jsonParser const &json_options,
                   jsonParser const &cli_options_as_json) const = 0;
};

/// Convert `casm import` CLI input to JSON
jsonParser &to_json(const Completer::ImportOption &import_opt,
                    jsonParser &json);

/// Combine --input / --settings JSON with CLI options
ParentInputParser make_import_parent_parser(
    Log &log, jsonParser const &json_options,
    jsonParser const &cli_options_as_json);

}  // namespace CASM

#endif
