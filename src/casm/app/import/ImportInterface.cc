#include "casm/app/import/ImportInterface.hh"

#include "casm/app/io/json_io.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {

/// Convert `casm import` CLI input to JSON
jsonParser &to_json(const Completer::ImportOption &import_opt,
                    jsonParser &json) {
  const auto &vm = import_opt.vm();

  json.put_obj();
  if (vm.count("desc")) {
    json["desc"] = import_opt.desc_vec();  // vector<std::string>
  }
  if (vm.count("help")) {
    json["help"] = static_cast<bool>(vm.count("help"));  // bool
  }
  if (vm.count("type")) {
    json["type"] = import_opt.configtype();  // str
  }
  if (vm.count("settings")) {
    json["settings"] = import_opt.settings_path().string();  // str
  }
  if (vm.count("input")) {
    json["input"] = import_opt.input_str();  // str
  }
  if (vm.count("properties") || vm.count("data")) {
    json["data"]["import_properties"] = true;
  }
  if (vm.count("copy-structure-files")) {
    json["data"]["copy_structure_files"] = true;
  }
  if (vm.count("copy-additional-files")) {
    json["data"]["copy_structure_files"] = true;
    json["data"]["copy_additional_files"] = true;
  }
  if (vm.count("batch")) {
    json["batch"] = import_opt.batch_path().string();  // fs::path
  }

  if (vm.count("pos") || vm.count("structures")) {
    json["structures"] = jsonParser::array();
  }
  if (vm.count("pos")) {
    for (fs::path const &path : import_opt.pos_vec()) {
      json["structures"].push_back(path.string());
    }
  }
  if (vm.count("structures")) {
    for (fs::path const &path : import_opt.structures_vec()) {
      json["structures"].push_back(path.string());
    }
  }

  if (vm.count("selection")) {
    json["selection"] = import_opt.selection_path().string();  // std::string
  }
  return json;
}

/// Combine --input / --settings JSON with CLI options
ParentInputParser make_import_parent_parser(
    Log &log, jsonParser const &json_options,
    jsonParser const &cli_options_as_json) {
  log.indent() << "Input from JSON (--input or --setings):\n"
               << json_options << std::endl
               << std::endl;
  log.indent() << "Input from `casm import` options:\n"
               << cli_options_as_json << std::endl
               << std::endl;

  // combine JSON options and CLI options
  std::map<std::string, std::string> cli_to_combined_keys{
      {"selection", "selection"},    // --selection
      {"structures", "structures"},  // --structures, --pos
      {"batch", "batch"},            // --batch
      {"data", "data"}               // --data, --properties,
                                     // --copy-structure-files,
                                     // --copy-additional-files
  };

  jsonParser json_combined{json_options};
  combine_json_options(cli_to_combined_keys, cli_options_as_json,
                       json_combined);

  log.indent() << "Combined Input:\n"
               << json_combined << std::endl
               << std::endl;
  return ParentInputParser{json_combined};
}

}  // namespace CASM
