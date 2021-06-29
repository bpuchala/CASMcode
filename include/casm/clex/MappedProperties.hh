#ifndef CASM_MappedProperties
#define CASM_MappedProperties

#include <ctime>
#include <string>

#include "casm/casm_io/FileData.hh"
#include "casm/casm_io/enum/io_traits.hh"
#include "casm/casm_io/enum/stream_io.hh"
#include "casm/global/eigen.hh"

namespace CASM {

/// Calculated properties of a Configuration
///
/// Associated with a particular Configuration as "calculated properties"
/// (std::map<std::string, MappedProperties>), MappedProperties are values of
/// properties (i.e. energy, relaxation displacements, relaxation strain, etc.)
/// that are dependent on the values of the DoF.
///
/// MappedProperties are "mapped" in the sense that they are not necessarily the
/// raw output values from a DFT calculation (though some properties are), but
/// are determined by "mapping" the final calculated structure's atomic
/// positions and lattice vectors to the ideal configuration they are most
/// similar to. In the case of significant atomic and lattice relaxations this
/// allows associating calculated properties with the correct Configuration. If
/// the relaxations are too severe, there may be no Configuration that is an
/// appropriate mapping. The user has control of criteria used to judge whether
/// any Configuration is an appropriate mapping and which is the best mapping.
///
/// Note: Multiple instances of MappedProperties may be associated with a single
/// Configuration. This may occur if multiple calculation types are performed
/// for a single Configuration. This may also occur if multiple ideal
/// Configuration are unstable and relax to the same other Configuration.
///
/// MappedProperties has:
/// - "global" (std::map<std::string, Eigen::MatrixXd>), a map of property name
/// to value for
///   global properties properties of the configuration such as energy, strain
///   metrics, or lattice vectors. Property names must follow CASM property
///   naming conventions as documented for AnisoValTraits. This map also holds
///   any scalar properties such as energy or mapping cost
///   ("lattice_deformation_cost", "atomic_deformation_cost", "total_cost"),
///   which can be queried and set using member functions. Properties are
///   expected to be vectors expressed in the \link DoFSet "standard
///   basis"\endlink.
/// - "site" (std::map<std::string, Eigen::MatrixXd>), a map of property name to
/// value for
///   properties of particular sites such as displacements, forces, or atomic
///   coordinates. Property names must follow CASM property naming conventions
///   as documented for AnisoValTraits. Site properties are expected to be
///   vectors (columns of the matrices) expressed in the \link DoFSet "standard
///   basis"\endlink.
///
/// MappedProperties also has:
/// - "origin" (std::string), indicates the source of the properties. This can
/// vary depending on
///   the context, but most typically it is a path to a `properties.calc.json`
///   type file. In some cases "prim:" is prepended to the file path to indicate
///   that the "to" configuration is the primitive of the structure read from
///   the origin file.
/// - "to" (std::string), the name of the Configuration that the properties are
/// mapped to. It must
///   be a Configuration name.
/// - "init_config" (std::string), used to indicate the Configuration that was
/// calculated to
///   generate these properties. Typically, it is the name of the Configuration
///   which was used to generate the initial structure input to the calculation.
///   The string "import" is typically used when the properties were imported
///   from outside the project and no initial ideal Configuration is known.
/// - "file_data" (FileData), used to store the last write time of the "origin"
/// file
///
/// Note: In the future, "origin", "to", "init_config", and "file_data" may be
/// moved out of MappedProperties and stored as part of the PropertiesDatabase
/// or as needed by `casm update` and `casm import`.
///
struct MappedProperties {
  /// \brief default construction sets init_config field to "import"
  MappedProperties() : init_config("import") {}

  /// Map of global property name to property value.
  /// Property names must follow CASM property naming conventions as documented
  /// for AnisoValTraits.
  ///
  /// Properties are expected to be vectors expressed in the \link xtal::DoFSet
  /// "standard basis"\endlink.
  std::map<std::string, Eigen::MatrixXd> global;

  /// Map of site property name to property value.
  /// Property names must follow CASM property naming conventions as documented
  /// for AnisoValTraits.
  ///
  /// Properties are expected to be vectors (columns of the matrices) expressed
  /// in the \link xtal::DoFSet "standard basis"\endlink.
  std::map<std::string, Eigen::MatrixXd> site;

  bool has_scalar(std::string const &_name) const;

  double const &scalar(std::string const &_name) const;

  double &scalar(std::string const &_name);

  // Note: In the future, the following may be moved out of MappedProperties

  /// Indicates the origin of the properties
  ///
  /// Typical usage: The path to the `properties.calc.json` type file containing
  /// the SimpleStructure with these properties.
  std::string origin;

  /// The name of the configuration that these properties are mapped to (the
  /// ending point of relaxation)
  std::string to;

  /// The name of the ideal configuration that served as the starting point for
  /// calculation
  ///
  /// Typical usage:
  /// - If generated by 'casm update', the name of the starting config.
  /// - If generated by 'casm import', "import"
  std::string init_config;

  /// Path to properties file and time of last write time corresponding to
  /// 'origin'
  FileData file_data;
};

jsonParser &to_json(const MappedProperties &prop, jsonParser &json);

jsonParser const &from_json(MappedProperties &prop, const jsonParser &json);

/// Resolve mapping conflicts by 'scoring' the MappedProperties structure
///
/// Options:
/// \code
/// {
///   "method" : "deformation_cost",
///   "lattice_weight" : number in range [0,1.0]
/// }
/// {
///   "method" : "minimum",
///   "property" : property name (i.e. "energy")
/// }
/// {
///   "method" : "maximum",
///   "property" : property name (i.e. "some property")
/// }
/// {
///   "method" : "direct_selection",
///   "name" : configname to force as 'best' (i.e. "SCEL3_1_1_3_0_0_0/4")
/// }
/// \endcode
///
class ScoreMappedProperties {
 public:
  enum class Method { deformation_cost, minimum, maximum, direct_selection };

  static std::string method_name(Method m) {
    if (m == Method::minimum)
      return "minimum";
    else if (m == Method::maximum)
      return "maximum";
    else if (m == Method::direct_selection)
      return "direct_selection";
    else if (m == Method::deformation_cost)
      return "deformation_cost";
    return "unknown";
  }

  struct Option {
    Option(Method _method = Method::minimum, std::string _name = "energy")
        : Option(_method, _name, -1.) {}

    Option(Method _method, double _lattice_weight = 0.5)
        : Option(_method, "", _lattice_weight) {}

    /// Method for scoring
    Method method;

    /// Property name or configname used for scoring
    std::string name;

    double lattice_weight;

   private:
    Option(Method _method, std::string _name, double _lattice_weight);
  };

  /// \brief Default uses minimum energy
  explicit ScoreMappedProperties(Option _opt = Option(Method::minimum,
                                                      "energy"));

  double operator()(const MappedProperties &obj) const;

  bool validate(const MappedProperties &obj) const;

  bool operator==(const ScoreMappedProperties &B) const;

  bool operator!=(const ScoreMappedProperties &B) const;

  const Option &option() const { return m_opt; }

 private:
  Option m_opt;
};

jsonParser &to_json(const ScoreMappedProperties &score, jsonParser &json);

jsonParser const &from_json(ScoreMappedProperties &score,
                            const jsonParser &json);
}  // namespace CASM

namespace CASM {
ENUM_IO_DECL(CASM::ScoreMappedProperties::Method)
ENUM_TRAITS(CASM::ScoreMappedProperties::Method)
}  // namespace CASM

#endif
