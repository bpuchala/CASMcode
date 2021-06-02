#include "casm/symmetry/io/json/SuperlatticeIO_json_io.hh"

#include <vector>

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/crystallography/CanonicalForm.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/symmetry/SupercellSymInfo.hh"

namespace CASM {

/// \brief Parse a Superlattice from JSON, allowing several possible input
/// specifications, including supercell name
///
/// This parser expects a JSON object, from which it will parse a Superlattice
/// from any of the following attributes, using first found when searching in
/// the order listed below. If none found, the prim lattice is used.
///
/// Input options:
///
///   transformation_matrix_to_super: 3x3 array of integer (optional)
///     Transformation matrix T, defining the supercell lattice vectors
///     S, in terms of the prim lattice vectors, P: `S = P * T`, where
///     S and P are column vector matrices.
///
///   supercell_lattice_row_vectors: 3x3 array of integer (optional)
///     Supercell lattice vectors, as a row vector matrix.
///
///   supercell_lattice_column_matrix: 3x3 array of integer (optional)
///     Supercell lattice vectors, as a column vector matrix.
///
///   supercell_name: string (optional)
///     A name given to all equivalent super lattices of the prim
///     lattice. For the canonical super lattice, the name is
///     constructed from the hermite normal form of
///     `transformation_matrix_to_super`. For a non-canonical super
///     lattice, the name is the constructed from the name of the
///     canonical super lattice and the index of the prim factor group
///     operation that (excluding the shift) transforms the canonical
///     super lattice into this super lattice.
///
///     Example 1: Canonical supercell name
///
///         \"SCEL8_4_2_1_1_3_2\"
///
///     Example 2: Non-canonical supercell name representing a
///     re-orientation by application of prim factor group operation
///     with index 2 (indexing starting at 0).
///
///         \"SCEL8_4_2_1_1_3_2.2\"
///
///   make_canonical: bool (optional, default=false)
///     If \"true\", the canonical equivalent supercell is used.
///
void parse(InputParser<SuperlatticeIO> &parser,
           std::shared_ptr<Structure const> shared_prim) {
  Log &log = CASM::log();

  // read "transformation_matrix_to_super"
  Eigen::Matrix3l T;
  if (parser.self.contains("transformation_matrix_to_super")) {
    parser.optional(T, "transformation_matrix_to_super");

    // or read "supercell_lattice_row_vectors"
  } else if (parser.self.contains("supercell_lattice_row_vectors")) {
    Eigen::Matrix3d L_transpose;
    parser.optional(L_transpose, "supercell_lattice_row_vectors");
    Lattice super_lattice{L_transpose.transpose()};
    try {
      T = make_transformation_matrix_to_super(shared_prim->lattice(),
                                              super_lattice, TOL);
    } catch (std::exception &e) {
      auto pair = is_superlattice(super_lattice, shared_prim->lattice(), TOL);
      log << "The transformation_matrix_to_super determined from "
             "\"supercell_lattice_row_vectors\" is not approximately integer. "
             "Found:\n"
          << pair.second << std::endl;
      parser.insert_error("supercell_lattice_row_vectors", e.what());
    }

    // or read "supercell_lattice_column_matrix"
  } else if (parser.self.contains("supercell_lattice_column_matrix")) {
    Eigen::Matrix3d L;
    parser.optional(L, "supercell_lattice_column_matrix");
    Lattice super_lattice{L};
    try {
      T = make_transformation_matrix_to_super(shared_prim->lattice(),
                                              super_lattice, TOL);
    } catch (std::exception &e) {
      auto pair = is_superlattice(super_lattice, shared_prim->lattice(), TOL);
      log << "The transformation_matrix_to_super determined from "
             "\"supercell_lattice_column_matrix\" is not approximately "
             "integer. Found:\n"
          << pair.second << std::endl;
      parser.insert_error("supercell_lattice_column_matrix", e.what());
    }

    // or read "supercell_name"
  } else if (parser.self.contains("supercell_name")) {
    std::string supercell_name;
    parser.optional(supercell_name, "supercell_name");
    try {
      xtal::Superlattice superlattice = make_superlattice_from_supercell_name(
          shared_prim->factor_group(), shared_prim->lattice(), supercell_name);
      T = superlattice.transformation_matrix_to_super();
    } catch (std::exception &e) {
      parser.insert_error("supercell_name", e.what());
    }

    // else use Identity (prim cell)
  } else {
    T = Eigen::Matrix3l::Identity();
  }

  // read "make_canonical"
  bool make_canonical = false;
  parser.optional(make_canonical, "make_canonical");

  if (!parser.valid()) {
    return;
  }

  Lattice super_lattice = make_superlattice(shared_prim->lattice(), T);
  if (make_canonical) {
    super_lattice =
        xtal::canonical::equivalent(super_lattice, shared_prim->point_group());
  }
  parser.value = notstd::make_unique<SuperlatticeIO>(
      xtal::Superlattice(shared_prim->lattice(), super_lattice));
}

} // namespace CASM
