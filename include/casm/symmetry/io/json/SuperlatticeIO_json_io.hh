#ifndef CASM_symmetry_io_json_SuperlatticeIO
#define CASM_symmetry_io_json_SuperlatticeIO

#include <memory>
#include "casm/crystallography/Superlattice.hh"

namespace CASM {

template <typename T>
class InputParser;
class Structure;

/// Used for InputParser parsing of Lattice using one of several possible input
/// specifications, including supercell name
struct SuperlatticeIO {
  SuperlatticeIO(xtal::Superlattice const &_superlattice)
      : superlattice(_superlattice) {}
  xtal::Superlattice superlattice;
};

/// \brief Parse a Superlattice from JSON, allowing several possible input
/// specifications, including supercell name
void parse(InputParser<SuperlatticeIO> &parser,
           std::shared_ptr<Structure const> shared_prim);

}  // namespace CASM

#endif
