#ifndef MOLECULE_HH
#define MOLECULE_HH

#include <iostream>
#include <array>

#include "casm/casm_io/json_io/container.hh"
#include "casm/basis_set/DoF.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM {

  /** \defgroup Molecule
   *  \ingroup Crystallography
   *  \brief Relates to Molecule
   *  @{
   */

  class Molecule;
  template <typename T>
  class OccupantDoF;
  typedef OccupantDoF<Molecule> MoleculeOccupant;
  //****************************************************

  /// \brief A vacancy is any Specie/Molecule with (name == "VA" || name == "va" || name == "Va")
  inline bool is_vacancy(const std::string &name) {
    return (name == "VA" || name == "va" || name == "Va");
  }

  class Specie : public Comparisons<Specie> {
  public:
    std::string name;
    //BP: for current Kinetics logic, Specie is essentially "Atom"
    //  (spherically symmetric, comparable by name only)
    //  double mass, magmom, U, J; //Changed 05/10/13 -- was int
    Specie() { };
    explicit Specie(std::string init_name) : name(init_name) { };

    bool is_vacancy() const {
      return CASM::is_vacancy(name);
    }

    bool operator<(const Specie &RHS) const {
      return name < RHS.name;
    };

  private:

    friend Comparisons<Specie>;
    bool _eq(const Specie &RHS) const {
      return (name == RHS.name);
    };

  };

  jsonParser &to_json(const Specie &specie, jsonParser &json);

  void from_json(Specie &specie, const jsonParser &json);

  //****************************************************

  class AtomPosition : public Coordinate {
  public:
    typedef std::array<bool, 3> sd_type;
    Specie specie;

    sd_type SD_flag;

    explicit AtomPosition(const Lattice &init_lattice) : Coordinate(0, 0, 0, init_lattice, CART) { };
    AtomPosition(double elem1,
                 double elem2,
                 double elem3,
                 std::string sp_name,
                 const Lattice &init_lattice,
                 COORD_TYPE mode,
    sd_type _SD_flag = sd_type {{false, false, false}}) :
      Coordinate(elem1, elem2, elem3, init_lattice, mode),
      specie(sp_name),
      SD_flag(_SD_flag) { };


    bool operator==(const AtomPosition &RHS) const;

    void print(std::ostream &stream, const Coordinate &trans, int spaces, bool SD_is_on = false) const;
    // TODO: If comparing coordinates alone does not suffice, add a == operator here.

    AtomPosition &apply_sym(const SymOp &op);
    AtomPosition &apply_sym_no_trans(const SymOp &op);


  };

  jsonParser &to_json(const AtomPosition &apos, jsonParser &json);

  // Lattice must be set already
  void from_json(AtomPosition &apos, const jsonParser &json);

  //****************************************************

  /** \defgroup Molecule
   *  \ingroup Crystallography
   *  \brief Relates to Molecule
   *  @{
   */

  class Molecule :
    public Array<AtomPosition>,
    public Comparisons<Molecule> {

    Lattice const *m_home;

    Array<SymOp> point_group;

  public:
    using Array<AtomPosition>::push_back;
    using Array<AtomPosition>::at;
    using Array<AtomPosition>::size;

    Coordinate center;
    std::string name;

    explicit Molecule(const Lattice &init_home) : m_home(&init_home), center(0, 0, 0, init_home, CART) {};

    Lattice const *home() const {
      return m_home;
    }

    bool is_vacancy() const {
      return CASM::is_vacancy(name);
    };

    /// \brief Check if Molecule is indivisible
    ///
    /// - Currently, always false
    bool is_indivisible() const {
      return false;
    }

    Molecule &apply_sym(const SymOp &op);
    Molecule &apply_sym_no_trans(const SymOp &op);

    void set_lattice(const Lattice &new_lat, COORD_TYPE invariant_mode);

    bool contains(const std::string &name) const;

    void read(std::istream &stream);
    void print(std::ostream &stream, const Coordinate &trans, int spaces, char delim, bool SD_is_on = false) const;

    jsonParser &to_json(jsonParser &json) const;

    // Lattice must be set already
    void from_json(const jsonParser &json);

    /// \brief Name comparison via '<', '>', '<=', '>='
    bool operator<(const Molecule &B) const;

    using Comparisons<Molecule>::operator==;
    using Comparisons<Molecule>::operator!=;
    using Comparisons<Molecule>::operator<=;
    using Comparisons<Molecule>::operator>;
    using Comparisons<Molecule>::operator>=;

  private:

    friend Comparisons<Molecule>;

    /// \brief center and AtomPosition comparison via '==', '!='
    bool _eq(const Molecule &B) const;
  };

  /// \brief Return an atomic Molecule with specified name and Lattice
  Molecule make_atom(std::string atom_name, const Lattice &lat);

  /// \brief Return an vacancy Molecule with specified Lattice
  Molecule make_vacancy(const Lattice &lat);

  /// \brief Return true if Molecule name matches 'name', including Va checks
  bool is_molecule_name(const Molecule &mol, std::string name);


  jsonParser &to_json(const Molecule &mol, jsonParser &json);

  // Lattice must be set already
  void from_json(Molecule &mol, const jsonParser &json);

  /** @} */
};
#endif
