#ifndef UNITCELLCOORD_HH
#define UNITCELLCOORD_HH

#include <iostream>

#include "casm/CASM_global_definitions.hh"
#include "casm/CASM_global_Eigen.hh"
#include "casm/misc/Comparisons.hh"

namespace CASM {

  class jsonParser;
  class Coordinate;
  class Site;
  class Structure;
  class SymOp;

  /** \ingroup Coordinate
   *  @{
   */

  /// \brief Unit Cell Indices
  ///
  /// - Integer vector to represent a particular unit cell using multiples of the unit cell vectors
  /// - Is a <a href="http://eigen.tuxfamily.org/dox/group__QuickRefPage.html">Eigen::Vector3i</a>.
  ///
  class UnitCell : public Eigen::Vector3l {
  public:

    UnitCell(void) : Eigen::Vector3l() {}

    UnitCell(Index a, Index b, Index c) : Eigen::Vector3l(a, b, c) {}

    // This constructor allows you to construct MyVectorType from Eigen expressions
    template<typename OtherDerived>
    UnitCell(const Eigen::MatrixBase<OtherDerived> &other) :
      Eigen::Vector3l(other) {}

    // This method allows you to assign Eigen expressions to MyVectorType
    template<typename OtherDerived>
    UnitCell &operator=(const Eigen::MatrixBase <OtherDerived> &other) {
      this->Eigen::Vector3l::operator=(other);
      return *this;
    }
  };

  template<typename Derived>
  struct Translatable {

    Derived operator+(UnitCell frac) {
      Derived tmp {derived()};
      return tmp += frac;
    }

    Derived operator-(UnitCell frac) {
      Derived tmp {derived()};
      return tmp -= frac;
    }

  protected:
    const Derived &derived() const {
      return *static_cast<const Derived *>(this);
    }
  };

  /* -- UnitCellCoord Declarations ------------------------------------- */

  /// \brief Unit Cell Coordinates
  ///
  /// - Represent a crystal site using UnitCell indices and sublattice index
  ///
  class UnitCellCoord : public Comparisons<UnitCellCoord>, public Translatable<UnitCellCoord> {

  public:

    /// \brief SymBasisPermute rep should be obtainable from UnitType
    typedef Structure UnitType;

    UnitCellCoord(const UnitType &unit);

    UnitCellCoord(const UnitType &unit, Index _sublat, const UnitCell &_unitcell);

    UnitCellCoord(const UnitType &unit, Index _sublat, Index i, Index j, Index k);

    UnitCellCoord(const UnitType &unit, const Coordinate &coord, double tol);


    UnitCellCoord(const UnitCellCoord &B) = default;

    UnitCellCoord &operator=(const UnitCellCoord &B) = default;

    UnitCellCoord(UnitCellCoord &&B) = default;

    UnitCellCoord &operator=(UnitCellCoord &&B) = default;


    /// \brief Get unit structure reference
    const UnitType &unit() const;

    /// \brief Change unit structure, keeping indices constant
    void set_unit(const UnitType &_unit);

    /// \brief Access the Lattice
    const Lattice &lattice() const;

    /// \brief Get corresponding coordinate
    Coordinate coordinate() const;

    /// \brief Get a copy of corresponding site
    Site site() const;

    /// \brief Get reference to corresponding sublattice site in the unit structure
    const Site &sublat_site() const;

    UnitCell &unitcell();
    const UnitCell &unitcell() const;

    Index &unitcell(Index i);
    const Index &unitcell(Index i) const;

    Index &sublat();
    const Index &sublat() const;

    Index &operator[](Index i);
    const Index &operator[](Index i) const;

    UnitCellCoord &operator+=(UnitCell frac);

    UnitCellCoord &operator-=(UnitCell frac);

    bool operator<(const UnitCellCoord &B) const;

    UnitCellCoord &apply_sym(const SymOp &op);

    UnitCellCoord copy_apply(const SymOp &op) const;

    operator Coordinate() const;

  private:

    /// make _eq accessible
    friend class Comparisons<UnitCellCoord>;

    bool _eq(const UnitCellCoord &B) const;

    const UnitType *m_unit;
    UnitCell m_unitcell;
    Index m_sublat;

  };

  inline std::ostream &operator<<(std::ostream &sout, const UnitCellCoord &site) {
    return sout << site.sublat() << ", " << site.unitcell().transpose();
  }

  /// \brief Print to json as [b, i, j, k]
  jsonParser &to_json(const UnitCellCoord &ucc_val, jsonParser &fill_json);

  template<typename T> struct jsonConstructor;

  template<>
  struct jsonConstructor<UnitCellCoord> {

    /// \brief Read from json [b, i, j, k], using 'unit' for UnitCellCoord::unit()
    static UnitCellCoord from_json(const jsonParser &json, const Structure &unit);
  };

  /// \brief Read from json [b, i, j, k]
  void from_json(UnitCellCoord &fill_value, const jsonParser &read_json);


  /* -- UnitCellCoord Definitions ------------------------------------- */

  inline UnitCellCoord::UnitCellCoord(const UnitType &unit) :
    m_unit(&unit) {}

  inline UnitCellCoord::UnitCellCoord(const UnitType &unit, Index _sublat, const UnitCell &_unitcell) :
    m_unit(&unit),
    m_unitcell(_unitcell),
    m_sublat(_sublat) {}

  inline UnitCellCoord::UnitCellCoord(const UnitType &unit, Index _sublat, Index i, Index j, Index k) :
    m_unit(&unit),
    m_unitcell(i, j, k),
    m_sublat(_sublat) {}

  inline const UnitCellCoord::UnitType &UnitCellCoord::unit() const {
    return *m_unit;
  }

  /// \brief Change unit structure, keeping indices constant
  inline void UnitCellCoord::set_unit(const UnitType &_unit) {
    m_unit = &_unit;
  }

  /// \brief Access the Lattice
  inline const Lattice &UnitCellCoord::lattice() const {
    return unit().lattice();
  }

  /// \brief Get corresponding coordinate
  inline Coordinate UnitCellCoord::coordinate() const {
    return site();
  }

  inline UnitCell &UnitCellCoord::unitcell() {
    return m_unitcell;
  }

  inline const UnitCell &UnitCellCoord::unitcell() const {
    return m_unitcell;
  }

  inline Index &UnitCellCoord::unitcell(Index i) {
    return m_unitcell[i];
  }

  inline const Index &UnitCellCoord::unitcell(Index i) const {
    return m_unitcell[i];
  }

  inline Index &UnitCellCoord::sublat() {
    return m_sublat;
  }

  inline const Index &UnitCellCoord::sublat() const {
    return m_sublat;
  }

  inline Index &UnitCellCoord::operator[](Index i) {
    if(i == 0) {
      return m_sublat;
    }
    return unitcell(i - 1);
  }

  inline const Index &UnitCellCoord::operator[](Index i) const {
    if(i == 0) {
      return m_sublat;
    }
    return unitcell(i - 1);
  }

  inline UnitCellCoord &UnitCellCoord::operator+=(UnitCell frac) {
    m_unitcell += frac;
    return *this;
  }

  inline UnitCellCoord &UnitCellCoord::operator-=(UnitCell frac) {
    m_unitcell -= frac;
    return *this;
  }

  /// \brief Compare UnitCellCoord
  inline bool UnitCellCoord::operator<(const UnitCellCoord &B) const {
    const auto &A = *this;
    for(Index i = 0; i < 3; i++) {
      if(A.unitcell()(i) < B.unitcell()(i)) {
        return true;
      }
      if(A.unitcell()(i) > B.unitcell()(i)) {
        return false;
      }
    }
    if(A.sublat() < B.sublat()) {
      return true;
    }

    return false;
  }

  inline bool UnitCellCoord::_eq(const UnitCellCoord &B) const {
    const auto &A = *this;
    return A.unitcell()(0) == B.unitcell()(0) &&
           A.unitcell()(1) == B.unitcell()(1) &&
           A.unitcell()(2) == B.unitcell()(2) &&
           A.sublat() == B.sublat();
  }


  /** @} */
}
#endif





