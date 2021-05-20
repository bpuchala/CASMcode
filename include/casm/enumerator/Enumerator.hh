#ifndef CASM_Enumerator
#define CASM_Enumerator

#include <string>

#include "casm/casm_io/json/jsonParser.hh"  // TODO remove this dependency
#include "casm/global/definitions.hh"
#include "casm/misc/CASM_TMP.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/misc/unique_cloneable_map.hh"

namespace CASM {

/** \defgroup Enumerator
 *
 *  \ingroup Container
 *  \brief Algorithms for enumerating objects such as Lattice, Supercell,
 Configuration
 *
 *  Enumerators are classes that implement methods to enumerate Supercells,
 *  Configurations, and other objects by providing iterators over a range of
 *  enumerated objects. The objects are not usually stored in the enumerators
 *  but constructed as iterators are modified or dereferenced. New enumerators
 *  can be created by inheriting from:
 *  - InputEnumeratorBase (for single-pass  enumerators)
 *  - RandomAccessEnumeratorBase (for multi-pass, random-access enumerators)
 *
 *  Some example enumerators are:
 *  - ScelEnumT, ScelEnumByNameT, ScelEnumByPropsT
 *  - ConfigEnumAllOccupations
 *  - ConfigEnumInterpolation
 *  - ConfigEnumStrain
 *  - SuperConfigEnum
 *
 *  Enumerators are required to "know" their own name by implementing a
 *  traits class with 'name' as a const std::string member. For enumerators
 *  meant to be accessible via the 'casm enum' API, the traits class must also
 *  have a 'help' const std::string member explaining the enumerator and the
 *  input parameters, and a 'CASM::EnumInterface<EnumMethod>::run' method which
 *  implements collecting input parameters from the CLI, executing the
 enumerator
 *  and saving the resulting objects. Some helper functions exist to make
 *  parsing input and saving results easier:
 *
 *  - To collect CLI input and construct an ScelEnum which enumerates canonical
 *  Supercell: ::make_enumerator_scel_enum
 *  - To collect CLI input and construct an SuperlatticeEnumerator which
 *  enumerates super-lattices which may not be canonical:
 ::make_enumerator_superlat_enum
 *  - To save results from enumerators of unique, primitive, canonical
 *  Configurations: ::insert_unique_canon_configs
 *  - To save results from enumerators of general Configurations:
 ::insert_configs
 *
 *
 *  For Enumerators only meant to be used internally the required members are:
 *
 *  Variables:
 *  - \code public: static const std::string enumerator_name; \endcode
 *
 *  Functions:
 *  - \code public: std::string name() const override { return enumerator_name;
 } \endcode
 *
 *
 *  For Enumerators to be added to the API the required members are:
 *
 *  Variables:
 *  - \code public: static const std::string enumerator_name; \endcode
 *  - \code public: static std::string interface_help(); \endcode
 *
 *  Functions:
 *  - public: std::string name() const override { return enumerator_name; }
 \endcode
 *  - public: static int run(const PrimClex &primclex, const jsonParser &kwargs,
 const Completer::EnumOption &enum_opt); \endcode
 *
 *  To enable use as a plugin:
 *  - \code
 *    extern "C" {
 *      CASM::EnumInterfaceBase *make_EnumMethod_interface() {
 *        return new CASM::EnumInterface<CASM::ConfigEnumAllOccupations>();
 *      }
 *    }
 *    \endcode
 *  - To use an enumerator as a plugin for an existing CASM project, place the
 *    source code in the `.casm/enumerators` directory in a file named
 *    `EnumMethodName.cc`, where `EnumMethod` is the name of the enumerator
 *    class.

    @{
*/

// ---- Enumerator ---------------------

/// \brief Abstract base class for enumerators
///
/// - To implement a new enumeration method do not inherit from this
///   directly. Instead, inherit from either InputEnum (for single-pass
///   enumerators) or RandomAccessEnum (for multi-pass, random-access
///   enumerators)
/// - The Enumerator base class only holds a pointer to the 'current'
///   object, the current 'step' index, and a 'valid' flag
///
class EnumeratorBase {
 public:
  typedef long step_type;

  /// Default constructor
  EnumeratorBase() : m_valid(false), m_step(0) {}

  virtual ~EnumeratorBase() {}

  /// Increments with each enumerated object
  step_type step() const { return m_step; }

  /// Returns false if enumeration is complete
  bool valid() const { return m_valid; }

  /// Default Object source just uses step#
  ///
  /// Returns:
  /// \code
  /// {
  ///   "enumerated_by": "<enumerator_type>",
  ///   "step": <step #>
  /// }
  /// \endcode
  virtual jsonParser source(step_type step) const {
    jsonParser src;
    src["enumerated_by"] = this->name();
    src["step"] = step;
    return src;
  }

  /// \brief Derived enumerators must implement name, via ENUM_MEMBERS
  virtual std::string name() const = 0;

 protected:
  /// Initialize
  ///
  /// - Sets step to 0
  /// - Sets valid to true
  void _initialize() {
    m_step = 0;
    m_valid = true;
  }

  /// Set current step value
  void _set_step(step_type val) { m_step = val; }

  /// Increment current step value
  void _increment_step() { ++m_step; }

  /// Decrement current step value
  void _decrement_step() { --m_step; }

  /// Call if enumeration complete
  void _invalidate() { m_valid = false; }

  /// Used if random access enumerator step is moved into valid range
  void _validate() { m_valid = true; }

 private:
  bool m_valid;

  step_type m_step;
};

template <typename ValueType, bool IsConst = true>
class ValEnumerator : public EnumeratorBase {
 public:
  typedef ValueType value_type;
  typedef CASM_TMP::ConstSwitch<IsConst, ValueType> &reference;
  using EnumeratorBase::step_type;

  ValEnumerator() : m_current_ptr(nullptr) {}

  virtual ~ValEnumerator() {}

  // -- from EnumeratorBase --

 public:
  using EnumeratorBase::name;    // pure virtual
  using EnumeratorBase::source;  // virtual
  using EnumeratorBase::step;
  using EnumeratorBase::valid;

 protected:
  using EnumeratorBase::_initialize;

  /// Initialize
  ///
  /// - Sets current to point at _initial
  /// - Sets step to 0
  /// - Sets valid to true
  void _initialize(CASM_TMP::ConstSwitch<IsConst, value_type> *_initial) {
    _set_current_ptr(_initial);
    EnumeratorBase::_initialize();
  }

  using EnumeratorBase::_decrement_step;
  using EnumeratorBase::_increment_step;
  using EnumeratorBase::_invalidate;
  using EnumeratorBase::_set_step;

  // -- For ValEnumerator --

 public:
  /// Access the current ObjectType by reference
  reference current() const { return *m_current_ptr; }

 protected:
  /// Change the pointer
  void _set_current_ptr(CASM_TMP::ConstSwitch<IsConst, value_type> *_new) {
    m_current_ptr = _new;
  }

 private:
  CASM_TMP::ConstSwitch<IsConst, value_type> *m_current_ptr;
};

class EnumIteratorBase {
 public:
  typedef EnumeratorBase::step_type step_type;

  /// Default Constructor
  EnumIteratorBase() : m_enum_ptr(nullptr) {}

  /// Constructor
  EnumIteratorBase(EnumeratorBase &enumerator) : m_enum_ptr(&enumerator) {}

  virtual ~EnumIteratorBase() {}

  /// Return current step number
  ///
  /// - Only valid if iterator refers to valid object (not end)
  virtual step_type step() const = 0;

  /// Uses 'step' and enumerator class 'source' implementation
  ///
  /// - Only valid if iterator refers to valid object (not end)
  jsonParser source() const { return _enum().source(this->step()); }

  /// Uses enumerator class 'name' implementation
  std::string name() const { return m_enum_ptr->name(); }

  /// Returns true if 'end' iterator
  virtual bool is_end() const = 0;

  std::unique_ptr<EnumIteratorBase> clone() const {
    return std::unique_ptr<EnumIteratorBase>(this->_clone());
  }

 protected:
  void _assert_same_ptr(const EnumIteratorBase &other) const {
    assert((m_enum_ptr == other.m_enum_ptr) &&
           "Comparing EnumIterator referring to different enumerators is not "
           "allowed!");
  }

  void _assert_ptr() const {
    assert(m_enum_ptr && "EnumIterator does not point to any enumerator!");
  }

  void _assert_valid() const {
    assert(m_enum_ptr->valid() &&
           "EnumIterator points to an invalid enumerator!");
  }

  /// \brief boost::iterator_facade implementation
  ///
  /// - Uses 'is_end' implementation to check if both iterators are 'end'
  /// - If both are not end, then compares iterator 'step'
  bool equal(const EnumIteratorBase &other) const {
    bool this_is_end = this->is_end();
    bool other_is_end = other.is_end();

    if (this_is_end != other_is_end) {
      return false;
    }

    if (this_is_end) {
      return true;
    }

    return (this->step() == other.step());
  }

  EnumeratorBase *_enum_ptr() { return m_enum_ptr; }

  EnumeratorBase *_enum_ptr() const { return m_enum_ptr; }

 private:
  virtual EnumIteratorBase *_clone() const = 0;

  EnumeratorBase &_enum() { return *m_enum_ptr; }

  EnumeratorBase &_enum() const { return *m_enum_ptr; }

  // pointer to Enumerator
  EnumeratorBase *m_enum_ptr;
};

template <typename ValueType, bool IsConst = true>
class ValEnumIterator : public EnumIteratorBase {
 public:
  using EnumIteratorBase::step_type;
  typedef ValueType value_type;
  typedef typename ValEnumerator<ValueType, IsConst>::reference reference;

  ValEnumIterator() {}

  ValEnumIterator(ValEnumerator<ValueType, IsConst> &enumerator)
      : EnumIteratorBase(enumerator) {}

  virtual ~ValEnumIterator() {}

  using EnumIteratorBase::is_end;  // pure virtual
  using EnumIteratorBase::name;
  using EnumIteratorBase::source;
  using EnumIteratorBase::step;  // pure virtual

  std::unique_ptr<EnumIteratorBase> clone() const {
    return std::unique_ptr<EnumIteratorBase>(this->_clone());
  }

 protected:
  using EnumIteratorBase::_assert_ptr;
  using EnumIteratorBase::_assert_same_ptr;
  using EnumIteratorBase::_assert_valid;
  using EnumIteratorBase::_enum_ptr;
  using EnumIteratorBase::equal;

 private:
  virtual EnumIteratorBase *_clone() const = 0;
  virtual reference dereference() const = 0;
};

/// Can be specialized to return true if enumerator output is guaranteed to be
/// canonical (and primitive for Configurations) as required for database
/// insertion
template <typename EnumeratorType>
bool is_guaranteed_for_database_insert(EnumeratorType const &enumerator) {
  return false;
}

/** @}*/

}  // namespace CASM

#endif
