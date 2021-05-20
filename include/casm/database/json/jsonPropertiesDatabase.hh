#ifndef CASM_jsonPropertiesDatabase
#define CASM_jsonPropertiesDatabase

#include <boost/filesystem/path.hpp>

#include "casm/database/PropertiesDatabase.hh"

namespace CASM {
namespace DB {

class jsonPropertiesDatabase;

class jsonPropertiesDatabaseIterator : public PropertiesDatabaseIteratorBase {
 public:
  jsonPropertiesDatabaseIterator() {}

  std::unique_ptr<jsonPropertiesDatabaseIterator> clone() const {
    return std::unique_ptr<jsonPropertiesDatabaseIterator>(this->_clone());
  }

 private:
  typedef typename std::map<std::string, MappedProperties>::const_iterator
      base_iterator;
  friend jsonPropertiesDatabase;

  jsonPropertiesDatabaseIterator(base_iterator _it) : m_it(_it) {}

  base_iterator base() const { return m_it; }

  bool equal(const PropertiesDatabaseIteratorBase &other) const override {
    return m_it ==
           static_cast<const jsonPropertiesDatabaseIterator &>(other).m_it;
  }

  void increment() override { ++m_it; }

  const MappedProperties &dereference() const override { return m_it->second; }

  long distance_to(const PropertiesDatabaseIteratorBase &other) const override {
    return std::distance(
        m_it, static_cast<const jsonPropertiesDatabaseIterator &>(other).m_it);
  }

  jsonPropertiesDatabaseIterator *_clone() const override {
    return new jsonPropertiesDatabaseIterator(*this);
  }

  base_iterator m_it;
};

/// An implementation of PropertiesDatabase for reading/writing JSON
class jsonPropertiesDatabase : public PropertiesDatabase {
 public:
  /// Constructor
  ///
  /// \param location Where the JSON is read from on "open", written to on
  /// "commit". Can be used all in memory with empty location.
  ///
  /// Note: "_primclex" and "calc_type" are unused and will be removed in the
  /// future
  jsonPropertiesDatabase(const PrimClex &_primclex, std::string calc_type,
                         fs::path location);

  DatabaseBase &open() override;

  void commit() override;

  void close() override;

  /// Clear all and read from JSON
  void from_json(jsonParser const &json);

  /// Export all to JSON
  jsonParser &to_json(jsonParser &json) const;

  /// \brief Begin iterator
  iterator begin() const override;

  /// \brief End iterator
  iterator end() const override;

  size_type size() const override;

  /// \brief Return iterator to MappedProperties that is the best mapping to
  /// specified config
  iterator find_via_to(std::string to_configname) const override;

  /// \brief Return iterator to MappedProperties that is from the specified
  /// config
  iterator find_via_origin(std::string origin) const override;

  /// \brief Names of all configurations that relaxed 'origin'->'to'
  std::set<std::string, Compare> all_origins(
      std::string to_configname) const override;

  /// \brief Change the score method for a single configuration
  void set_score_method(std::string to_configname,
                        const ScoreMappedProperties &score) override;

  /// \brief Change the default score method
  void set_default_score_method(const ScoreMappedProperties &score) override;

  /// \brief Get default score method
  ScoreMappedProperties default_score_method() const override;

 private:
  iterator _iterator(jsonPropertiesDatabaseIterator::base_iterator _it) const;

  /// \brief Private _insert MappedProperties, without modifying 'm_origins'
  std::pair<iterator, bool> _insert(const MappedProperties &value) override;

  /// \brief Private _erase MappedProperties, without modifying 'm_origins'
  iterator _erase(iterator pos) override;

  /// \brief Names of all configurations that relaxed 'origin'->'to'
  void _set_all_origins(std::string to_configname,
                        const std::set<std::string, Compare> &_set) override;

  std::set<std::string, Compare> _make_set(
      std::string to_configname, const ScoreMappedProperties &score) const;

  bool m_is_open;

  std::string m_calc_type;
  fs::path m_location;

  ScoreMappedProperties m_default_score;

  // the MappedProperties container, origin -> MappedProperties
  std::map<std::string, MappedProperties> m_data;

  // to -> {origin, origin, ...}, used to find best mapping
  std::map<std::string, std::set<std::string, Compare> > m_origins;
};

}  // namespace DB
}  // namespace CASM

#endif
