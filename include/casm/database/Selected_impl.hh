#include "casm/database/Selected.hh"
#include "casm/database/DatabaseTypeDefs.hh"

namespace CASM {
  namespace DB {

    template<typename ObjType>
    void Selected<ObjType>::init(const ObjType &_tmplt) const {
      if(!m_selection) {
        if(m_selection_name.empty()) {
          m_selection_name = "MASTER";
        }
        m_selection = notstd::make_cloneable<Selection<ObjType> >(
                        _tmplt.primclex().template db<ObjType>(),
                        m_selection_name);
      }
      else if(m_selection_name.empty()) {
        m_selection_name = m_selection->name();
        if(m_selection_name.empty()) {
          m_selection_name = "unknown";
        }
      }
    }

    template<typename ObjType>
    std::unique_ptr<Selected<ObjType> > Selected<ObjType>::clone() const {
      return std::unique_ptr<Selected>(this->_clone());
    }

    template<typename ObjType>
    Selected<ObjType> *Selected<ObjType>::_clone() const {
      return new Selected(*this);
    }

    template<typename ObjType>
    std::string Selected<ObjType>::short_header(const ObjType &_obj) const {
      return this->name() + "(" + m_selection_name + ")";
    }

    template<typename ObjType>
    bool Selected<ObjType>::evaluate(const ObjType &_obj) const {
      return m_selection->is_selected(_obj.name());
    }

    template<typename ObjType>
    bool Selected<ObjType>::parse_args(const std::string &args) {
      if(m_selection->data().size() || m_selection_name.size()) {
        return false;
      }

      m_selection_name = args;
      return true;
    }
  }
}