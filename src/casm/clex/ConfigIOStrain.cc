#include "casm/clex/ConfigIOStrain.hh"

#include <functional>

#include "casm/clex/ConfigIO.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

namespace ConfigIO {
bool RelaxationStrain::parse_args(const std::string &args) {
  std::vector<std::string> splt_vec;
  boost::split(splt_vec, args, boost::is_any_of(","), boost::token_compress_on);
  std::string tmetric_name = "GL";
  std::string index_expr;

  if (splt_vec.size()) {
    tmetric_name = splt_vec[0];
  }
  if (splt_vec.size() >= 2) {
    index_expr = splt_vec[1];
  }
  if (splt_vec.size() == 3) {
    m_calctype = splt_vec[2];
  }
  if (splt_vec.size() > 3) {
    std::stringstream ss;
    ss << "Too many arguments for 'relaxation_strain'.  Received: " << args
       << "\n";
    throw std::runtime_error(ss.str());
  }
  // if(m_metric_name.size() > 0 && tmetric_name != m_metric_name) {
  //   return false;
  // }
  m_metric_name = tmetric_name;
  if (index_expr.size() > 0) {
    _parse_index_expression(index_expr);
  }

  return false;
}

//****************************************************************************************

bool RelaxationStrain::init(const Configuration &_tmplt) const {
  if (m_metric_name.size() == 0) m_metric_name = "GL";
  m_straincalc.set_mode(m_metric_name);
  if (_index_rules().size() > 0) return true;

  for (Index i = 0; i < 6; i++) _add_rule(std::vector<Index>({i}));
  return true;
}
//****************************************************************************************

bool RelaxationStrain::validate(const Configuration &_config) const {
  return has_strain_property(_config.calc_properties(m_calctype));
}

//****************************************************************************************

std::vector<std::string> RelaxationStrain::col_header(
    const Configuration &_tmplt) const {
  std::vector<std::string> col;
  auto it(_index_rules().cbegin()), end_it(_index_rules().cend());
  // Index s = max(8 - int(name().size()), 0);
  for (; it != end_it; ++it) {
    std::stringstream t_ss;
    t_ss << "    " << name() << '(' << m_metric_name << ',' << (*it)[0] << ','
         << m_calctype << ')';
    col.push_back(t_ss.str());
  }
  return col;
}

//****************************************************************************************

std::string RelaxationStrain::short_header(const Configuration &_tmplt) const {
  return name() + "(" + m_metric_name + ")";
}

//****************************************************************************************
Eigen::VectorXd RelaxationStrain::evaluate(const Configuration &_config) const {
  MappedProperties const &properties = _config.calc_properties(m_calctype);
  DoFKey strain_key = get_strain_property_key(properties);
  Eigen::VectorXd unrolled_metric = properties.global.at(strain_key);
  StrainConverter c(xtal::get_strain_metric(strain_key));
  Eigen::Matrix3d F = c.unrolled_strain_metric_to_F(unrolled_metric);
  return m_straincalc.unrolled_strain_metric(F);
}

//****************************************************************************************

bool DoFStrain::parse_args(const std::string &args) {
  std::vector<std::string> splt_vec;
  boost::split(splt_vec, args, boost::is_any_of(","), boost::token_compress_on);
  std::string tmetric_name = "GL";
  std::string index_expr;

  if (splt_vec.size()) {
    tmetric_name = splt_vec[0];
  }
  if (splt_vec.size() == 2) {
    index_expr = splt_vec[1];
  }
  if (splt_vec.size() > 2) {
    std::stringstream ss;
    ss << "Too many arguments for 'dof_strain'.  Received: " << args << "\n";
    throw std::runtime_error(ss.str());
  }
  // if(m_metric_name.size() > 0 && tmetric_name != m_metric_name) {
  //   return false;
  // }
  m_metric_name = tmetric_name;
  if (index_expr.size() > 0) {
    _parse_index_expression(index_expr);
  }

  return false;
}

bool DoFStrain::init(const Configuration &_tmplt) const {
  if (m_metric_name.size() == 0) m_metric_name = "GL";
  m_straincalc.set_mode(m_metric_name);
  if (_index_rules().size() > 0) return true;

  Index size = evaluate(_tmplt).size();
  for (Index i = 0; i < size; i++) _add_rule(std::vector<Index>({i}));
  return true;
}

std::vector<std::string> DoFStrain::col_header(
    const Configuration &_tmplt) const {
  std::vector<std::string> col;
  auto it(_index_rules().cbegin()), end_it(_index_rules().cend());
  // Index s = max(8 - int(name().size()), 0);
  for (; it != end_it; ++it) {
    std::stringstream t_ss;
    t_ss << name() << '(' << m_metric_name << ',' << (*it)[0] << ')';
    col.push_back(t_ss.str());
  }
  return col;
}

std::string DoFStrain::short_header(const Configuration &_tmplt) const {
  return name() + "(" + m_metric_name + ")";
}

Eigen::VectorXd DoFStrain::evaluate(const Configuration &_config) const {
  if (!m_prim_straincalc ||
      (_config.supercell().shared_prim() != m_shared_prim)) {
    m_shared_prim = _config.supercell().shared_prim();
    if (!has_strain_dof(m_shared_prim->structure())) {
      std::stringstream msg;
      msg << "Error in DoFStrain: Prim does not have strain DoF.";
      throw std::runtime_error(msg.str());
    }
    m_dof_key = get_strain_dof_key(m_shared_prim->structure());
    m_prim_straincalc = notstd::make_unique<StrainConverter>(
        xtal::get_strain_metric(m_dof_key));
  }
  Eigen::VectorXd unrolled_metric =
      _config.configdof().global_dof(m_dof_key).standard_values();
  Eigen::Matrix3d F =
      m_prim_straincalc->unrolled_strain_metric_to_F(unrolled_metric);
  return m_straincalc.unrolled_strain_metric(F);
}

}  // namespace ConfigIO

}  // namespace CASM
