// Copyright (c) 2005-2016 Code Synthesis Tools CC
//
// This program was generated by CodeSynthesis XSD/e, an XML Schema
// to C++ data binding compiler for embedded systems.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License version 2 as
// published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
//
// Furthermore, Code Synthesis Tools CC makes a special exception for
// the Free/Libre and Open Source Software (FLOSS) which is described
// in the accompanying FLOSSE file.
//

// Begin prologue.
//
//
// End prologue.

#include "variables-pimpl.hpp"

#include <xsde/cxx/parser/validating/string-common.hxx>

namespace variables
{
  // amount_type_pimpl
  //

  amount_type_pimpl::
  amount_type_pimpl ()
  : amount_type_pskel (0)
  {
  }

  void amount_type_pimpl::
  pre ()
  {
    this->amount_type_pimpl_state_.str_.clear ();
  }

  void amount_type_pimpl::
  _characters (const ::xsde::cxx::ro_string& s)
  {
    if (this->_facets ().whitespace_ == 2 &&
        this->amount_type_pimpl_state_.str_.size () == 0)
    {
      ::xsde::cxx::ro_string tmp (s.data (), s.size ());

      if (::xsde::cxx::trim_left (tmp) != 0)
      {
        this->amount_type_pimpl_state_.str_ += tmp;
      }
    }
    else
      this->amount_type_pimpl_state_.str_ += s;
  }

  void amount_type_pimpl::
  _post ()
  {
    ::xsde::cxx::parser::validating::string_common::validate_facets (
      this->amount_type_pimpl_state_.str_,
      this->_facets (),
      this->_context ());
  }

  ::variables::amount_type amount_type_pimpl::
  post_amount_type ()
  {
    ::variables::amount_type::value_type v =
    static_cast< ::variables::amount_type::value_type > (0);
    const char* s = this->amount_type_pimpl_state_.str_.c_str ();

    if (strcmp (s, "concentration") == 0)
      v = ::variables::amount_type::concentration;
    else if (strcmp (s, "density") == 0)
      v = ::variables::amount_type::density;
    else if (strcmp (s, "volume_fraction") == 0)
      v = ::variables::amount_type::volume_fraction;
    else if (strcmp (s, "volume_percent") == 0)
      v = ::variables::amount_type::volume_percent;
    else if (strcmp (s, "volume_percentage") == 0)
      v = ::variables::amount_type::volume_percentage;
    else if (strcmp (s, "surface_density") == 0)
      v = ::variables::amount_type::surface_density;
    else if (strcmp (s, "area_fraction") == 0)
      v = ::variables::amount_type::area_fraction;
    else if (strcmp (s, "area_percent") == 0)
      v = ::variables::amount_type::area_percent;
    else if (strcmp (s, "area_percentage") == 0)
      v = ::variables::amount_type::area_percentage;
    else if (strcmp (s, "count") == 0)
      v = ::variables::amount_type::count;
    else if (strcmp (s, "partial_pressure") == 0)
      v = ::variables::amount_type::partial_pressure;
    else if (strcmp (s, "surface") == 0)
      v = ::variables::amount_type::surface;

    ::variables::amount_type r (v);
    return r;
  }

  // variable_pimpl
  //

  variable_pimpl::
  variable_pimpl (bool b)
  {
    this->variable_pimpl_base_ = b;
    this->variable_pimpl_state_.variable_ = 0;
  }

  variable_pimpl::
  ~variable_pimpl ()
  {
    if (!this->variable_pimpl_base_ && this->variable_pimpl_state_.variable_)
      delete this->variable_pimpl_state_.variable_;
  }

  void variable_pimpl::
  _reset ()
  {
    variable_pskel::_reset ();

    if (!this->variable_pimpl_base_ && this->variable_pimpl_state_.variable_)
    {
      delete this->variable_pimpl_state_.variable_;
      this->variable_pimpl_state_.variable_ = 0;
    }
  }

  void variable_pimpl::
  pre_impl (::variables::variable* x)
  {
    this->variable_pimpl_state_.variable_ = x;
  }

  void variable_pimpl::
  pre ()
  {
    ::variables::variable* x = new ::variables::variable;
    this->pre_impl (x);
  }

  void variable_pimpl::
  name (const ::std::string& x)
  {
    this->variable_pimpl_state_.variable_->name (x);
  }

  void variable_pimpl::
  units (const ::std::string& x)
  {
    this->variable_pimpl_state_.variable_->units (x);
  }

  void variable_pimpl::
  ID (unsigned long long x)
  {
    this->variable_pimpl_state_.variable_->ID (x);
  }

  void variable_pimpl::
  type (const ::variables::amount_type& x)
  {
    this->variable_pimpl_state_.variable_->type (x);
  }

  void variable_pimpl::
  ChEBI_ID (const ::std::string& x)
  {
    this->variable_pimpl_state_.variable_->ChEBI_ID (x);
  }

  void variable_pimpl::
  MeSH_ID (const ::std::string& x)
  {
    this->variable_pimpl_state_.variable_->MeSH_ID (x);
  }

  void variable_pimpl::
  DrugBank_ID (const ::std::string& x)
  {
    this->variable_pimpl_state_.variable_->DrugBank_ID (x);
  }

  void variable_pimpl::
  GMO_ID (const ::std::string& x)
  {
    this->variable_pimpl_state_.variable_->GMO_ID (x);
  }

  void variable_pimpl::
  GO_ID (const ::std::string& x)
  {
    this->variable_pimpl_state_.variable_->GO_ID (x);
  }

  void variable_pimpl::
  UniProt_ID (const ::std::string& x)
  {
    this->variable_pimpl_state_.variable_->UniProt_ID (x);
  }

  void variable_pimpl::
  PR_ID (const ::std::string& x)
  {
    this->variable_pimpl_state_.variable_->PR_ID (x);
  }

  void variable_pimpl::
  material_amount (::variables::material_amount* x)
  {
    this->variable_pimpl_state_.variable_->material_amount (x);
  }

  void variable_pimpl::
  physical_parameter_set (::variables::physical_parameter_set* x)
  {
    this->variable_pimpl_state_.variable_->physical_parameter_set (x);
  }

  ::variables::variable* variable_pimpl::
  post_variable ()
  {
    ::variables::variable* r = this->variable_pimpl_state_.variable_;
    this->variable_pimpl_state_.variable_ = 0;
    return r;
  }

  // material_amount_pimpl
  //

  material_amount_pimpl::
  material_amount_pimpl (bool b)
  : material_amount_pskel (&base_impl_),
    base_impl_ (true)
  {
    this->material_amount_pimpl_base_ = b;
    this->material_amount_pimpl_state_.material_amount_ = 0;
  }

  material_amount_pimpl::
  ~material_amount_pimpl ()
  {
    if (!this->material_amount_pimpl_base_ && this->material_amount_pimpl_state_.material_amount_)
      delete this->material_amount_pimpl_state_.material_amount_;
  }

  void material_amount_pimpl::
  _reset ()
  {
    material_amount_pskel::_reset ();

    if (!this->material_amount_pimpl_base_ && this->material_amount_pimpl_state_.material_amount_)
    {
      delete this->material_amount_pimpl_state_.material_amount_;
      this->material_amount_pimpl_state_.material_amount_ = 0;
    }
  }

  void material_amount_pimpl::
  pre_impl (::variables::material_amount* x)
  {
    this->material_amount_pimpl_state_.material_amount_ = x;
    this->base_impl_.pre_impl (x);
  }

  void material_amount_pimpl::
  pre ()
  {
    ::variables::material_amount* x = new ::variables::material_amount;
    this->pre_impl (x);
  }

  void material_amount_pimpl::
  type (const ::variables::amount_type& x)
  {
    this->material_amount_pimpl_state_.material_amount_->type (x);
  }

  void material_amount_pimpl::
  scale_units (const ::std::string& x)
  {
    this->material_amount_pimpl_state_.material_amount_->scale_units (x);
  }

  ::variables::material_amount* material_amount_pimpl::
  post_material_amount ()
  {
    this->base_impl_.post_units_decimal ();
    ::variables::material_amount* r = this->material_amount_pimpl_state_.material_amount_;
    this->material_amount_pimpl_state_.material_amount_ = 0;
    return r;
  }

  // physical_parameter_set_pimpl
  //

  physical_parameter_set_pimpl::
  physical_parameter_set_pimpl (bool b)
  {
    this->physical_parameter_set_pimpl_base_ = b;
    this->physical_parameter_set_pimpl_state_.physical_parameter_set_ = 0;
  }

  physical_parameter_set_pimpl::
  ~physical_parameter_set_pimpl ()
  {
    if (!this->physical_parameter_set_pimpl_base_ && this->physical_parameter_set_pimpl_state_.physical_parameter_set_)
      delete this->physical_parameter_set_pimpl_state_.physical_parameter_set_;
  }

  void physical_parameter_set_pimpl::
  _reset ()
  {
    physical_parameter_set_pskel::_reset ();

    if (!this->physical_parameter_set_pimpl_base_ && this->physical_parameter_set_pimpl_state_.physical_parameter_set_)
    {
      delete this->physical_parameter_set_pimpl_state_.physical_parameter_set_;
      this->physical_parameter_set_pimpl_state_.physical_parameter_set_ = 0;
    }
  }

  void physical_parameter_set_pimpl::
  pre_impl (::variables::physical_parameter_set* x)
  {
    this->physical_parameter_set_pimpl_state_.physical_parameter_set_ = x;
  }

  void physical_parameter_set_pimpl::
  pre ()
  {
    ::variables::physical_parameter_set* x = new ::variables::physical_parameter_set;
    this->pre_impl (x);
  }

  void physical_parameter_set_pimpl::
  conditions (::variables::physical_conditions* x)
  {
    this->physical_parameter_set_pimpl_state_.physical_parameter_set_->conditions (x);
  }

  void physical_parameter_set_pimpl::
  diffusion_coefficient (::common::units_decimal* x)
  {
    this->physical_parameter_set_pimpl_state_.physical_parameter_set_->diffusion_coefficient (x);
  }

  void physical_parameter_set_pimpl::
  decay_rate (::common::units_decimal* x)
  {
    this->physical_parameter_set_pimpl_state_.physical_parameter_set_->decay_rate (x);
  }

  void physical_parameter_set_pimpl::
  custom (::common::custom* x)
  {
    this->physical_parameter_set_pimpl_state_.physical_parameter_set_->custom (x);
  }

  ::variables::physical_parameter_set* physical_parameter_set_pimpl::
  post_physical_parameter_set ()
  {
    ::variables::physical_parameter_set* r = this->physical_parameter_set_pimpl_state_.physical_parameter_set_;
    this->physical_parameter_set_pimpl_state_.physical_parameter_set_ = 0;
    return r;
  }

  // physical_conditions_pimpl
  //

  physical_conditions_pimpl::
  physical_conditions_pimpl (bool b)
  {
    this->physical_conditions_pimpl_base_ = b;
    this->physical_conditions_pimpl_state_.physical_conditions_ = 0;
  }

  physical_conditions_pimpl::
  ~physical_conditions_pimpl ()
  {
    if (!this->physical_conditions_pimpl_base_ && this->physical_conditions_pimpl_state_.physical_conditions_)
      delete this->physical_conditions_pimpl_state_.physical_conditions_;
  }

  void physical_conditions_pimpl::
  _reset ()
  {
    physical_conditions_pskel::_reset ();

    if (!this->physical_conditions_pimpl_base_ && this->physical_conditions_pimpl_state_.physical_conditions_)
    {
      delete this->physical_conditions_pimpl_state_.physical_conditions_;
      this->physical_conditions_pimpl_state_.physical_conditions_ = 0;
    }
  }

  void physical_conditions_pimpl::
  pre_impl (::variables::physical_conditions* x)
  {
    this->physical_conditions_pimpl_state_.physical_conditions_ = x;
  }

  void physical_conditions_pimpl::
  pre ()
  {
    ::variables::physical_conditions* x = new ::variables::physical_conditions;
    this->pre_impl (x);
  }

  void physical_conditions_pimpl::
  temperature (::common::units_decimal* x)
  {
    this->physical_conditions_pimpl_state_.physical_conditions_->temperature (x);
  }

  void physical_conditions_pimpl::
  mechanical_pressure (::common::units_decimal* x)
  {
    this->physical_conditions_pimpl_state_.physical_conditions_->mechanical_pressure (x);
  }

  void physical_conditions_pimpl::
  acidity (::common::units_decimal* x)
  {
    this->physical_conditions_pimpl_state_.physical_conditions_->acidity (x);
  }

  void physical_conditions_pimpl::
  pH (::common::units_decimal* x)
  {
    this->physical_conditions_pimpl_state_.physical_conditions_->pH (x);
  }

  void physical_conditions_pimpl::
  custom (::common::custom* x)
  {
    this->physical_conditions_pimpl_state_.physical_conditions_->custom (x);
  }

  ::variables::physical_conditions* physical_conditions_pimpl::
  post_physical_conditions ()
  {
    ::variables::physical_conditions* r = this->physical_conditions_pimpl_state_.physical_conditions_;
    this->physical_conditions_pimpl_state_.physical_conditions_ = 0;
    return r;
  }

  // system_pimpl
  //

  system_pimpl::
  system_pimpl ()
  : system_pskel (0)
  {
  }

  void system_pimpl::
  pre ()
  {
    this->system_pimpl_state_.str_.clear ();
  }

  void system_pimpl::
  _characters (const ::xsde::cxx::ro_string& s)
  {
    if (this->_facets ().whitespace_ == 2 &&
        this->system_pimpl_state_.str_.size () == 0)
    {
      ::xsde::cxx::ro_string tmp (s.data (), s.size ());

      if (::xsde::cxx::trim_left (tmp) != 0)
      {
        this->system_pimpl_state_.str_ += tmp;
      }
    }
    else
      this->system_pimpl_state_.str_ += s;
  }

  void system_pimpl::
  _post ()
  {
    ::xsde::cxx::parser::validating::string_common::validate_facets (
      this->system_pimpl_state_.str_,
      this->_facets (),
      this->_context ());
  }

  ::variables::system system_pimpl::
  post_system ()
  {
    ::variables::system::value_type v =
    static_cast< ::variables::system::value_type > (0);
    const char* s = this->system_pimpl_state_.str_.c_str ();

    if (strcmp (s, "in vivo") == 0)
      v = ::variables::system::in_vivo;
    else if (strcmp (s, "in vitro") == 0)
      v = ::variables::system::in_vitro;
    else if (strcmp (s, "ex vivo") == 0)
      v = ::variables::system::ex_vivo;
    else if (strcmp (s, "in silico") == 0)
      v = ::variables::system::in_silico;

    ::variables::system r (v);
    return r;
  }

  // conditions_pimpl
  //

  conditions_pimpl::
  conditions_pimpl ()
  : conditions_pskel (0)
  {
  }

  void conditions_pimpl::
  pre ()
  {
    this->conditions_pimpl_state_.str_.clear ();
  }

  void conditions_pimpl::
  _characters (const ::xsde::cxx::ro_string& s)
  {
    if (this->_facets ().whitespace_ == 2 &&
        this->conditions_pimpl_state_.str_.size () == 0)
    {
      ::xsde::cxx::ro_string tmp (s.data (), s.size ());

      if (::xsde::cxx::trim_left (tmp) != 0)
      {
        this->conditions_pimpl_state_.str_ += tmp;
      }
    }
    else
      this->conditions_pimpl_state_.str_ += s;
  }

  void conditions_pimpl::
  _post ()
  {
    ::xsde::cxx::parser::validating::string_common::validate_facets (
      this->conditions_pimpl_state_.str_,
      this->_facets (),
      this->_context ());
  }

  ::variables::conditions conditions_pimpl::
  post_conditions ()
  {
    ::variables::conditions::value_type v =
    static_cast< ::variables::conditions::value_type > (0);
    const char* s = this->conditions_pimpl_state_.str_.c_str ();

    if (strcmp (s, "surface") == 0)
      v = ::variables::conditions::surface;
    else if (strcmp (s, "suspension") == 0)
      v = ::variables::conditions::suspension;
    else if (strcmp (s, "spheroid") == 0)
      v = ::variables::conditions::spheroid;

    ::variables::conditions r (v);
    return r;
  }

  // experimental_conditions_pimpl
  //

  experimental_conditions_pimpl::
  experimental_conditions_pimpl (bool b)
  {
    this->experimental_conditions_pimpl_base_ = b;
    this->experimental_conditions_pimpl_state_.experimental_conditions_ = 0;
  }

  experimental_conditions_pimpl::
  ~experimental_conditions_pimpl ()
  {
    if (!this->experimental_conditions_pimpl_base_ && this->experimental_conditions_pimpl_state_.experimental_conditions_)
      delete this->experimental_conditions_pimpl_state_.experimental_conditions_;
  }

  void experimental_conditions_pimpl::
  _reset ()
  {
    experimental_conditions_pskel::_reset ();

    if (!this->experimental_conditions_pimpl_base_ && this->experimental_conditions_pimpl_state_.experimental_conditions_)
    {
      delete this->experimental_conditions_pimpl_state_.experimental_conditions_;
      this->experimental_conditions_pimpl_state_.experimental_conditions_ = 0;
    }
  }

  void experimental_conditions_pimpl::
  pre_impl (::variables::experimental_conditions* x)
  {
    this->experimental_conditions_pimpl_state_.experimental_conditions_ = x;
  }

  void experimental_conditions_pimpl::
  pre ()
  {
    ::variables::experimental_conditions* x = new ::variables::experimental_conditions;
    this->pre_impl (x);
  }

  void experimental_conditions_pimpl::
  type (const ::std::string& x)
  {
    this->experimental_conditions_pimpl_state_.experimental_conditions_->type (x);
  }

  void experimental_conditions_pimpl::
  dimensionality (unsigned short x)
  {
    this->experimental_conditions_pimpl_state_.experimental_conditions_->dimensionality (x);
  }

  void experimental_conditions_pimpl::
  system (const ::variables::system& x)
  {
    this->experimental_conditions_pimpl_state_.experimental_conditions_->system (x);
  }

  void experimental_conditions_pimpl::
  conditions (const ::variables::conditions& x)
  {
    this->experimental_conditions_pimpl_state_.experimental_conditions_->conditions (x);
  }

  void experimental_conditions_pimpl::
  surface_variable (::variables::variable* x)
  {
    this->experimental_conditions_pimpl_state_.experimental_conditions_->surface_variable ().push_back (x);
  }

  ::variables::experimental_conditions* experimental_conditions_pimpl::
  post_experimental_conditions ()
  {
    ::variables::experimental_conditions* r = this->experimental_conditions_pimpl_state_.experimental_conditions_;
    this->experimental_conditions_pimpl_state_.experimental_conditions_ = 0;
    return r;
  }

  // data_vector_pimpl
  //

  data_vector_pimpl::
  data_vector_pimpl (bool b)
  : data_vector_pskel (&base_impl_),
    base_impl_ (true)
  {
    this->data_vector_pimpl_base_ = b;
    this->data_vector_pimpl_state_.data_vector_ = 0;
  }

  data_vector_pimpl::
  ~data_vector_pimpl ()
  {
    if (!this->data_vector_pimpl_base_ && this->data_vector_pimpl_state_.data_vector_)
      delete this->data_vector_pimpl_state_.data_vector_;
  }

  void data_vector_pimpl::
  _reset ()
  {
    data_vector_pskel::_reset ();

    if (!this->data_vector_pimpl_base_ && this->data_vector_pimpl_state_.data_vector_)
    {
      delete this->data_vector_pimpl_state_.data_vector_;
      this->data_vector_pimpl_state_.data_vector_ = 0;
    }
  }

  void data_vector_pimpl::
  pre_impl (::variables::data_vector* x)
  {
    this->data_vector_pimpl_state_.data_vector_ = x;
    this->base_impl_.pre_impl (x);
  }

  void data_vector_pimpl::
  pre ()
  {
    ::variables::data_vector* x = new ::variables::data_vector;
    this->pre_impl (x);
  }

  void data_vector_pimpl::
  voxel_ID (::common::unsigned_int_list* x)
  {
    this->data_vector_pimpl_state_.data_vector_->voxel_ID (x);
  }

  ::variables::data_vector* data_vector_pimpl::
  post_data_vector ()
  {
    this->base_impl_.post_units_double_list ();
    ::variables::data_vector* r = this->data_vector_pimpl_state_.data_vector_;
    this->data_vector_pimpl_state_.data_vector_ = 0;
    return r;
  }

  // data_pimpl
  //

  data_pimpl::
  data_pimpl (bool b)
  {
    this->data_pimpl_base_ = b;
    this->data_pimpl_state_.data_ = 0;
  }

  data_pimpl::
  ~data_pimpl ()
  {
    if (!this->data_pimpl_base_ && this->data_pimpl_state_.data_)
      delete this->data_pimpl_state_.data_;
  }

  void data_pimpl::
  _reset ()
  {
    data_pskel::_reset ();

    if (!this->data_pimpl_base_ && this->data_pimpl_state_.data_)
    {
      delete this->data_pimpl_state_.data_;
      this->data_pimpl_state_.data_ = 0;
    }
  }

  void data_pimpl::
  pre_impl (::variables::data* x)
  {
    this->data_pimpl_state_.data_ = x;
  }

  void data_pimpl::
  pre ()
  {
    ::variables::data* x = new ::variables::data;
    this->pre_impl (x);
  }

  void data_pimpl::
  type (const ::common::data_storage_formats& x)
  {
    this->data_pimpl_state_.data_->type (x);
  }

  void data_pimpl::
  filename (const ::std::string& x)
  {
    this->data_pimpl_state_.data_->filename (x);
  }

  void data_pimpl::
  data_vector (::variables::data_vector* x)
  {
    this->data_pimpl_state_.data_->data_vector ().push_back (x);
  }

  void data_pimpl::
  custom (::common::custom* x)
  {
    this->data_pimpl_state_.data_->custom (x);
  }

  ::variables::data* data_pimpl::
  post_data ()
  {
    ::variables::data* r = this->data_pimpl_state_.data_;
    this->data_pimpl_state_.data_ = 0;
    return r;
  }

  // list_of_variables_pimpl
  //

  list_of_variables_pimpl::
  list_of_variables_pimpl (bool b)
  {
    this->list_of_variables_pimpl_base_ = b;
    this->list_of_variables_pimpl_state_.list_of_variables_ = 0;
  }

  list_of_variables_pimpl::
  ~list_of_variables_pimpl ()
  {
    if (!this->list_of_variables_pimpl_base_ && this->list_of_variables_pimpl_state_.list_of_variables_)
      delete this->list_of_variables_pimpl_state_.list_of_variables_;
  }

  void list_of_variables_pimpl::
  _reset ()
  {
    list_of_variables_pskel::_reset ();

    if (!this->list_of_variables_pimpl_base_ && this->list_of_variables_pimpl_state_.list_of_variables_)
    {
      delete this->list_of_variables_pimpl_state_.list_of_variables_;
      this->list_of_variables_pimpl_state_.list_of_variables_ = 0;
    }
  }

  void list_of_variables_pimpl::
  pre_impl (::variables::list_of_variables* x)
  {
    this->list_of_variables_pimpl_state_.list_of_variables_ = x;
  }

  void list_of_variables_pimpl::
  pre ()
  {
    ::variables::list_of_variables* x = new ::variables::list_of_variables;
    this->pre_impl (x);
  }

  void list_of_variables_pimpl::
  variable (::variables::variable* x)
  {
    this->list_of_variables_pimpl_state_.list_of_variables_->variable ().push_back (x);
  }

  void list_of_variables_pimpl::
  physical_parameter_set (::variables::physical_parameter_set* x)
  {
    this->list_of_variables_pimpl_state_.list_of_variables_->physical_parameter_set (x);
  }

  void list_of_variables_pimpl::
  custom (::common::custom* x)
  {
    this->list_of_variables_pimpl_state_.list_of_variables_->custom (x);
  }

  ::variables::list_of_variables* list_of_variables_pimpl::
  post_list_of_variables ()
  {
    ::variables::list_of_variables* r = this->list_of_variables_pimpl_state_.list_of_variables_;
    this->list_of_variables_pimpl_state_.list_of_variables_ = 0;
    return r;
  }

  // transition_threshold_pimpl
  //

  transition_threshold_pimpl::
  transition_threshold_pimpl (bool b)
  : transition_threshold_pskel (&base_impl_),
    base_impl_ (true)
  {
    this->transition_threshold_pimpl_base_ = b;
    this->transition_threshold_pimpl_state_.transition_threshold_ = 0;
  }

  transition_threshold_pimpl::
  ~transition_threshold_pimpl ()
  {
    if (!this->transition_threshold_pimpl_base_ && this->transition_threshold_pimpl_state_.transition_threshold_)
      delete this->transition_threshold_pimpl_state_.transition_threshold_;
  }

  void transition_threshold_pimpl::
  _reset ()
  {
    transition_threshold_pskel::_reset ();

    if (!this->transition_threshold_pimpl_base_ && this->transition_threshold_pimpl_state_.transition_threshold_)
    {
      delete this->transition_threshold_pimpl_state_.transition_threshold_;
      this->transition_threshold_pimpl_state_.transition_threshold_ = 0;
    }
  }

  void transition_threshold_pimpl::
  pre_impl (::variables::transition_threshold* x)
  {
    this->transition_threshold_pimpl_state_.transition_threshold_ = x;
    this->base_impl_.pre_impl (x);
  }

  void transition_threshold_pimpl::
  pre ()
  {
    ::variables::transition_threshold* x = new ::variables::transition_threshold;
    this->pre_impl (x);
  }

  void transition_threshold_pimpl::
  ChEBI_ID (const ::std::string& x)
  {
    this->transition_threshold_pimpl_state_.transition_threshold_->ChEBI_ID (x);
  }

  void transition_threshold_pimpl::
  MeSH_ID (const ::std::string& x)
  {
    this->transition_threshold_pimpl_state_.transition_threshold_->MeSH_ID (x);
  }

  void transition_threshold_pimpl::
  DrugBank_ID (const ::std::string& x)
  {
    this->transition_threshold_pimpl_state_.transition_threshold_->DrugBank_ID (x);
  }

  void transition_threshold_pimpl::
  GMO_ID (const ::std::string& x)
  {
    this->transition_threshold_pimpl_state_.transition_threshold_->GMO_ID (x);
  }

  void transition_threshold_pimpl::
  GO_ID (const ::std::string& x)
  {
    this->transition_threshold_pimpl_state_.transition_threshold_->GO_ID (x);
  }

  void transition_threshold_pimpl::
  UniProt_ID (const ::std::string& x)
  {
    this->transition_threshold_pimpl_state_.transition_threshold_->UniProt_ID (x);
  }

  void transition_threshold_pimpl::
  PR_ID (const ::std::string& x)
  {
    this->transition_threshold_pimpl_state_.transition_threshold_->PR_ID (x);
  }

  ::variables::transition_threshold* transition_threshold_pimpl::
  post_transition_threshold1 ()
  {
    this->base_impl_.post_transition_threshold ();
    ::variables::transition_threshold* r = this->transition_threshold_pimpl_state_.transition_threshold_;
    this->transition_threshold_pimpl_state_.transition_threshold_ = 0;
    return r;
  }
}

// Begin epilogue.
//
//
// End epilogue.

