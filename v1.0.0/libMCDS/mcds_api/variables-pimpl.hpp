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

#ifndef VARIABLES_PIMPL_HPP
#define VARIABLES_PIMPL_HPP

#include <xsde/cxx/pre.hxx>

// Begin prologue.
//
//
// End prologue.

#ifndef XSDE_OMIT_PAGGR
#  define XSDE_OMIT_PAGGR
#  define VARIABLES_PIMPL_HPP_CLEAR_OMIT_PAGGR
#endif

#include "variables-pskel.hpp"

#include <xsde/cxx/stack.hxx>

#include "common-pimpl.hpp"

namespace variables
{
  class amount_type_pimpl: public amount_type_pskel
  {
    public:
    amount_type_pimpl ();

    virtual void
    pre ();

    virtual void
    _characters (const ::xsde::cxx::ro_string&);

    virtual void
    _post ();

    virtual ::variables::amount_type
    post_amount_type ();

    public:
    struct amount_type_pimpl_state
    {
      ::std::string str_;
    };

    amount_type_pimpl_state amount_type_pimpl_state_;
  };

  class variable_pimpl: public variable_pskel
  {
    public:
    variable_pimpl (bool = false);

    ~variable_pimpl ();

    virtual void
    _reset ();

    virtual void
    pre ();

    // Attributes.
    //
    virtual void
    name (const ::std::string&);

    virtual void
    units (const ::std::string&);

    virtual void
    ID (unsigned long long);

    virtual void
    type (const ::variables::amount_type&);

    virtual void
    ChEBI_ID (const ::std::string&);

    virtual void
    MeSH_ID (const ::std::string&);

    virtual void
    DrugBank_ID (const ::std::string&);

    virtual void
    GMO_ID (const ::std::string&);

    virtual void
    GO_ID (const ::std::string&);

    virtual void
    UniProt_ID (const ::std::string&);

    virtual void
    PR_ID (const ::std::string&);

    // Elements.
    //
    virtual void
    material_amount (::variables::material_amount*);

    virtual void
    physical_parameter_set (::variables::physical_parameter_set*);

    virtual ::variables::variable*
    post_variable ();

    public:
    void
    pre_impl (::variables::variable*);

    public:
    struct variable_pimpl_state
    {
      ::variables::variable* variable_;
    };

    variable_pimpl_state variable_pimpl_state_;
    bool variable_pimpl_base_;
  };

  class material_amount_pimpl: public material_amount_pskel
  {
    public:
    material_amount_pimpl (bool = false);

    ~material_amount_pimpl ();

    virtual void
    _reset ();

    virtual void
    pre ();

    // Attributes.
    //
    virtual void
    type (const ::variables::amount_type&);

    virtual void
    scale_units (const ::std::string&);

    virtual ::variables::material_amount*
    post_material_amount ();

    public:
    void
    pre_impl (::variables::material_amount*);

    public:
    ::common::units_decimal_pimpl base_impl_;

    public:
    struct material_amount_pimpl_state
    {
      ::variables::material_amount* material_amount_;
    };

    material_amount_pimpl_state material_amount_pimpl_state_;
    bool material_amount_pimpl_base_;
  };

  class physical_parameter_set_pimpl: public physical_parameter_set_pskel
  {
    public:
    physical_parameter_set_pimpl (bool = false);

    ~physical_parameter_set_pimpl ();

    virtual void
    _reset ();

    virtual void
    pre ();

    // Elements.
    //
    virtual void
    conditions (::variables::physical_conditions*);

    virtual void
    diffusion_coefficient (::common::units_decimal*);

    virtual void
    decay_rate (::common::units_decimal*);

    virtual void
    custom (::common::custom*);

    virtual ::variables::physical_parameter_set*
    post_physical_parameter_set ();

    public:
    void
    pre_impl (::variables::physical_parameter_set*);

    public:
    struct physical_parameter_set_pimpl_state
    {
      ::variables::physical_parameter_set* physical_parameter_set_;
    };

    physical_parameter_set_pimpl_state physical_parameter_set_pimpl_state_;
    bool physical_parameter_set_pimpl_base_;
  };

  class physical_conditions_pimpl: public physical_conditions_pskel
  {
    public:
    physical_conditions_pimpl (bool = false);

    ~physical_conditions_pimpl ();

    virtual void
    _reset ();

    virtual void
    pre ();

    // Elements.
    //
    virtual void
    temperature (::common::units_decimal*);

    virtual void
    mechanical_pressure (::common::units_decimal*);

    virtual void
    acidity (::common::units_decimal*);

    virtual void
    pH (::common::units_decimal*);

    virtual void
    custom (::common::custom*);

    virtual ::variables::physical_conditions*
    post_physical_conditions ();

    public:
    void
    pre_impl (::variables::physical_conditions*);

    public:
    struct physical_conditions_pimpl_state
    {
      ::variables::physical_conditions* physical_conditions_;
    };

    physical_conditions_pimpl_state physical_conditions_pimpl_state_;
    bool physical_conditions_pimpl_base_;
  };

  class system_pimpl: public system_pskel
  {
    public:
    system_pimpl ();

    virtual void
    pre ();

    virtual void
    _characters (const ::xsde::cxx::ro_string&);

    virtual void
    _post ();

    virtual ::variables::system
    post_system ();

    public:
    struct system_pimpl_state
    {
      ::std::string str_;
    };

    system_pimpl_state system_pimpl_state_;
  };

  class conditions_pimpl: public conditions_pskel
  {
    public:
    conditions_pimpl ();

    virtual void
    pre ();

    virtual void
    _characters (const ::xsde::cxx::ro_string&);

    virtual void
    _post ();

    virtual ::variables::conditions
    post_conditions ();

    public:
    struct conditions_pimpl_state
    {
      ::std::string str_;
    };

    conditions_pimpl_state conditions_pimpl_state_;
  };

  class experimental_conditions_pimpl: public experimental_conditions_pskel
  {
    public:
    experimental_conditions_pimpl (bool = false);

    ~experimental_conditions_pimpl ();

    virtual void
    _reset ();

    virtual void
    pre ();

    // Attributes.
    //
    virtual void
    type (const ::std::string&);

    // Elements.
    //
    virtual void
    dimensionality (unsigned short);

    virtual void
    system (const ::variables::system&);

    virtual void
    conditions (const ::variables::conditions&);

    virtual void
    surface_variable (::variables::variable*);

    virtual ::variables::experimental_conditions*
    post_experimental_conditions ();

    public:
    void
    pre_impl (::variables::experimental_conditions*);

    public:
    struct experimental_conditions_pimpl_state
    {
      ::variables::experimental_conditions* experimental_conditions_;
    };

    experimental_conditions_pimpl_state experimental_conditions_pimpl_state_;
    bool experimental_conditions_pimpl_base_;
  };

  class data_vector_pimpl: public data_vector_pskel
  {
    public:
    data_vector_pimpl (bool = false);

    ~data_vector_pimpl ();

    virtual void
    _reset ();

    virtual void
    pre ();

    // Attributes.
    //
    virtual void
    voxel_ID (::common::unsigned_int_list*);

    virtual ::variables::data_vector*
    post_data_vector ();

    public:
    void
    pre_impl (::variables::data_vector*);

    public:
    ::common::units_double_list_pimpl base_impl_;

    public:
    struct data_vector_pimpl_state
    {
      ::variables::data_vector* data_vector_;
    };

    data_vector_pimpl_state data_vector_pimpl_state_;
    bool data_vector_pimpl_base_;
  };

  class data_pimpl: public data_pskel
  {
    public:
    data_pimpl (bool = false);

    ~data_pimpl ();

    virtual void
    _reset ();

    virtual void
    pre ();

    // Attributes.
    //
    virtual void
    type (const ::common::data_storage_formats&);

    // Elements.
    //
    virtual void
    filename (const ::std::string&);

    virtual void
    data_vector (::variables::data_vector*);

    virtual void
    custom (::common::custom*);

    virtual ::variables::data*
    post_data ();

    public:
    void
    pre_impl (::variables::data*);

    public:
    struct data_pimpl_state
    {
      ::variables::data* data_;
    };

    data_pimpl_state data_pimpl_state_;
    bool data_pimpl_base_;
  };

  class list_of_variables_pimpl: public list_of_variables_pskel
  {
    public:
    list_of_variables_pimpl (bool = false);

    ~list_of_variables_pimpl ();

    virtual void
    _reset ();

    virtual void
    pre ();

    // Elements.
    //
    virtual void
    variable (::variables::variable*);

    virtual void
    physical_parameter_set (::variables::physical_parameter_set*);

    virtual void
    custom (::common::custom*);

    virtual ::variables::list_of_variables*
    post_list_of_variables ();

    public:
    void
    pre_impl (::variables::list_of_variables*);

    public:
    struct list_of_variables_pimpl_state
    {
      ::variables::list_of_variables* list_of_variables_;
    };

    list_of_variables_pimpl_state list_of_variables_pimpl_state_;
    bool list_of_variables_pimpl_base_;
  };

  class transition_threshold_pimpl: public transition_threshold_pskel
  {
    public:
    transition_threshold_pimpl (bool = false);

    ~transition_threshold_pimpl ();

    virtual void
    _reset ();

    virtual void
    pre ();

    // Attributes.
    //
    virtual void
    ChEBI_ID (const ::std::string&);

    virtual void
    MeSH_ID (const ::std::string&);

    virtual void
    DrugBank_ID (const ::std::string&);

    virtual void
    GMO_ID (const ::std::string&);

    virtual void
    GO_ID (const ::std::string&);

    virtual void
    UniProt_ID (const ::std::string&);

    virtual void
    PR_ID (const ::std::string&);

    virtual ::variables::transition_threshold*
    post_transition_threshold1 ();

    public:
    void
    pre_impl (::variables::transition_threshold*);

    public:
    ::common::transition_threshold_pimpl base_impl_;

    public:
    struct transition_threshold_pimpl_state
    {
      ::variables::transition_threshold* transition_threshold_;
    };

    transition_threshold_pimpl_state transition_threshold_pimpl_state_;
    bool transition_threshold_pimpl_base_;
  };
}

#ifdef VARIABLES_PIMPL_HPP_CLEAR_OMIT_PAGGR
#  undef XSDE_OMIT_PAGGR
#endif

#ifndef XSDE_OMIT_PAGGR

#endif // XSDE_OMIT_PAGGR

// Begin epilogue.
//
//
// End epilogue.

#include <xsde/cxx/post.hxx>

#endif // VARIABLES_PIMPL_HPP
