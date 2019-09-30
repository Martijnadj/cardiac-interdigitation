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

#ifndef CELL_SIMPL_HPP
#define CELL_SIMPL_HPP

#include <xsde/cxx/pre.hxx>

// Begin prologue.
//
//
// End prologue.

#ifndef XSDE_OMIT_SAGGR
#  define XSDE_OMIT_SAGGR
#  define CELL_SIMPL_HPP_CLEAR_OMIT_SAGGR
#endif

#include "cell-sskel.hpp"

#include <xsde/cxx/stack.hxx>

#include "common-simpl.hpp"

#include "phenotype_dataset-simpl.hpp"

#include "mesh-simpl.hpp"

#include "cell_line-simpl.hpp"

#include "state-simpl.hpp"

namespace cell
{
  class population_definition_simpl: public population_definition_sskel
  {
    public:
    virtual void
    pre (const ::cell::population_definition&);

    // Attributes.
    //
    virtual bool
    ID_present ();

    virtual unsigned int
    ID ();

    virtual bool
    name_present ();

    virtual ::std::string
    name ();

    virtual bool
    units_present ();

    virtual ::std::string
    units ();

    // Elements.
    //
    virtual bool
    phenotype_dataset_present ();

    virtual const ::phenotype_dataset::phenotype_dataset&
    phenotype_dataset ();

    virtual bool
    custom_present ();

    virtual const ::common::custom&
    custom ();

    public:
    struct population_definition_simpl_state
    {
      const ::cell::population_definition* population_definition_;
    };

    population_definition_simpl_state population_definition_simpl_state_;
  };

  class population_definitions_simpl: public population_definitions_sskel
  {
    public:
    virtual void
    pre (const ::cell::population_definitions&);

    // Elements.
    //
    virtual bool
    population_definition_next ();

    virtual const ::cell::population_definition&
    population_definition ();

    virtual bool
    custom_present ();

    virtual const ::common::custom&
    custom ();

    public:
    struct population_definitions_simpl_state
    {
      const ::cell::population_definitions* population_definitions_;
      ::cell::population_definitions::population_definition_const_iterator population_definition_;
      ::cell::population_definitions::population_definition_const_iterator population_definition_end_;
    };

    population_definitions_simpl_state population_definitions_simpl_state_;
  };

  class cell_simpl: public cell_sskel
  {
    public:
    virtual void
    pre (const ::cell::cell&);

    // Attributes.
    //
    virtual bool
    ID_present ();

    virtual unsigned int
    ID ();

    // Elements.
    //
    virtual bool
    phenotype_dataset_present ();

    virtual const ::phenotype_dataset::phenotype_dataset&
    phenotype_dataset ();

    virtual bool
    state_present ();

    virtual const ::state::state&
    state ();

    virtual bool
    custom_present ();

    virtual const ::common::custom&
    custom ();

    public:
    struct cell_simpl_state
    {
      const ::cell::cell* cell_;
    };

    cell_simpl_state cell_simpl_state_;
  };

  class cell_population_individual_simpl: public cell_population_individual_sskel
  {
    public:
    virtual void
    pre (const ::cell::cell_population_individual&);

    // Attributes.
    //
    virtual bool
    type_present ();

    virtual ::std::string
    type ();

    virtual bool
    population_ID_present ();

    virtual unsigned int
    population_ID ();

    // Elements.
    //
    virtual bool
    cell_next ();

    virtual const ::cell::cell&
    cell ();

    virtual bool
    custom_present ();

    virtual const ::common::custom&
    custom ();

    public:
    struct cell_population_individual_simpl_state
    {
      const ::cell::cell_population_individual* cell_population_individual_;
      ::cell::cell_population_individual::cell_const_iterator cell_;
      ::cell::cell_population_individual::cell_const_iterator cell_end_;
    };

    cell_population_individual_simpl_state cell_population_individual_simpl_state_;
  };

  class cell_population_aggregate_simpl: public cell_population_aggregate_sskel
  {
    public:
    virtual void
    pre (const ::cell::cell_population_aggregate&);

    // Attributes.
    //
    virtual bool
    type_present ();

    virtual ::std::string
    type ();

    virtual bool
    population_ID_present ();

    virtual unsigned int
    population_ID ();

    // Elements.
    //
    virtual bool
    value_present ();

    virtual const ::common::units_decimal&
    value ();

    virtual bool
    sequence_present ();

    virtual bool
    phenotype_dataset_present ();

    virtual const ::phenotype_dataset::phenotype_dataset&
    phenotype_dataset ();

    virtual bool
    state_present ();

    virtual const ::state::state&
    state ();

    virtual bool
    custom_present ();

    virtual const ::common::custom&
    custom ();

    public:
    struct cell_population_aggregate_simpl_state
    {
      const ::cell::cell_population_aggregate* cell_population_aggregate_;
    };

    cell_population_aggregate_simpl_state cell_population_aggregate_simpl_state_;
  };

  class population_vector_simpl: public population_vector_sskel
  {
    public:
    virtual void
    pre (const ::cell::population_vector&);

    // Attributes.
    //
    virtual bool
    voxel_ID_present ();

    virtual const ::common::unsigned_int_list&
    voxel_ID ();

    // Elements.
    //
    virtual bool
    value_present ();

    virtual const ::common::units_double_list&
    value ();

    virtual bool
    cell_population_next ();

    virtual const ::cell::cell_population_aggregate&
    cell_population ();

    virtual bool
    custom_present ();

    virtual const ::common::custom&
    custom ();

    public:
    struct population_vector_simpl_state
    {
      const ::cell::population_vector* population_vector_;
      ::cell::population_vector::cell_population_const_iterator cell_population_;
      ::cell::population_vector::cell_population_const_iterator cell_population_end_;
    };

    population_vector_simpl_state population_vector_simpl_state_;
  };

  class cell_populations_simpl: public cell_populations_sskel
  {
    public:
    virtual void
    pre (const ::cell::cell_populations&);

    // Elements.
    //
    virtual bool
    population_vector_next ();

    virtual const ::cell::population_vector&
    population_vector ();

    virtual bool
    cell_population_present ();

    virtual const ::cell::cell_population_individual&
    cell_population ();

    public:
    struct cell_populations_simpl_state
    {
      const ::cell::cell_populations* cell_populations_;
      ::cell::cell_populations::population_vector_const_iterator population_vector_;
      ::cell::cell_populations::population_vector_const_iterator population_vector_end_;
    };

    cell_populations_simpl_state cell_populations_simpl_state_;
  };

  class cellular_information_simpl: public cellular_information_sskel
  {
    public:
    virtual void
    pre (const ::cell::cellular_information&);

    // Elements.
    //
    virtual bool
    DCLs_present ();

    virtual const ::cell_line::DCLs&
    DCLs ();

    virtual bool
    population_definitions_present ();

    virtual const ::cell::population_definitions&
    population_definitions ();

    virtual bool
    mesh_present ();

    virtual const ::mesh::mesh&
    mesh ();

    virtual bool
    cell_populations_present ();

    virtual const ::cell::cell_populations&
    cell_populations ();

    public:
    struct cellular_information_simpl_state
    {
      const ::cell::cellular_information* cellular_information_;
    };

    cellular_information_simpl_state cellular_information_simpl_state_;
  };
}

#ifdef CELL_SIMPL_HPP_CLEAR_OMIT_SAGGR
#  undef XSDE_OMIT_SAGGR
#endif

#ifndef XSDE_OMIT_SAGGR

#endif // XSDE_OMIT_SAGGR

// Begin epilogue.
//
//
// End epilogue.

#include <xsde/cxx/post.hxx>

#endif // CELL_SIMPL_HPP
