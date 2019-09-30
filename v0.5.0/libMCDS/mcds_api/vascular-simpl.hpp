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

#ifndef VASCULAR_SIMPL_HPP
#define VASCULAR_SIMPL_HPP

#include <xsde/cxx/pre.hxx>

// Begin prologue.
//
//
// End prologue.

#ifndef XSDE_OMIT_SAGGR
#  define XSDE_OMIT_SAGGR
#  define VASCULAR_SIMPL_HPP_CLEAR_OMIT_SAGGR
#endif

#include "vascular-sskel.hpp"

#include <xsde/cxx/stack.hxx>

#include "common-simpl.hpp"

#include "mesh-simpl.hpp"

#include "variables-simpl.hpp"

#include "phenotype_common-simpl.hpp"

namespace vascular
{
  class vascular_node_simpl: public vascular_node_sskel
  {
    public:
    vascular_node_simpl ();

    virtual void
    pre (const ::mesh::node&);

    virtual void
    pre (const ::vascular::vascular_node&);

    // Attributes.
    //
    virtual bool
    boundary_node_present ();

    virtual bool
    boundary_node ();

    public:
    ::mesh::node_simpl base_impl_;

    public:
    struct vascular_node_simpl_state
    {
      const ::vascular::vascular_node* vascular_node_;
    };

    vascular_node_simpl_state vascular_node_simpl_state_;
  };

  class list_of_vascular_nodes_simpl: public list_of_vascular_nodes_sskel
  {
    public:
    list_of_vascular_nodes_simpl ();

    virtual void
    pre (const ::vascular::list_of_vascular_nodes&);

    // Elements.
    //
    virtual bool
    vascular_node_next ();

    virtual const ::vascular::vascular_node&
    vascular_node ();

    virtual bool
    custom_present ();

    virtual const ::common::custom&
    custom ();

    virtual void
    post ();

    virtual void
    _reset ();

    public:
    struct list_of_vascular_nodes_simpl_state
    {
      const ::vascular::list_of_vascular_nodes* list_of_vascular_nodes_;
      ::vascular::list_of_vascular_nodes::vascular_node_const_iterator vascular_node_;
      ::vascular::list_of_vascular_nodes::vascular_node_const_iterator vascular_node_end_;
    };

    list_of_vascular_nodes_simpl_state list_of_vascular_nodes_simpl_state_first_;
    ::xsde::cxx::stack list_of_vascular_nodes_simpl_state_;
  };

  class boundary_node_simpl: public boundary_node_sskel
  {
    public:
    virtual void
    pre (const ::vascular::boundary_node&);

    // Attributes.
    //
    virtual bool
    node_ID_present ();

    virtual unsigned int
    node_ID ();

    // Elements.
    //
    virtual bool
    fluid_flow_velocity_present ();

    virtual const ::common::units_decimal&
    fluid_flow_velocity ();

    virtual bool
    variables_present ();

    virtual const ::variables::list_of_variables&
    variables ();

    virtual bool
    boundary_conditions_present ();

    virtual const ::vascular::boundary_conditions&
    boundary_conditions ();

    virtual bool
    custom_present ();

    virtual const ::common::custom&
    custom ();

    public:
    struct boundary_node_simpl_state
    {
      const ::vascular::boundary_node* boundary_node_;
    };

    boundary_node_simpl_state boundary_node_simpl_state_;
  };

  class list_of_boundary_nodes_simpl: public list_of_boundary_nodes_sskel
  {
    public:
    virtual void
    pre (const ::vascular::list_of_boundary_nodes&);

    // Elements.
    //
    virtual bool
    boundary_node_next ();

    virtual const ::vascular::boundary_node&
    boundary_node ();

    virtual bool
    custom_present ();

    virtual const ::common::custom&
    custom ();

    public:
    struct list_of_boundary_nodes_simpl_state
    {
      const ::vascular::list_of_boundary_nodes* list_of_boundary_nodes_;
      ::vascular::list_of_boundary_nodes::boundary_node_const_iterator boundary_node_;
      ::vascular::list_of_boundary_nodes::boundary_node_const_iterator boundary_node_end_;
    };

    list_of_boundary_nodes_simpl_state list_of_boundary_nodes_simpl_state_;
  };

  class boundary_conditions_simpl: public boundary_conditions_sskel
  {
    public:
    virtual void
    pre (const ::vascular::boundary_conditions&);

    // Attributes.
    //
    virtual bool
    ID_present ();

    virtual unsigned int
    ID ();

    // Elements.
    //
    virtual bool
    boundary_condition_next ();

    virtual const ::vascular::boundary_condition&
    boundary_condition ();

    virtual bool
    custom_present ();

    virtual const ::common::custom&
    custom ();

    public:
    struct boundary_conditions_simpl_state
    {
      const ::vascular::boundary_conditions* boundary_conditions_;
      ::vascular::boundary_conditions::boundary_condition_const_iterator boundary_condition_;
      ::vascular::boundary_conditions::boundary_condition_const_iterator boundary_condition_end_;
    };

    boundary_conditions_simpl_state boundary_conditions_simpl_state_;
  };

  class boundary_type_simpl: public boundary_type_sskel
  {
    public:
    boundary_type_simpl ();

    virtual void
    pre (const ::vascular::boundary_type&);

    virtual void
    _serialize_content ();

    public:
    const ::vascular::boundary_type* boundary_type_simpl_state_;
  };

  class boundary_condition_simpl: public boundary_condition_sskel
  {
    public:
    virtual void
    pre (const ::vascular::boundary_condition&);

    // Attributes.
    //
    virtual bool
    ID_present ();

    virtual unsigned int
    ID ();

    virtual unsigned int
    variable_ID ();

    // Elements.
    //
    virtual const ::vascular::boundary_type&
    boundary_type ();

    virtual bool
    value_present ();

    virtual const ::common::units_decimal&
    value ();

    virtual bool
    direction_present ();

    virtual ::std::string
    direction ();

    virtual bool
    custom_present ();

    virtual const ::common::custom&
    custom ();

    public:
    struct boundary_condition_simpl_state
    {
      const ::vascular::boundary_condition* boundary_condition_;
    };

    boundary_condition_simpl_state boundary_condition_simpl_state_;
  };

  class vascular_segments_simpl: public vascular_segments_sskel
  {
    public:
    virtual void
    pre (const ::vascular::vascular_segments&);

    // Elements.
    //
    virtual bool
    vascular_segment_next ();

    virtual const ::vascular::vascular_segment&
    vascular_segment ();

    virtual bool
    custom_present ();

    virtual const ::common::custom&
    custom ();

    public:
    struct vascular_segments_simpl_state
    {
      const ::vascular::vascular_segments* vascular_segments_;
      ::vascular::vascular_segments::vascular_segment_const_iterator vascular_segment_;
      ::vascular::vascular_segments::vascular_segment_const_iterator vascular_segment_end_;
    };

    vascular_segments_simpl_state vascular_segments_simpl_state_;
  };

  class vascular_segment_simpl: public vascular_segment_sskel
  {
    public:
    virtual void
    pre (const ::vascular::vascular_segment&);

    // Elements.
    //
    virtual const ::vascular::endpoint&
    endpoint_1 ();

    virtual const ::vascular::endpoint&
    endpoint_2 ();

    virtual bool
    surface_present ();

    virtual const ::vascular::surface_properties&
    surface ();

    virtual bool
    interior_present ();

    virtual const ::vascular::volume_properties&
    interior ();

    virtual bool
    custom_present ();

    virtual const ::common::custom&
    custom ();

    public:
    struct vascular_segment_simpl_state
    {
      const ::vascular::vascular_segment* vascular_segment_;
    };

    vascular_segment_simpl_state vascular_segment_simpl_state_;
  };

  class endpoint_simpl: public endpoint_sskel
  {
    public:
    virtual void
    pre (const ::vascular::endpoint&);

    // Attributes.
    //
    virtual bool
    node_ID_present ();

    virtual unsigned int
    node_ID ();

    // Elements.
    //
    virtual bool
    lengths_present ();

    virtual const ::phenotype_common::lengths&
    lengths ();

    virtual bool
    areas_present ();

    virtual const ::phenotype_common::areas_2D&
    areas ();

    virtual bool
    fluid_flow_velocity_present ();

    virtual const ::common::units_decimal&
    fluid_flow_velocity ();

    virtual bool
    custom_present ();

    virtual const ::common::custom&
    custom ();

    public:
    struct endpoint_simpl_state
    {
      const ::vascular::endpoint* endpoint_;
    };

    endpoint_simpl_state endpoint_simpl_state_;
  };

  class surface_properties_simpl: public surface_properties_sskel
  {
    public:
    virtual void
    pre (const ::vascular::surface_properties&);

    // Elements.
    //
    virtual bool
    areas_present ();

    virtual const ::phenotype_common::areas_3D&
    areas ();

    virtual bool
    fluid_flow_velocity_present ();

    virtual const ::common::units_decimal&
    fluid_flow_velocity ();

    virtual bool
    mechanics_present ();

    virtual const ::phenotype_common::mechanics&
    mechanics ();

    virtual bool
    permeability_present ();

    virtual const ::common::units_decimal&
    permeability ();

    virtual bool
    surface_proteins_present ();

    virtual const ::variables::list_of_variables&
    surface_proteins ();

    virtual bool
    custom_present ();

    virtual const ::common::custom&
    custom ();

    public:
    struct surface_properties_simpl_state
    {
      const ::vascular::surface_properties* surface_properties_;
    };

    surface_properties_simpl_state surface_properties_simpl_state_;
  };

  class volume_properties_simpl: public volume_properties_sskel
  {
    public:
    virtual void
    pre (const ::vascular::volume_properties&);

    // Elements.
    //
    virtual bool
    fluid_flow_velocity_present ();

    virtual const ::common::units_decimal&
    fluid_flow_velocity ();

    virtual bool
    variables_present ();

    virtual const ::variables::list_of_variables&
    variables ();

    virtual bool
    volumes_present ();

    virtual const ::phenotype_common::volumes&
    volumes ();

    virtual bool
    custom_present ();

    virtual const ::common::custom&
    custom ();

    public:
    struct volume_properties_simpl_state
    {
      const ::vascular::volume_properties* volume_properties_;
    };

    volume_properties_simpl_state volume_properties_simpl_state_;
  };

  class vascular_network_simpl: public vascular_network_sskel
  {
    public:
    virtual void
    pre (const ::vascular::vascular_network&);

    // Attributes.
    //
    virtual bool
    ID_present ();

    virtual unsigned int
    ID ();

    virtual bool
    keywords_present ();

    virtual ::std::string
    keywords ();

    virtual bool
    name_present ();

    virtual ::std::string
    name ();

    // Elements.
    //
    virtual bool
    vascular_nodes_present ();

    virtual const ::vascular::list_of_vascular_nodes&
    vascular_nodes ();

    virtual bool
    boundary_nodes_present ();

    virtual const ::vascular::list_of_boundary_nodes&
    boundary_nodes ();

    virtual bool
    vascular_segments_present ();

    virtual const ::vascular::vascular_segments&
    vascular_segments ();

    virtual bool
    voxels_present ();

    virtual const ::mesh::int_list_xpath&
    voxels ();

    virtual bool
    custom_present ();

    virtual const ::common::custom&
    custom ();

    public:
    struct vascular_network_simpl_state
    {
      const ::vascular::vascular_network* vascular_network_;
    };

    vascular_network_simpl_state vascular_network_simpl_state_;
  };
}

#ifdef VASCULAR_SIMPL_HPP_CLEAR_OMIT_SAGGR
#  undef XSDE_OMIT_SAGGR
#endif

#ifndef XSDE_OMIT_SAGGR

#endif // XSDE_OMIT_SAGGR

// Begin epilogue.
//
//
// End epilogue.

#include <xsde/cxx/post.hxx>

#endif // VASCULAR_SIMPL_HPP
