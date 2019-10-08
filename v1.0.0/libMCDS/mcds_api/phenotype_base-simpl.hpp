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

#ifndef PHENOTYPE_BASE_SIMPL_HPP
#define PHENOTYPE_BASE_SIMPL_HPP

#include <xsde/cxx/pre.hxx>

// Begin prologue.
//
//
// End prologue.

#ifndef XSDE_OMIT_SAGGR
#  define XSDE_OMIT_SAGGR
#  define PHENOTYPE_BASE_SIMPL_HPP_CLEAR_OMIT_SAGGR
#endif

#include "phenotype_base-sskel.hpp"

#include <xsde/cxx/stack.hxx>

#include "phenotype_common-simpl.hpp"

#include "common-simpl.hpp"

#include "pkpd-simpl.hpp"

#include "variables-simpl.hpp"

namespace phenotype_base
{
  class phenotype_type_simpl: public phenotype_type_sskel
  {
    public:
    phenotype_type_simpl ();

    virtual void
    pre (const ::phenotype_base::phenotype_type&);

    virtual void
    _serialize_content ();

    public:
    const ::phenotype_base::phenotype_type* phenotype_type_simpl_state_;
  };

  class phenotype_base_simpl: public phenotype_base_sskel
  {
    public:
    virtual void
    pre (const ::phenotype_base::phenotype_base&);

    // Attributes.
    //
    virtual bool
    type_present ();

    virtual const ::phenotype_base::phenotype_type&
    type ();

    // Elements.
    //
    virtual bool
    adhesion_present ();

    virtual const ::phenotype_common::adhesion&
    adhesion ();

    virtual bool
    geometrical_properties_present ();

    virtual const ::phenotype_common::geometrical_properties&
    geometrical_properties ();

    virtual bool
    mass_present ();

    virtual const ::phenotype_common::mass&
    mass ();

    virtual bool
    mechanics_present ();

    virtual const ::phenotype_common::mechanics&
    mechanics ();

    virtual bool
    motility_present ();

    virtual const ::phenotype_common::motility&
    motility ();

    virtual bool
    PKPD_present ();

    virtual const ::pkpd::PKPD&
    PKPD ();

    virtual bool
    timescale_present ();

    virtual const ::phenotype_base::expected_timescale&
    timescale ();

    virtual bool
    transport_processes_present ();

    virtual const ::phenotype_common::transport_processes&
    transport_processes ();

    virtual bool
    custom_present ();

    virtual const ::common::custom&
    custom ();

    public:
    struct phenotype_base_simpl_state
    {
      const ::phenotype_base::phenotype_base* phenotype_base_;
    };

    phenotype_base_simpl_state phenotype_base_simpl_state_;
  };

  class expected_timescale_simpl: public expected_timescale_sskel
  {
    public:
    expected_timescale_simpl ();

    virtual void
    pre (const ::phenotype_base::expected_timescale&);

    // Attributes.
    //
    virtual bool
    cell_cycle_ID_present ();

    virtual unsigned int
    cell_cycle_ID ();

    virtual bool
    cell_cycle_phase_ID_present ();

    virtual unsigned int
    cell_cycle_phase_ID ();

    virtual bool
    cell_death_ID_present ();

    virtual unsigned int
    cell_death_ID ();

    public:
    ::common::units_decimal_nonnegative_simpl base_impl_;

    public:
    struct expected_timescale_simpl_state
    {
      const ::phenotype_base::expected_timescale* expected_timescale_;
    };

    expected_timescale_simpl_state expected_timescale_simpl_state_;
  };

  class cell_parts_simpl: public cell_parts_sskel
  {
    public:
    cell_parts_simpl ();

    virtual void
    pre (const ::phenotype_base::cell_parts&);

    // Attributes.
    //
    virtual ::std::string
    name ();

    virtual bool
    ID_present ();

    virtual unsigned int
    ID ();

    // Elements.
    //
    virtual bool
    phenotype_next ();

    virtual const ::phenotype_base::phenotype_base&
    phenotype ();

    virtual bool
    cell_part_next ();

    virtual const ::phenotype_base::cell_parts&
    cell_part ();

    virtual bool
    custom_present ();

    virtual const ::common::custom&
    custom ();

    virtual void
    post ();

    virtual void
    _reset ();

    public:
    struct cell_parts_simpl_state
    {
      const ::phenotype_base::cell_parts* cell_parts_;
      ::phenotype_base::cell_parts::phenotype_const_iterator phenotype_;
      ::phenotype_base::cell_parts::phenotype_const_iterator phenotype_end_;
      ::phenotype_base::cell_parts::cell_part_const_iterator cell_part_;
      ::phenotype_base::cell_parts::cell_part_const_iterator cell_part_end_;
    };

    cell_parts_simpl_state cell_parts_simpl_state_first_;
    ::xsde::cxx::stack cell_parts_simpl_state_;
  };
}

#ifdef PHENOTYPE_BASE_SIMPL_HPP_CLEAR_OMIT_SAGGR
#  undef XSDE_OMIT_SAGGR
#endif

#ifndef XSDE_OMIT_SAGGR

#endif // XSDE_OMIT_SAGGR

// Begin epilogue.
//
//
// End epilogue.

#include <xsde/cxx/post.hxx>

#endif // PHENOTYPE_BASE_SIMPL_HPP
