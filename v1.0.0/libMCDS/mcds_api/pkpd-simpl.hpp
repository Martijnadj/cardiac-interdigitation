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

#ifndef PKPD_SIMPL_HPP
#define PKPD_SIMPL_HPP

#include <xsde/cxx/pre.hxx>

// Begin prologue.
//
//
// End prologue.

#ifndef XSDE_OMIT_SAGGR
#  define XSDE_OMIT_SAGGR
#  define PKPD_SIMPL_HPP_CLEAR_OMIT_SAGGR
#endif

#include "pkpd-sskel.hpp"

#include <xsde/cxx/stack.hxx>

#include "common-simpl.hpp"

#include "variables-simpl.hpp"

namespace pkpd
{
  class pharmacokinetics_simpl: public pharmacokinetics_sskel
  {
    public:
    virtual void
    pre (const ::pkpd::pharmacokinetics&);

    // Elements.
    //
    virtual bool
    inactivation_rate_present ();

    virtual const ::common::units_decimal&
    inactivation_rate ();

    virtual bool
    half_life_present ();

    virtual const ::common::units_decimal&
    half_life ();

    public:
    struct pharmacokinetics_simpl_state
    {
      const ::pkpd::pharmacokinetics* pharmacokinetics_;
    };

    pharmacokinetics_simpl_state pharmacokinetics_simpl_state_;
  };

  class drug_simpl: public drug_sskel
  {
    public:
    virtual void
    pre (const ::pkpd::drug&);

    // Attributes.
    //
    // Elements.
    //
    virtual choice_arm_tag
    choice_arm ();

    virtual const ::pkpd::dose&
    dose ();

    virtual const ::pkpd::pharmacokinetics&
    pharmacokinetics ();

    public:
    struct drug_simpl_state
    {
      const ::pkpd::drug* drug_;
    };

    drug_simpl_state drug_simpl_state_;
  };

  class drug_dose_simpl: public drug_dose_sskel
  {
    public:
    virtual void
    pre (const ::pkpd::drug_dose&);

    // Attributes.
    //
    virtual bool
    ID_present ();

    virtual unsigned short
    ID ();

    virtual bool
    ChEBI_ID_present ();

    virtual ::std::string
    ChEBI_ID ();

    virtual bool
    MeSH_ID_present ();

    virtual ::std::string
    MeSH_ID ();

    virtual bool
    DrugBank_ID_present ();

    virtual ::std::string
    DrugBank_ID ();

    virtual bool
    GMO_ID_present ();

    virtual ::std::string
    GMO_ID ();

    virtual bool
    GO_ID_present ();

    virtual ::std::string
    GO_ID ();

    virtual bool
    UniProt_ID_present ();

    virtual ::std::string
    UniProt_ID ();

    virtual bool
    PR_ID_present ();

    virtual ::std::string
    PR_ID ();

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
    virtual const ::pkpd::dose&
    dose ();

    public:
    struct drug_dose_simpl_state
    {
      const ::pkpd::drug_dose* drug_dose_;
    };

    drug_dose_simpl_state drug_dose_simpl_state_;
  };

  class drug_pk_simpl: public drug_pk_sskel
  {
    public:
    virtual void
    pre (const ::pkpd::drug_pk&);

    // Attributes.
    //
    virtual bool
    ID_present ();

    virtual unsigned short
    ID ();

    virtual bool
    ChEBI_ID_present ();

    virtual ::std::string
    ChEBI_ID ();

    virtual bool
    MeSH_ID_present ();

    virtual ::std::string
    MeSH_ID ();

    virtual bool
    DrugBank_ID_present ();

    virtual ::std::string
    DrugBank_ID ();

    virtual bool
    GMO_ID_present ();

    virtual ::std::string
    GMO_ID ();

    virtual bool
    GO_ID_present ();

    virtual ::std::string
    GO_ID ();

    virtual bool
    UniProt_ID_present ();

    virtual ::std::string
    UniProt_ID ();

    virtual bool
    PR_ID_present ();

    virtual ::std::string
    PR_ID ();

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
    pharmacokinetics_present ();

    virtual const ::pkpd::pharmacokinetics&
    pharmacokinetics ();

    public:
    struct drug_pk_simpl_state
    {
      const ::pkpd::drug_pk* drug_pk_;
    };

    drug_pk_simpl_state drug_pk_simpl_state_;
  };

  class dose_simpl: public dose_sskel
  {
    public:
    dose_simpl ();

    virtual void
    pre (const ::pkpd::dose&);

    // Attributes.
    //
    virtual bool
    type_present ();

    virtual ::std::string
    type ();

    public:
    ::common::units_decimal_simpl base_impl_;

    public:
    struct dose_simpl_state
    {
      const ::pkpd::dose* dose_;
    };

    dose_simpl_state dose_simpl_state_;
  };

  class therapy_simpl: public therapy_sskel
  {
    public:
    virtual void
    pre (const ::pkpd::therapy&);

    // Elements.
    //
    virtual bool
    drug_next ();

    virtual const ::pkpd::drug_dose&
    drug ();

    public:
    struct therapy_simpl_state
    {
      const ::pkpd::therapy* therapy_;
      ::pkpd::therapy::drug_const_iterator drug_;
      ::pkpd::therapy::drug_const_iterator drug_end_;
    };

    therapy_simpl_state therapy_simpl_state_;
  };

  class response_simpl: public response_sskel
  {
    public:
    virtual void
    pre (const ::pkpd::response&);

    // Elements.
    //
    virtual bool
    maximum_birth_inhibition_present ();

    virtual const ::common::units_decimal&
    maximum_birth_inhibition ();

    virtual bool
    maximum_birth_inhibition_time_present ();

    virtual const ::common::units_decimal&
    maximum_birth_inhibition_time ();

    virtual bool
    birth_inhibition_recovery_rate_present ();

    virtual const ::common::units_decimal&
    birth_inhibition_recovery_rate ();

    virtual bool
    maximum_death_rate_present ();

    virtual const ::common::units_decimal&
    maximum_death_rate ();

    virtual bool
    maximum_death_time_present ();

    virtual const ::common::units_decimal&
    maximum_death_time ();

    virtual bool
    death_recovery_rate_present ();

    virtual const ::common::units_decimal&
    death_recovery_rate ();

    virtual bool
    response_observation_next ();

    virtual const ::pkpd::response_observation&
    response_observation ();

    virtual bool
    custom_present ();

    virtual const ::common::custom&
    custom ();

    public:
    struct response_simpl_state
    {
      const ::pkpd::response* response_;
      ::pkpd::response::response_observation_const_iterator response_observation_;
      ::pkpd::response::response_observation_const_iterator response_observation_end_;
    };

    response_simpl_state response_simpl_state_;
  };

  class response_observation_simpl: public response_observation_sskel
  {
    public:
    virtual void
    pre (const ::pkpd::response_observation&);

    // Elements.
    //
    virtual bool
    time_present ();

    virtual const ::common::units_decimal&
    time ();

    virtual bool
    birth_rate_present ();

    virtual const ::common::units_decimal_nonnegative&
    birth_rate ();

    virtual bool
    death_rate_present ();

    virtual const ::common::units_decimal_nonnegative&
    death_rate ();

    virtual bool
    net_birth_rate_present ();

    virtual const ::common::units_decimal&
    net_birth_rate ();

    virtual bool
    net_death_rate_present ();

    virtual const ::common::units_decimal&
    net_death_rate ();

    virtual bool
    apoptotic_duration_present ();

    virtual const ::common::units_decimal&
    apoptotic_duration ();

    virtual bool
    percent_cell_viability_present ();

    virtual const ::common::units_decimal&
    percent_cell_viability ();

    virtual bool
    custom_present ();

    virtual const ::common::custom&
    custom ();

    public:
    struct response_observation_simpl_state
    {
      const ::pkpd::response_observation* response_observation_;
    };

    response_observation_simpl_state response_observation_simpl_state_;
  };

  class pharmacodynamics_simpl: public pharmacodynamics_sskel
  {
    public:
    virtual void
    pre (const ::pkpd::pharmacodynamics&);

    // Elements.
    //
    virtual bool
    therapy_measurement_set_next ();

    virtual const ::pkpd::therapy_measurement_set&
    therapy_measurement_set ();

    public:
    struct pharmacodynamics_simpl_state
    {
      const ::pkpd::pharmacodynamics* pharmacodynamics_;
      ::pkpd::pharmacodynamics::therapy_measurement_set_const_iterator therapy_measurement_set_;
      ::pkpd::pharmacodynamics::therapy_measurement_set_const_iterator therapy_measurement_set_end_;
    };

    pharmacodynamics_simpl_state pharmacodynamics_simpl_state_;
  };

  class therapy_measurement_set_simpl: public therapy_measurement_set_sskel
  {
    public:
    virtual void
    pre (const ::pkpd::therapy_measurement_set&);

    // Attributes.
    //
    virtual bool
    ID_present ();

    virtual unsigned short
    ID ();

    // Elements.
    //
    virtual const ::pkpd::therapy&
    therapy ();

    virtual const ::pkpd::response&
    response ();

    public:
    struct therapy_measurement_set_simpl_state
    {
      const ::pkpd::therapy_measurement_set* therapy_measurement_set_;
    };

    therapy_measurement_set_simpl_state therapy_measurement_set_simpl_state_;
  };

  class PKPD_simpl: public PKPD_sskel
  {
    public:
    virtual void
    pre (const ::pkpd::PKPD&);

    // Elements.
    //
    virtual bool
    drug_next ();

    virtual const ::pkpd::drug_pk&
    drug ();

    virtual bool
    pharmacodynamics_present ();

    virtual const ::pkpd::pharmacodynamics&
    pharmacodynamics ();

    public:
    struct PKPD_simpl_state
    {
      const ::pkpd::PKPD* PKPD_;
      ::pkpd::PKPD::drug_const_iterator drug_;
      ::pkpd::PKPD::drug_const_iterator drug_end_;
    };

    PKPD_simpl_state PKPD_simpl_state_;
  };
}

#ifdef PKPD_SIMPL_HPP_CLEAR_OMIT_SAGGR
#  undef XSDE_OMIT_SAGGR
#endif

#ifndef XSDE_OMIT_SAGGR

#endif // XSDE_OMIT_SAGGR

// Begin epilogue.
//
//
// End epilogue.

#include <xsde/cxx/post.hxx>

#endif // PKPD_SIMPL_HPP
