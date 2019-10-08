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

#ifndef PKPD_HPP
#define PKPD_HPP

#include <xsde/cxx/version.hxx>

#if (XSDE_INT_VERSION != 3020000L)
#error XSD/e runtime version mismatch
#endif

#include <xsde/cxx/config.hxx>

#ifndef XSDE_ENCODING_UTF8
#error the generated code uses the UTF-8 encodingwhile the XSD/e runtime does not (reconfigure the runtime or change the --char-encoding value)
#endif

#ifndef XSDE_STL
#error the generated code uses STL while the XSD/e runtime does not (reconfigure the runtime or add --no-stl)
#endif

#ifndef XSDE_EXCEPTIONS
#error the generated code uses exceptions while the XSD/e runtime does not (reconfigure the runtime or add --no-exceptions)
#endif

#ifndef XSDE_LONGLONG
#error the generated code uses long long while the XSD/e runtime does not (reconfigure the runtime or add --no-long-long)
#endif

#ifdef XSDE_CUSTOM_ALLOCATOR
#error the XSD/e runtime uses custom allocator while the generated code does not (reconfigure the runtime or add --custom-allocator)
#endif

#include <xsde/cxx/pre.hxx>

// Begin prologue.
//
//
// End prologue.

#include "pkpd-fwd.hpp"

#ifndef XSDE_DONT_INCLUDE_INLINE
#define XSDE_DONT_INCLUDE_INLINE

#include "common.hpp"

#include "variables.hpp"

#undef XSDE_DONT_INCLUDE_INLINE
#else

#include "common.hpp"

#include "variables.hpp"

#endif // XSDE_DONT_INCLUDE_INLINE

namespace pkpd
{
  // pharmacokinetics (variable-length)
  //
  class pharmacokinetics
  {
    private:
    pharmacokinetics (const pharmacokinetics&);
    pharmacokinetics& operator= (const pharmacokinetics&);

    public:
    pharmacokinetics ();

    pharmacokinetics*
    _clone () const;

    ~pharmacokinetics ();

    // inactivation_rate
    //
    bool
    inactivation_rate_present () const;

    const ::common::units_decimal&
    inactivation_rate () const;

    ::common::units_decimal&
    inactivation_rate ();

    void
    inactivation_rate (::common::units_decimal*);

    ::common::units_decimal*
    inactivation_rate_detach ();

    // half_life
    //
    bool
    half_life_present () const;

    const ::common::units_decimal&
    half_life () const;

    ::common::units_decimal&
    half_life ();

    void
    half_life (::common::units_decimal*);

    ::common::units_decimal*
    half_life_detach ();

    void
    _copy (pharmacokinetics&) const;

    private:
    ::common::units_decimal* inactivation_rate_;
    ::common::units_decimal* half_life_;
  };

  // drug (variable-length)
  //
  class drug
  {
    private:
    drug (const drug&);
    drug& operator= (const drug&);

    public:
    drug ();

    drug*
    _clone () const;

    ~drug ();

    // choice
    //
    enum choice_arm_tag
    {
      dose_tag,
      pharmacokinetics_tag
    };

    choice_arm_tag
    choice_arm () const;

    void
    choice_arm (choice_arm_tag);

    // dose
    //
    const ::pkpd::dose&
    dose () const;

    ::pkpd::dose&
    dose ();

    void
    dose (::pkpd::dose*);

    ::pkpd::dose*
    dose_detach ();

    // pharmacokinetics
    //
    const ::pkpd::pharmacokinetics&
    pharmacokinetics () const;

    ::pkpd::pharmacokinetics&
    pharmacokinetics ();

    void
    pharmacokinetics (::pkpd::pharmacokinetics*);

    ::pkpd::pharmacokinetics*
    pharmacokinetics_detach ();

    void
    _copy (drug&) const;

    private:
    union
    {
      ::pkpd::dose* dose_;
      ::pkpd::pharmacokinetics* pharmacokinetics_;
    } choice_;
    choice_arm_tag choice_arm_;
  };

  // drug_dose (variable-length)
  //
  class drug_dose
  {
    private:
    drug_dose (const drug_dose&);
    drug_dose& operator= (const drug_dose&);

    public:
    drug_dose ();

    drug_dose*
    _clone () const;

    ~drug_dose ();

    // ID
    //
    bool
    ID_present () const;

    void
    ID_present (bool);

    unsigned short
    ID () const;

    unsigned short&
    ID ();

    void
    ID (unsigned short);

    // ChEBI_ID
    //
    bool
    ChEBI_ID_present () const;

    void
    ChEBI_ID_present (bool);

    const ::std::string&
    ChEBI_ID () const;

    ::std::string&
    ChEBI_ID ();

    void
    ChEBI_ID (const ::std::string&);

    // MeSH_ID
    //
    bool
    MeSH_ID_present () const;

    void
    MeSH_ID_present (bool);

    const ::std::string&
    MeSH_ID () const;

    ::std::string&
    MeSH_ID ();

    void
    MeSH_ID (const ::std::string&);

    // DrugBank_ID
    //
    bool
    DrugBank_ID_present () const;

    void
    DrugBank_ID_present (bool);

    const ::std::string&
    DrugBank_ID () const;

    ::std::string&
    DrugBank_ID ();

    void
    DrugBank_ID (const ::std::string&);

    // GMO_ID
    //
    bool
    GMO_ID_present () const;

    void
    GMO_ID_present (bool);

    const ::std::string&
    GMO_ID () const;

    ::std::string&
    GMO_ID ();

    void
    GMO_ID (const ::std::string&);

    // GO_ID
    //
    bool
    GO_ID_present () const;

    void
    GO_ID_present (bool);

    const ::std::string&
    GO_ID () const;

    ::std::string&
    GO_ID ();

    void
    GO_ID (const ::std::string&);

    // UniProt_ID
    //
    bool
    UniProt_ID_present () const;

    void
    UniProt_ID_present (bool);

    const ::std::string&
    UniProt_ID () const;

    ::std::string&
    UniProt_ID ();

    void
    UniProt_ID (const ::std::string&);

    // PR_ID
    //
    bool
    PR_ID_present () const;

    void
    PR_ID_present (bool);

    const ::std::string&
    PR_ID () const;

    ::std::string&
    PR_ID ();

    void
    PR_ID (const ::std::string&);

    // name
    //
    bool
    name_present () const;

    void
    name_present (bool);

    const ::std::string&
    name () const;

    ::std::string&
    name ();

    void
    name (const ::std::string&);

    // units
    //
    bool
    units_present () const;

    void
    units_present (bool);

    const ::std::string&
    units () const;

    ::std::string&
    units ();

    void
    units (const ::std::string&);

    // dose
    //
    const ::pkpd::dose&
    dose () const;

    ::pkpd::dose&
    dose ();

    void
    dose (::pkpd::dose*);

    ::pkpd::dose*
    dose_detach ();

    void
    _copy (drug_dose&) const;

    private:
    unsigned short ID_;
    unsigned char ID_present_;
    ::std::string ChEBI_ID_;
    unsigned char ChEBI_ID_present_;
    ::std::string MeSH_ID_;
    unsigned char MeSH_ID_present_;
    ::std::string DrugBank_ID_;
    unsigned char DrugBank_ID_present_;
    ::std::string GMO_ID_;
    unsigned char GMO_ID_present_;
    ::std::string GO_ID_;
    unsigned char GO_ID_present_;
    ::std::string UniProt_ID_;
    unsigned char UniProt_ID_present_;
    ::std::string PR_ID_;
    unsigned char PR_ID_present_;
    ::std::string name_;
    unsigned char name_present_;
    ::std::string units_;
    unsigned char units_present_;
    ::pkpd::dose* dose_;
  };

  // drug_pk (variable-length)
  //
  class drug_pk
  {
    private:
    drug_pk (const drug_pk&);
    drug_pk& operator= (const drug_pk&);

    public:
    drug_pk ();

    drug_pk*
    _clone () const;

    ~drug_pk ();

    // ID
    //
    bool
    ID_present () const;

    void
    ID_present (bool);

    unsigned short
    ID () const;

    unsigned short&
    ID ();

    void
    ID (unsigned short);

    // ChEBI_ID
    //
    bool
    ChEBI_ID_present () const;

    void
    ChEBI_ID_present (bool);

    const ::std::string&
    ChEBI_ID () const;

    ::std::string&
    ChEBI_ID ();

    void
    ChEBI_ID (const ::std::string&);

    // MeSH_ID
    //
    bool
    MeSH_ID_present () const;

    void
    MeSH_ID_present (bool);

    const ::std::string&
    MeSH_ID () const;

    ::std::string&
    MeSH_ID ();

    void
    MeSH_ID (const ::std::string&);

    // DrugBank_ID
    //
    bool
    DrugBank_ID_present () const;

    void
    DrugBank_ID_present (bool);

    const ::std::string&
    DrugBank_ID () const;

    ::std::string&
    DrugBank_ID ();

    void
    DrugBank_ID (const ::std::string&);

    // GMO_ID
    //
    bool
    GMO_ID_present () const;

    void
    GMO_ID_present (bool);

    const ::std::string&
    GMO_ID () const;

    ::std::string&
    GMO_ID ();

    void
    GMO_ID (const ::std::string&);

    // GO_ID
    //
    bool
    GO_ID_present () const;

    void
    GO_ID_present (bool);

    const ::std::string&
    GO_ID () const;

    ::std::string&
    GO_ID ();

    void
    GO_ID (const ::std::string&);

    // UniProt_ID
    //
    bool
    UniProt_ID_present () const;

    void
    UniProt_ID_present (bool);

    const ::std::string&
    UniProt_ID () const;

    ::std::string&
    UniProt_ID ();

    void
    UniProt_ID (const ::std::string&);

    // PR_ID
    //
    bool
    PR_ID_present () const;

    void
    PR_ID_present (bool);

    const ::std::string&
    PR_ID () const;

    ::std::string&
    PR_ID ();

    void
    PR_ID (const ::std::string&);

    // name
    //
    bool
    name_present () const;

    void
    name_present (bool);

    const ::std::string&
    name () const;

    ::std::string&
    name ();

    void
    name (const ::std::string&);

    // units
    //
    bool
    units_present () const;

    void
    units_present (bool);

    const ::std::string&
    units () const;

    ::std::string&
    units ();

    void
    units (const ::std::string&);

    // pharmacokinetics
    //
    bool
    pharmacokinetics_present () const;

    const ::pkpd::pharmacokinetics&
    pharmacokinetics () const;

    ::pkpd::pharmacokinetics&
    pharmacokinetics ();

    void
    pharmacokinetics (::pkpd::pharmacokinetics*);

    ::pkpd::pharmacokinetics*
    pharmacokinetics_detach ();

    void
    _copy (drug_pk&) const;

    private:
    unsigned short ID_;
    unsigned char ID_present_;
    ::std::string ChEBI_ID_;
    unsigned char ChEBI_ID_present_;
    ::std::string MeSH_ID_;
    unsigned char MeSH_ID_present_;
    ::std::string DrugBank_ID_;
    unsigned char DrugBank_ID_present_;
    ::std::string GMO_ID_;
    unsigned char GMO_ID_present_;
    ::std::string GO_ID_;
    unsigned char GO_ID_present_;
    ::std::string UniProt_ID_;
    unsigned char UniProt_ID_present_;
    ::std::string PR_ID_;
    unsigned char PR_ID_present_;
    ::std::string name_;
    unsigned char name_present_;
    ::std::string units_;
    unsigned char units_present_;
    ::pkpd::pharmacokinetics* pharmacokinetics_;
  };

  // dose (variable-length)
  //
  class dose: public ::common::units_decimal
  {
    private:
    dose (const dose&);
    dose& operator= (const dose&);

    public:
    dose ();

    dose*
    _clone () const;

    ~dose ();

    // type
    //
    bool
    type_present () const;

    void
    type_present (bool);

    const ::std::string&
    type () const;

    ::std::string&
    type ();

    void
    type (const ::std::string&);

    void
    _copy (dose&) const;

    private:
    ::std::string type_;
    unsigned char type_present_;
  };

  // therapy (variable-length)
  //
  class therapy
  {
    private:
    therapy (const therapy&);
    therapy& operator= (const therapy&);

    public:
    therapy ();

    therapy*
    _clone () const;

    ~therapy ();

    // drug
    //
    typedef ::xsde::cxx::hybrid::var_sequence< ::pkpd::drug_dose > drug_sequence;
    typedef drug_sequence::iterator drug_iterator;
    typedef drug_sequence::const_iterator drug_const_iterator;

    const drug_sequence&
    drug () const;

    drug_sequence&
    drug ();

    void
    _copy (therapy&) const;

    private:
    drug_sequence drug_;
  };

  // response (variable-length)
  //
  class response
  {
    private:
    response (const response&);
    response& operator= (const response&);

    public:
    response ();

    response*
    _clone () const;

    ~response ();

    // maximum_birth_inhibition
    //
    bool
    maximum_birth_inhibition_present () const;

    const ::common::units_decimal&
    maximum_birth_inhibition () const;

    ::common::units_decimal&
    maximum_birth_inhibition ();

    void
    maximum_birth_inhibition (::common::units_decimal*);

    ::common::units_decimal*
    maximum_birth_inhibition_detach ();

    // maximum_birth_inhibition_time
    //
    bool
    maximum_birth_inhibition_time_present () const;

    const ::common::units_decimal&
    maximum_birth_inhibition_time () const;

    ::common::units_decimal&
    maximum_birth_inhibition_time ();

    void
    maximum_birth_inhibition_time (::common::units_decimal*);

    ::common::units_decimal*
    maximum_birth_inhibition_time_detach ();

    // birth_inhibition_recovery_rate
    //
    bool
    birth_inhibition_recovery_rate_present () const;

    const ::common::units_decimal&
    birth_inhibition_recovery_rate () const;

    ::common::units_decimal&
    birth_inhibition_recovery_rate ();

    void
    birth_inhibition_recovery_rate (::common::units_decimal*);

    ::common::units_decimal*
    birth_inhibition_recovery_rate_detach ();

    // maximum_death_rate
    //
    bool
    maximum_death_rate_present () const;

    const ::common::units_decimal&
    maximum_death_rate () const;

    ::common::units_decimal&
    maximum_death_rate ();

    void
    maximum_death_rate (::common::units_decimal*);

    ::common::units_decimal*
    maximum_death_rate_detach ();

    // maximum_death_time
    //
    bool
    maximum_death_time_present () const;

    const ::common::units_decimal&
    maximum_death_time () const;

    ::common::units_decimal&
    maximum_death_time ();

    void
    maximum_death_time (::common::units_decimal*);

    ::common::units_decimal*
    maximum_death_time_detach ();

    // death_recovery_rate
    //
    bool
    death_recovery_rate_present () const;

    const ::common::units_decimal&
    death_recovery_rate () const;

    ::common::units_decimal&
    death_recovery_rate ();

    void
    death_recovery_rate (::common::units_decimal*);

    ::common::units_decimal*
    death_recovery_rate_detach ();

    // response_observation
    //
    typedef ::xsde::cxx::hybrid::var_sequence< ::pkpd::response_observation > response_observation_sequence;
    typedef response_observation_sequence::iterator response_observation_iterator;
    typedef response_observation_sequence::const_iterator response_observation_const_iterator;

    const response_observation_sequence&
    response_observation () const;

    response_observation_sequence&
    response_observation ();

    // custom
    //
    bool
    custom_present () const;

    const ::common::custom&
    custom () const;

    ::common::custom&
    custom ();

    void
    custom (::common::custom*);

    ::common::custom*
    custom_detach ();

    void
    _copy (response&) const;

    private:
    ::common::units_decimal* maximum_birth_inhibition_;
    ::common::units_decimal* maximum_birth_inhibition_time_;
    ::common::units_decimal* birth_inhibition_recovery_rate_;
    ::common::units_decimal* maximum_death_rate_;
    ::common::units_decimal* maximum_death_time_;
    ::common::units_decimal* death_recovery_rate_;
    response_observation_sequence response_observation_;
    ::common::custom* custom_;
  };

  // response_observation (variable-length)
  //
  class response_observation
  {
    private:
    response_observation (const response_observation&);
    response_observation& operator= (const response_observation&);

    public:
    response_observation ();

    response_observation*
    _clone () const;

    ~response_observation ();

    // time
    //
    bool
    time_present () const;

    const ::common::units_decimal&
    time () const;

    ::common::units_decimal&
    time ();

    void
    time (::common::units_decimal*);

    ::common::units_decimal*
    time_detach ();

    // birth_rate
    //
    bool
    birth_rate_present () const;

    const ::common::units_decimal_nonnegative&
    birth_rate () const;

    ::common::units_decimal_nonnegative&
    birth_rate ();

    void
    birth_rate (::common::units_decimal_nonnegative*);

    ::common::units_decimal_nonnegative*
    birth_rate_detach ();

    // death_rate
    //
    bool
    death_rate_present () const;

    const ::common::units_decimal_nonnegative&
    death_rate () const;

    ::common::units_decimal_nonnegative&
    death_rate ();

    void
    death_rate (::common::units_decimal_nonnegative*);

    ::common::units_decimal_nonnegative*
    death_rate_detach ();

    // net_birth_rate
    //
    bool
    net_birth_rate_present () const;

    const ::common::units_decimal&
    net_birth_rate () const;

    ::common::units_decimal&
    net_birth_rate ();

    void
    net_birth_rate (::common::units_decimal*);

    ::common::units_decimal*
    net_birth_rate_detach ();

    // net_death_rate
    //
    bool
    net_death_rate_present () const;

    const ::common::units_decimal&
    net_death_rate () const;

    ::common::units_decimal&
    net_death_rate ();

    void
    net_death_rate (::common::units_decimal*);

    ::common::units_decimal*
    net_death_rate_detach ();

    // apoptotic_duration
    //
    bool
    apoptotic_duration_present () const;

    const ::common::units_decimal&
    apoptotic_duration () const;

    ::common::units_decimal&
    apoptotic_duration ();

    void
    apoptotic_duration (::common::units_decimal*);

    ::common::units_decimal*
    apoptotic_duration_detach ();

    // percent_cell_viability
    //
    bool
    percent_cell_viability_present () const;

    const ::common::units_decimal&
    percent_cell_viability () const;

    ::common::units_decimal&
    percent_cell_viability ();

    void
    percent_cell_viability (::common::units_decimal*);

    ::common::units_decimal*
    percent_cell_viability_detach ();

    // custom
    //
    bool
    custom_present () const;

    const ::common::custom&
    custom () const;

    ::common::custom&
    custom ();

    void
    custom (::common::custom*);

    ::common::custom*
    custom_detach ();

    void
    _copy (response_observation&) const;

    private:
    ::common::units_decimal* time_;
    ::common::units_decimal_nonnegative* birth_rate_;
    ::common::units_decimal_nonnegative* death_rate_;
    ::common::units_decimal* net_birth_rate_;
    ::common::units_decimal* net_death_rate_;
    ::common::units_decimal* apoptotic_duration_;
    ::common::units_decimal* percent_cell_viability_;
    ::common::custom* custom_;
  };

  // pharmacodynamics (variable-length)
  //
  class pharmacodynamics
  {
    private:
    pharmacodynamics (const pharmacodynamics&);
    pharmacodynamics& operator= (const pharmacodynamics&);

    public:
    pharmacodynamics ();

    pharmacodynamics*
    _clone () const;

    ~pharmacodynamics ();

    // therapy_measurement_set
    //
    typedef ::xsde::cxx::hybrid::var_sequence< ::pkpd::therapy_measurement_set > therapy_measurement_set_sequence;
    typedef therapy_measurement_set_sequence::iterator therapy_measurement_set_iterator;
    typedef therapy_measurement_set_sequence::const_iterator therapy_measurement_set_const_iterator;

    const therapy_measurement_set_sequence&
    therapy_measurement_set () const;

    therapy_measurement_set_sequence&
    therapy_measurement_set ();

    void
    _copy (pharmacodynamics&) const;

    private:
    therapy_measurement_set_sequence therapy_measurement_set_;
  };

  // therapy_measurement_set (variable-length)
  //
  class therapy_measurement_set
  {
    private:
    therapy_measurement_set (const therapy_measurement_set&);
    therapy_measurement_set& operator= (const therapy_measurement_set&);

    public:
    therapy_measurement_set ();

    therapy_measurement_set*
    _clone () const;

    ~therapy_measurement_set ();

    // ID
    //
    bool
    ID_present () const;

    void
    ID_present (bool);

    unsigned short
    ID () const;

    unsigned short&
    ID ();

    void
    ID (unsigned short);

    // therapy
    //
    const ::pkpd::therapy&
    therapy () const;

    ::pkpd::therapy&
    therapy ();

    void
    therapy (::pkpd::therapy*);

    ::pkpd::therapy*
    therapy_detach ();

    // response
    //
    const ::pkpd::response&
    response () const;

    ::pkpd::response&
    response ();

    void
    response (::pkpd::response*);

    ::pkpd::response*
    response_detach ();

    void
    _copy (therapy_measurement_set&) const;

    private:
    unsigned short ID_;
    unsigned char ID_present_;
    ::pkpd::therapy* therapy_;
    ::pkpd::response* response_;
  };

  // PKPD (variable-length)
  //
  class PKPD
  {
    private:
    PKPD (const PKPD&);
    PKPD& operator= (const PKPD&);

    public:
    PKPD ();

    PKPD*
    _clone () const;

    ~PKPD ();

    // drug
    //
    typedef ::xsde::cxx::hybrid::var_sequence< ::pkpd::drug_pk > drug_sequence;
    typedef drug_sequence::iterator drug_iterator;
    typedef drug_sequence::const_iterator drug_const_iterator;

    const drug_sequence&
    drug () const;

    drug_sequence&
    drug ();

    // pharmacodynamics
    //
    bool
    pharmacodynamics_present () const;

    const ::pkpd::pharmacodynamics&
    pharmacodynamics () const;

    ::pkpd::pharmacodynamics&
    pharmacodynamics ();

    void
    pharmacodynamics (::pkpd::pharmacodynamics*);

    ::pkpd::pharmacodynamics*
    pharmacodynamics_detach ();

    void
    _copy (PKPD&) const;

    private:
    drug_sequence drug_;
    ::pkpd::pharmacodynamics* pharmacodynamics_;
  };
}

#ifndef XSDE_DONT_INCLUDE_INLINE

#include "common.ipp"

#include "variables.ipp"

#endif // XSDE_DONT_INCLUDE_INLINE

#ifndef XSDE_DONT_INCLUDE_INLINE
#include "pkpd.ipp"
#endif // XSDE_DONT_INCLUDE_INLINE

// Begin epilogue.
//
//
// End epilogue.

#include <xsde/cxx/post.hxx>

#endif // PKPD_HPP
