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

namespace pkpd
{
  // pharmacokinetics_sskel
  //

  inline
  void pharmacokinetics_sskel::
  inactivation_rate_serializer (::common::units_decimal_sskel& s)
  {
    this->inactivation_rate_serializer_ = &s;
  }

  inline
  void pharmacokinetics_sskel::
  inactivation_rate_serializer (::xml_schema::serializer_map& m)
  {
    this->inactivation_rate_serializer_map_ = &m;
  }

  inline
  void pharmacokinetics_sskel::
  half_life_serializer (::common::units_decimal_sskel& s)
  {
    this->half_life_serializer_ = &s;
  }

  inline
  void pharmacokinetics_sskel::
  half_life_serializer (::xml_schema::serializer_map& m)
  {
    this->half_life_serializer_map_ = &m;
  }

  inline
  void pharmacokinetics_sskel::
  serializers (::common::units_decimal_sskel& inactivation_rate,
               ::common::units_decimal_sskel& half_life)
  {
    this->inactivation_rate_serializer_ = &inactivation_rate;
    this->half_life_serializer_ = &half_life;
  }

  inline
  void pharmacokinetics_sskel::
  serializer_maps (::xml_schema::serializer_map& inactivation_rate,
                   ::xml_schema::serializer_map& half_life)
  {
    this->inactivation_rate_serializer_map_ = &inactivation_rate;
    this->half_life_serializer_map_ = &half_life;
  }

  inline
  pharmacokinetics_sskel::
  pharmacokinetics_sskel ()
  : pharmacokinetics_impl_ (0),
    inactivation_rate_serializer_ (0),
    inactivation_rate_serializer_map_ (0),
    half_life_serializer_ (0),
    half_life_serializer_map_ (0)
  {
  }

  inline
  pharmacokinetics_sskel::
  pharmacokinetics_sskel (pharmacokinetics_sskel* impl, void*)
  : ::xsde::cxx::serializer::validating::complex_content (impl, 0),
    pharmacokinetics_impl_ (impl),
    inactivation_rate_serializer_ (0),
    inactivation_rate_serializer_map_ (0),
    half_life_serializer_ (0),
    half_life_serializer_map_ (0)
  {
  }

  // drug_sskel
  //

  inline
  void drug_sskel::
  dose_serializer (::pkpd::dose_sskel& s)
  {
    this->dose_serializer_ = &s;
  }

  inline
  void drug_sskel::
  dose_serializer (::xml_schema::serializer_map& m)
  {
    this->dose_serializer_map_ = &m;
  }

  inline
  void drug_sskel::
  pharmacokinetics_serializer (::pkpd::pharmacokinetics_sskel& s)
  {
    this->pharmacokinetics_serializer_ = &s;
  }

  inline
  void drug_sskel::
  pharmacokinetics_serializer (::xml_schema::serializer_map& m)
  {
    this->pharmacokinetics_serializer_map_ = &m;
  }

  inline
  void drug_sskel::
  serializers (::pkpd::dose_sskel& dose,
               ::pkpd::pharmacokinetics_sskel& pharmacokinetics)
  {
    this->dose_serializer_ = &dose;
    this->pharmacokinetics_serializer_ = &pharmacokinetics;
  }

  inline
  void drug_sskel::
  serializer_maps (::xml_schema::serializer_map& dose,
                   ::xml_schema::serializer_map& pharmacokinetics)
  {
    this->dose_serializer_map_ = &dose;
    this->pharmacokinetics_serializer_map_ = &pharmacokinetics;
  }

  inline
  drug_sskel::
  drug_sskel ()
  : drug_impl_ (0),
    dose_serializer_ (0),
    dose_serializer_map_ (0),
    pharmacokinetics_serializer_ (0),
    pharmacokinetics_serializer_map_ (0)
  {
  }

  inline
  drug_sskel::
  drug_sskel (drug_sskel* impl, void*)
  : ::xsde::cxx::serializer::validating::complex_content (impl, 0),
    drug_impl_ (impl),
    dose_serializer_ (0),
    dose_serializer_map_ (0),
    pharmacokinetics_serializer_ (0),
    pharmacokinetics_serializer_map_ (0)
  {
  }

  // drug_dose_sskel
  //

  inline
  void drug_dose_sskel::
  ID_serializer (::xml_schema::unsigned_short_sskel& ID)
  {
    this->ID_serializer_ = &ID;
  }

  inline
  void drug_dose_sskel::
  ChEBI_ID_serializer (::xml_schema::string_sskel& ChEBI_ID)
  {
    this->ChEBI_ID_serializer_ = &ChEBI_ID;
  }

  inline
  void drug_dose_sskel::
  MeSH_ID_serializer (::xml_schema::string_sskel& MeSH_ID)
  {
    this->MeSH_ID_serializer_ = &MeSH_ID;
  }

  inline
  void drug_dose_sskel::
  DrugBank_ID_serializer (::xml_schema::string_sskel& DrugBank_ID)
  {
    this->DrugBank_ID_serializer_ = &DrugBank_ID;
  }

  inline
  void drug_dose_sskel::
  GMO_ID_serializer (::xml_schema::string_sskel& GMO_ID)
  {
    this->GMO_ID_serializer_ = &GMO_ID;
  }

  inline
  void drug_dose_sskel::
  GO_ID_serializer (::xml_schema::string_sskel& GO_ID)
  {
    this->GO_ID_serializer_ = &GO_ID;
  }

  inline
  void drug_dose_sskel::
  UniProt_ID_serializer (::xml_schema::string_sskel& UniProt_ID)
  {
    this->UniProt_ID_serializer_ = &UniProt_ID;
  }

  inline
  void drug_dose_sskel::
  PR_ID_serializer (::xml_schema::string_sskel& PR_ID)
  {
    this->PR_ID_serializer_ = &PR_ID;
  }

  inline
  void drug_dose_sskel::
  name_serializer (::xml_schema::string_sskel& name)
  {
    this->name_serializer_ = &name;
  }

  inline
  void drug_dose_sskel::
  units_serializer (::xml_schema::string_sskel& units)
  {
    this->units_serializer_ = &units;
  }

  inline
  void drug_dose_sskel::
  dose_serializer (::pkpd::dose_sskel& s)
  {
    this->dose_serializer_ = &s;
  }

  inline
  void drug_dose_sskel::
  dose_serializer (::xml_schema::serializer_map& m)
  {
    this->dose_serializer_map_ = &m;
  }

  inline
  void drug_dose_sskel::
  serializers (::xml_schema::unsigned_short_sskel& ID,
               ::xml_schema::string_sskel& ChEBI_ID,
               ::xml_schema::string_sskel& MeSH_ID,
               ::xml_schema::string_sskel& DrugBank_ID,
               ::xml_schema::string_sskel& GMO_ID,
               ::xml_schema::string_sskel& GO_ID,
               ::xml_schema::string_sskel& UniProt_ID,
               ::xml_schema::string_sskel& PR_ID,
               ::xml_schema::string_sskel& name,
               ::xml_schema::string_sskel& units,
               ::pkpd::dose_sskel& dose)
  {
    this->ID_serializer_ = &ID;
    this->ChEBI_ID_serializer_ = &ChEBI_ID;
    this->MeSH_ID_serializer_ = &MeSH_ID;
    this->DrugBank_ID_serializer_ = &DrugBank_ID;
    this->GMO_ID_serializer_ = &GMO_ID;
    this->GO_ID_serializer_ = &GO_ID;
    this->UniProt_ID_serializer_ = &UniProt_ID;
    this->PR_ID_serializer_ = &PR_ID;
    this->name_serializer_ = &name;
    this->units_serializer_ = &units;
    this->dose_serializer_ = &dose;
  }

  inline
  void drug_dose_sskel::
  serializer_maps (::xml_schema::serializer_map& dose)
  {
    this->dose_serializer_map_ = &dose;
  }

  inline
  drug_dose_sskel::
  drug_dose_sskel ()
  : drug_dose_impl_ (0),
    ID_serializer_ (0),
    ChEBI_ID_serializer_ (0),
    MeSH_ID_serializer_ (0),
    DrugBank_ID_serializer_ (0),
    GMO_ID_serializer_ (0),
    GO_ID_serializer_ (0),
    UniProt_ID_serializer_ (0),
    PR_ID_serializer_ (0),
    name_serializer_ (0),
    units_serializer_ (0),
    dose_serializer_ (0),
    dose_serializer_map_ (0)
  {
  }

  inline
  drug_dose_sskel::
  drug_dose_sskel (drug_dose_sskel* impl, void*)
  : ::xsde::cxx::serializer::validating::complex_content (impl, 0),
    drug_dose_impl_ (impl),
    ID_serializer_ (0),
    ChEBI_ID_serializer_ (0),
    MeSH_ID_serializer_ (0),
    DrugBank_ID_serializer_ (0),
    GMO_ID_serializer_ (0),
    GO_ID_serializer_ (0),
    UniProt_ID_serializer_ (0),
    PR_ID_serializer_ (0),
    name_serializer_ (0),
    units_serializer_ (0),
    dose_serializer_ (0),
    dose_serializer_map_ (0)
  {
  }

  // drug_pk_sskel
  //

  inline
  void drug_pk_sskel::
  ID_serializer (::xml_schema::unsigned_short_sskel& ID)
  {
    this->ID_serializer_ = &ID;
  }

  inline
  void drug_pk_sskel::
  ChEBI_ID_serializer (::xml_schema::string_sskel& ChEBI_ID)
  {
    this->ChEBI_ID_serializer_ = &ChEBI_ID;
  }

  inline
  void drug_pk_sskel::
  MeSH_ID_serializer (::xml_schema::string_sskel& MeSH_ID)
  {
    this->MeSH_ID_serializer_ = &MeSH_ID;
  }

  inline
  void drug_pk_sskel::
  DrugBank_ID_serializer (::xml_schema::string_sskel& DrugBank_ID)
  {
    this->DrugBank_ID_serializer_ = &DrugBank_ID;
  }

  inline
  void drug_pk_sskel::
  GMO_ID_serializer (::xml_schema::string_sskel& GMO_ID)
  {
    this->GMO_ID_serializer_ = &GMO_ID;
  }

  inline
  void drug_pk_sskel::
  GO_ID_serializer (::xml_schema::string_sskel& GO_ID)
  {
    this->GO_ID_serializer_ = &GO_ID;
  }

  inline
  void drug_pk_sskel::
  UniProt_ID_serializer (::xml_schema::string_sskel& UniProt_ID)
  {
    this->UniProt_ID_serializer_ = &UniProt_ID;
  }

  inline
  void drug_pk_sskel::
  PR_ID_serializer (::xml_schema::string_sskel& PR_ID)
  {
    this->PR_ID_serializer_ = &PR_ID;
  }

  inline
  void drug_pk_sskel::
  name_serializer (::xml_schema::string_sskel& name)
  {
    this->name_serializer_ = &name;
  }

  inline
  void drug_pk_sskel::
  units_serializer (::xml_schema::string_sskel& units)
  {
    this->units_serializer_ = &units;
  }

  inline
  void drug_pk_sskel::
  pharmacokinetics_serializer (::pkpd::pharmacokinetics_sskel& s)
  {
    this->pharmacokinetics_serializer_ = &s;
  }

  inline
  void drug_pk_sskel::
  pharmacokinetics_serializer (::xml_schema::serializer_map& m)
  {
    this->pharmacokinetics_serializer_map_ = &m;
  }

  inline
  void drug_pk_sskel::
  serializers (::xml_schema::unsigned_short_sskel& ID,
               ::xml_schema::string_sskel& ChEBI_ID,
               ::xml_schema::string_sskel& MeSH_ID,
               ::xml_schema::string_sskel& DrugBank_ID,
               ::xml_schema::string_sskel& GMO_ID,
               ::xml_schema::string_sskel& GO_ID,
               ::xml_schema::string_sskel& UniProt_ID,
               ::xml_schema::string_sskel& PR_ID,
               ::xml_schema::string_sskel& name,
               ::xml_schema::string_sskel& units,
               ::pkpd::pharmacokinetics_sskel& pharmacokinetics)
  {
    this->ID_serializer_ = &ID;
    this->ChEBI_ID_serializer_ = &ChEBI_ID;
    this->MeSH_ID_serializer_ = &MeSH_ID;
    this->DrugBank_ID_serializer_ = &DrugBank_ID;
    this->GMO_ID_serializer_ = &GMO_ID;
    this->GO_ID_serializer_ = &GO_ID;
    this->UniProt_ID_serializer_ = &UniProt_ID;
    this->PR_ID_serializer_ = &PR_ID;
    this->name_serializer_ = &name;
    this->units_serializer_ = &units;
    this->pharmacokinetics_serializer_ = &pharmacokinetics;
  }

  inline
  void drug_pk_sskel::
  serializer_maps (::xml_schema::serializer_map& pharmacokinetics)
  {
    this->pharmacokinetics_serializer_map_ = &pharmacokinetics;
  }

  inline
  drug_pk_sskel::
  drug_pk_sskel ()
  : drug_pk_impl_ (0),
    ID_serializer_ (0),
    ChEBI_ID_serializer_ (0),
    MeSH_ID_serializer_ (0),
    DrugBank_ID_serializer_ (0),
    GMO_ID_serializer_ (0),
    GO_ID_serializer_ (0),
    UniProt_ID_serializer_ (0),
    PR_ID_serializer_ (0),
    name_serializer_ (0),
    units_serializer_ (0),
    pharmacokinetics_serializer_ (0),
    pharmacokinetics_serializer_map_ (0)
  {
  }

  inline
  drug_pk_sskel::
  drug_pk_sskel (drug_pk_sskel* impl, void*)
  : ::xsde::cxx::serializer::validating::complex_content (impl, 0),
    drug_pk_impl_ (impl),
    ID_serializer_ (0),
    ChEBI_ID_serializer_ (0),
    MeSH_ID_serializer_ (0),
    DrugBank_ID_serializer_ (0),
    GMO_ID_serializer_ (0),
    GO_ID_serializer_ (0),
    UniProt_ID_serializer_ (0),
    PR_ID_serializer_ (0),
    name_serializer_ (0),
    units_serializer_ (0),
    pharmacokinetics_serializer_ (0),
    pharmacokinetics_serializer_map_ (0)
  {
  }

  // dose_sskel
  //

  inline
  void dose_sskel::
  type_serializer (::xml_schema::string_sskel& type)
  {
    this->type_serializer_ = &type;
  }

  inline
  void dose_sskel::
  serializers (::xml_schema::string_sskel& units,
               ::xml_schema::string_sskel& measurement_type,
               ::xml_schema::double_sskel& uncertainty,
               ::xml_schema::double_sskel& negative_uncertainty,
               ::xml_schema::double_sskel& positive_uncertainty,
               ::xml_schema::double_sskel& uncertainty_percentage,
               ::xml_schema::double_sskel& negative_uncertainty_percentage,
               ::xml_schema::double_sskel& positive_uncertainty_percentage,
               ::xml_schema::double_sskel& median,
               ::xml_schema::double_sskel& standard_deviation,
               ::common::two_doubles_sskel& interquartile_range,
               ::common::two_doubles_sskel& range,
               ::xml_schema::double_sskel& min,
               ::xml_schema::double_sskel& max,
               ::xml_schema::double_sskel& standard_error,
               ::xml_schema::double_sskel& standard_error_of_the_mean,
               ::xml_schema::int_sskel& number_obs,
               ::xml_schema::double_sskel& skewnesss,
               ::xml_schema::double_sskel& kurtosis,
               ::xml_schema::string_sskel& type)
  {
    this->units_serializer_ = &units;
    this->measurement_type_serializer_ = &measurement_type;
    this->uncertainty_serializer_ = &uncertainty;
    this->negative_uncertainty_serializer_ = &negative_uncertainty;
    this->positive_uncertainty_serializer_ = &positive_uncertainty;
    this->uncertainty_percentage_serializer_ = &uncertainty_percentage;
    this->negative_uncertainty_percentage_serializer_ = &negative_uncertainty_percentage;
    this->positive_uncertainty_percentage_serializer_ = &positive_uncertainty_percentage;
    this->median_serializer_ = &median;
    this->standard_deviation_serializer_ = &standard_deviation;
    this->interquartile_range_serializer_ = &interquartile_range;
    this->range_serializer_ = &range;
    this->min_serializer_ = &min;
    this->max_serializer_ = &max;
    this->standard_error_serializer_ = &standard_error;
    this->standard_error_of_the_mean_serializer_ = &standard_error_of_the_mean;
    this->number_obs_serializer_ = &number_obs;
    this->skewnesss_serializer_ = &skewnesss;
    this->kurtosis_serializer_ = &kurtosis;
    this->type_serializer_ = &type;
  }

  inline
  dose_sskel::
  dose_sskel (::common::units_decimal_sskel* tiein)
  : ::common::units_decimal_sskel (tiein, 0),
    dose_impl_ (0),
    type_serializer_ (0)
  {
  }

  inline
  dose_sskel::
  dose_sskel (dose_sskel* impl, void*)
  : ::common::units_decimal_sskel (impl, 0),
    dose_impl_ (impl),
    type_serializer_ (0)
  {
  }

  // therapy_sskel
  //

  inline
  void therapy_sskel::
  drug_serializer (::pkpd::drug_dose_sskel& s)
  {
    this->drug_serializer_ = &s;
  }

  inline
  void therapy_sskel::
  drug_serializer (::xml_schema::serializer_map& m)
  {
    this->drug_serializer_map_ = &m;
  }

  inline
  void therapy_sskel::
  serializers (::pkpd::drug_dose_sskel& drug)
  {
    this->drug_serializer_ = &drug;
  }

  inline
  void therapy_sskel::
  serializer_maps (::xml_schema::serializer_map& drug)
  {
    this->drug_serializer_map_ = &drug;
  }

  inline
  therapy_sskel::
  therapy_sskel ()
  : therapy_impl_ (0),
    drug_serializer_ (0),
    drug_serializer_map_ (0)
  {
  }

  inline
  therapy_sskel::
  therapy_sskel (therapy_sskel* impl, void*)
  : ::xsde::cxx::serializer::validating::complex_content (impl, 0),
    therapy_impl_ (impl),
    drug_serializer_ (0),
    drug_serializer_map_ (0)
  {
  }

  // response_sskel
  //

  inline
  void response_sskel::
  maximum_birth_inhibition_serializer (::common::units_decimal_sskel& s)
  {
    this->maximum_birth_inhibition_serializer_ = &s;
  }

  inline
  void response_sskel::
  maximum_birth_inhibition_serializer (::xml_schema::serializer_map& m)
  {
    this->maximum_birth_inhibition_serializer_map_ = &m;
  }

  inline
  void response_sskel::
  maximum_birth_inhibition_time_serializer (::common::units_decimal_sskel& s)
  {
    this->maximum_birth_inhibition_time_serializer_ = &s;
  }

  inline
  void response_sskel::
  maximum_birth_inhibition_time_serializer (::xml_schema::serializer_map& m)
  {
    this->maximum_birth_inhibition_time_serializer_map_ = &m;
  }

  inline
  void response_sskel::
  birth_inhibition_recovery_rate_serializer (::common::units_decimal_sskel& s)
  {
    this->birth_inhibition_recovery_rate_serializer_ = &s;
  }

  inline
  void response_sskel::
  birth_inhibition_recovery_rate_serializer (::xml_schema::serializer_map& m)
  {
    this->birth_inhibition_recovery_rate_serializer_map_ = &m;
  }

  inline
  void response_sskel::
  maximum_death_rate_serializer (::common::units_decimal_sskel& s)
  {
    this->maximum_death_rate_serializer_ = &s;
  }

  inline
  void response_sskel::
  maximum_death_rate_serializer (::xml_schema::serializer_map& m)
  {
    this->maximum_death_rate_serializer_map_ = &m;
  }

  inline
  void response_sskel::
  maximum_death_time_serializer (::common::units_decimal_sskel& s)
  {
    this->maximum_death_time_serializer_ = &s;
  }

  inline
  void response_sskel::
  maximum_death_time_serializer (::xml_schema::serializer_map& m)
  {
    this->maximum_death_time_serializer_map_ = &m;
  }

  inline
  void response_sskel::
  death_recovery_rate_serializer (::common::units_decimal_sskel& s)
  {
    this->death_recovery_rate_serializer_ = &s;
  }

  inline
  void response_sskel::
  death_recovery_rate_serializer (::xml_schema::serializer_map& m)
  {
    this->death_recovery_rate_serializer_map_ = &m;
  }

  inline
  void response_sskel::
  response_observation_serializer (::pkpd::response_observation_sskel& s)
  {
    this->response_observation_serializer_ = &s;
  }

  inline
  void response_sskel::
  response_observation_serializer (::xml_schema::serializer_map& m)
  {
    this->response_observation_serializer_map_ = &m;
  }

  inline
  void response_sskel::
  custom_serializer (::common::custom_sskel& s)
  {
    this->custom_serializer_ = &s;
  }

  inline
  void response_sskel::
  custom_serializer (::xml_schema::serializer_map& m)
  {
    this->custom_serializer_map_ = &m;
  }

  inline
  void response_sskel::
  serializers (::common::units_decimal_sskel& maximum_birth_inhibition,
               ::common::units_decimal_sskel& maximum_birth_inhibition_time,
               ::common::units_decimal_sskel& birth_inhibition_recovery_rate,
               ::common::units_decimal_sskel& maximum_death_rate,
               ::common::units_decimal_sskel& maximum_death_time,
               ::common::units_decimal_sskel& death_recovery_rate,
               ::pkpd::response_observation_sskel& response_observation,
               ::common::custom_sskel& custom)
  {
    this->maximum_birth_inhibition_serializer_ = &maximum_birth_inhibition;
    this->maximum_birth_inhibition_time_serializer_ = &maximum_birth_inhibition_time;
    this->birth_inhibition_recovery_rate_serializer_ = &birth_inhibition_recovery_rate;
    this->maximum_death_rate_serializer_ = &maximum_death_rate;
    this->maximum_death_time_serializer_ = &maximum_death_time;
    this->death_recovery_rate_serializer_ = &death_recovery_rate;
    this->response_observation_serializer_ = &response_observation;
    this->custom_serializer_ = &custom;
  }

  inline
  void response_sskel::
  serializer_maps (::xml_schema::serializer_map& maximum_birth_inhibition,
                   ::xml_schema::serializer_map& maximum_birth_inhibition_time,
                   ::xml_schema::serializer_map& birth_inhibition_recovery_rate,
                   ::xml_schema::serializer_map& maximum_death_rate,
                   ::xml_schema::serializer_map& maximum_death_time,
                   ::xml_schema::serializer_map& death_recovery_rate,
                   ::xml_schema::serializer_map& response_observation,
                   ::xml_schema::serializer_map& custom)
  {
    this->maximum_birth_inhibition_serializer_map_ = &maximum_birth_inhibition;
    this->maximum_birth_inhibition_time_serializer_map_ = &maximum_birth_inhibition_time;
    this->birth_inhibition_recovery_rate_serializer_map_ = &birth_inhibition_recovery_rate;
    this->maximum_death_rate_serializer_map_ = &maximum_death_rate;
    this->maximum_death_time_serializer_map_ = &maximum_death_time;
    this->death_recovery_rate_serializer_map_ = &death_recovery_rate;
    this->response_observation_serializer_map_ = &response_observation;
    this->custom_serializer_map_ = &custom;
  }

  inline
  response_sskel::
  response_sskel ()
  : response_impl_ (0),
    maximum_birth_inhibition_serializer_ (0),
    maximum_birth_inhibition_serializer_map_ (0),
    maximum_birth_inhibition_time_serializer_ (0),
    maximum_birth_inhibition_time_serializer_map_ (0),
    birth_inhibition_recovery_rate_serializer_ (0),
    birth_inhibition_recovery_rate_serializer_map_ (0),
    maximum_death_rate_serializer_ (0),
    maximum_death_rate_serializer_map_ (0),
    maximum_death_time_serializer_ (0),
    maximum_death_time_serializer_map_ (0),
    death_recovery_rate_serializer_ (0),
    death_recovery_rate_serializer_map_ (0),
    response_observation_serializer_ (0),
    response_observation_serializer_map_ (0),
    custom_serializer_ (0),
    custom_serializer_map_ (0)
  {
  }

  inline
  response_sskel::
  response_sskel (response_sskel* impl, void*)
  : ::xsde::cxx::serializer::validating::complex_content (impl, 0),
    response_impl_ (impl),
    maximum_birth_inhibition_serializer_ (0),
    maximum_birth_inhibition_serializer_map_ (0),
    maximum_birth_inhibition_time_serializer_ (0),
    maximum_birth_inhibition_time_serializer_map_ (0),
    birth_inhibition_recovery_rate_serializer_ (0),
    birth_inhibition_recovery_rate_serializer_map_ (0),
    maximum_death_rate_serializer_ (0),
    maximum_death_rate_serializer_map_ (0),
    maximum_death_time_serializer_ (0),
    maximum_death_time_serializer_map_ (0),
    death_recovery_rate_serializer_ (0),
    death_recovery_rate_serializer_map_ (0),
    response_observation_serializer_ (0),
    response_observation_serializer_map_ (0),
    custom_serializer_ (0),
    custom_serializer_map_ (0)
  {
  }

  // response_observation_sskel
  //

  inline
  void response_observation_sskel::
  time_serializer (::common::units_decimal_sskel& s)
  {
    this->time_serializer_ = &s;
  }

  inline
  void response_observation_sskel::
  time_serializer (::xml_schema::serializer_map& m)
  {
    this->time_serializer_map_ = &m;
  }

  inline
  void response_observation_sskel::
  birth_rate_serializer (::common::units_decimal_nonnegative_sskel& s)
  {
    this->birth_rate_serializer_ = &s;
  }

  inline
  void response_observation_sskel::
  birth_rate_serializer (::xml_schema::serializer_map& m)
  {
    this->birth_rate_serializer_map_ = &m;
  }

  inline
  void response_observation_sskel::
  death_rate_serializer (::common::units_decimal_nonnegative_sskel& s)
  {
    this->death_rate_serializer_ = &s;
  }

  inline
  void response_observation_sskel::
  death_rate_serializer (::xml_schema::serializer_map& m)
  {
    this->death_rate_serializer_map_ = &m;
  }

  inline
  void response_observation_sskel::
  net_birth_rate_serializer (::common::units_decimal_sskel& s)
  {
    this->net_birth_rate_serializer_ = &s;
  }

  inline
  void response_observation_sskel::
  net_birth_rate_serializer (::xml_schema::serializer_map& m)
  {
    this->net_birth_rate_serializer_map_ = &m;
  }

  inline
  void response_observation_sskel::
  net_death_rate_serializer (::common::units_decimal_sskel& s)
  {
    this->net_death_rate_serializer_ = &s;
  }

  inline
  void response_observation_sskel::
  net_death_rate_serializer (::xml_schema::serializer_map& m)
  {
    this->net_death_rate_serializer_map_ = &m;
  }

  inline
  void response_observation_sskel::
  apoptotic_duration_serializer (::common::units_decimal_sskel& s)
  {
    this->apoptotic_duration_serializer_ = &s;
  }

  inline
  void response_observation_sskel::
  apoptotic_duration_serializer (::xml_schema::serializer_map& m)
  {
    this->apoptotic_duration_serializer_map_ = &m;
  }

  inline
  void response_observation_sskel::
  percent_cell_viability_serializer (::common::units_decimal_sskel& s)
  {
    this->percent_cell_viability_serializer_ = &s;
  }

  inline
  void response_observation_sskel::
  percent_cell_viability_serializer (::xml_schema::serializer_map& m)
  {
    this->percent_cell_viability_serializer_map_ = &m;
  }

  inline
  void response_observation_sskel::
  custom_serializer (::common::custom_sskel& s)
  {
    this->custom_serializer_ = &s;
  }

  inline
  void response_observation_sskel::
  custom_serializer (::xml_schema::serializer_map& m)
  {
    this->custom_serializer_map_ = &m;
  }

  inline
  void response_observation_sskel::
  serializers (::common::units_decimal_sskel& time,
               ::common::units_decimal_nonnegative_sskel& birth_rate,
               ::common::units_decimal_nonnegative_sskel& death_rate,
               ::common::units_decimal_sskel& net_birth_rate,
               ::common::units_decimal_sskel& net_death_rate,
               ::common::units_decimal_sskel& apoptotic_duration,
               ::common::units_decimal_sskel& percent_cell_viability,
               ::common::custom_sskel& custom)
  {
    this->time_serializer_ = &time;
    this->birth_rate_serializer_ = &birth_rate;
    this->death_rate_serializer_ = &death_rate;
    this->net_birth_rate_serializer_ = &net_birth_rate;
    this->net_death_rate_serializer_ = &net_death_rate;
    this->apoptotic_duration_serializer_ = &apoptotic_duration;
    this->percent_cell_viability_serializer_ = &percent_cell_viability;
    this->custom_serializer_ = &custom;
  }

  inline
  void response_observation_sskel::
  serializer_maps (::xml_schema::serializer_map& time,
                   ::xml_schema::serializer_map& birth_rate,
                   ::xml_schema::serializer_map& death_rate,
                   ::xml_schema::serializer_map& net_birth_rate,
                   ::xml_schema::serializer_map& net_death_rate,
                   ::xml_schema::serializer_map& apoptotic_duration,
                   ::xml_schema::serializer_map& percent_cell_viability,
                   ::xml_schema::serializer_map& custom)
  {
    this->time_serializer_map_ = &time;
    this->birth_rate_serializer_map_ = &birth_rate;
    this->death_rate_serializer_map_ = &death_rate;
    this->net_birth_rate_serializer_map_ = &net_birth_rate;
    this->net_death_rate_serializer_map_ = &net_death_rate;
    this->apoptotic_duration_serializer_map_ = &apoptotic_duration;
    this->percent_cell_viability_serializer_map_ = &percent_cell_viability;
    this->custom_serializer_map_ = &custom;
  }

  inline
  response_observation_sskel::
  response_observation_sskel ()
  : response_observation_impl_ (0),
    time_serializer_ (0),
    time_serializer_map_ (0),
    birth_rate_serializer_ (0),
    birth_rate_serializer_map_ (0),
    death_rate_serializer_ (0),
    death_rate_serializer_map_ (0),
    net_birth_rate_serializer_ (0),
    net_birth_rate_serializer_map_ (0),
    net_death_rate_serializer_ (0),
    net_death_rate_serializer_map_ (0),
    apoptotic_duration_serializer_ (0),
    apoptotic_duration_serializer_map_ (0),
    percent_cell_viability_serializer_ (0),
    percent_cell_viability_serializer_map_ (0),
    custom_serializer_ (0),
    custom_serializer_map_ (0)
  {
  }

  inline
  response_observation_sskel::
  response_observation_sskel (response_observation_sskel* impl, void*)
  : ::xsde::cxx::serializer::validating::complex_content (impl, 0),
    response_observation_impl_ (impl),
    time_serializer_ (0),
    time_serializer_map_ (0),
    birth_rate_serializer_ (0),
    birth_rate_serializer_map_ (0),
    death_rate_serializer_ (0),
    death_rate_serializer_map_ (0),
    net_birth_rate_serializer_ (0),
    net_birth_rate_serializer_map_ (0),
    net_death_rate_serializer_ (0),
    net_death_rate_serializer_map_ (0),
    apoptotic_duration_serializer_ (0),
    apoptotic_duration_serializer_map_ (0),
    percent_cell_viability_serializer_ (0),
    percent_cell_viability_serializer_map_ (0),
    custom_serializer_ (0),
    custom_serializer_map_ (0)
  {
  }

  // pharmacodynamics_sskel
  //

  inline
  void pharmacodynamics_sskel::
  therapy_measurement_set_serializer (::pkpd::therapy_measurement_set_sskel& s)
  {
    this->therapy_measurement_set_serializer_ = &s;
  }

  inline
  void pharmacodynamics_sskel::
  therapy_measurement_set_serializer (::xml_schema::serializer_map& m)
  {
    this->therapy_measurement_set_serializer_map_ = &m;
  }

  inline
  void pharmacodynamics_sskel::
  serializers (::pkpd::therapy_measurement_set_sskel& therapy_measurement_set)
  {
    this->therapy_measurement_set_serializer_ = &therapy_measurement_set;
  }

  inline
  void pharmacodynamics_sskel::
  serializer_maps (::xml_schema::serializer_map& therapy_measurement_set)
  {
    this->therapy_measurement_set_serializer_map_ = &therapy_measurement_set;
  }

  inline
  pharmacodynamics_sskel::
  pharmacodynamics_sskel ()
  : pharmacodynamics_impl_ (0),
    therapy_measurement_set_serializer_ (0),
    therapy_measurement_set_serializer_map_ (0)
  {
  }

  inline
  pharmacodynamics_sskel::
  pharmacodynamics_sskel (pharmacodynamics_sskel* impl, void*)
  : ::xsde::cxx::serializer::validating::complex_content (impl, 0),
    pharmacodynamics_impl_ (impl),
    therapy_measurement_set_serializer_ (0),
    therapy_measurement_set_serializer_map_ (0)
  {
  }

  // therapy_measurement_set_sskel
  //

  inline
  void therapy_measurement_set_sskel::
  ID_serializer (::xml_schema::unsigned_short_sskel& ID)
  {
    this->ID_serializer_ = &ID;
  }

  inline
  void therapy_measurement_set_sskel::
  therapy_serializer (::pkpd::therapy_sskel& s)
  {
    this->therapy_serializer_ = &s;
  }

  inline
  void therapy_measurement_set_sskel::
  therapy_serializer (::xml_schema::serializer_map& m)
  {
    this->therapy_serializer_map_ = &m;
  }

  inline
  void therapy_measurement_set_sskel::
  response_serializer (::pkpd::response_sskel& s)
  {
    this->response_serializer_ = &s;
  }

  inline
  void therapy_measurement_set_sskel::
  response_serializer (::xml_schema::serializer_map& m)
  {
    this->response_serializer_map_ = &m;
  }

  inline
  void therapy_measurement_set_sskel::
  serializers (::xml_schema::unsigned_short_sskel& ID,
               ::pkpd::therapy_sskel& therapy,
               ::pkpd::response_sskel& response)
  {
    this->ID_serializer_ = &ID;
    this->therapy_serializer_ = &therapy;
    this->response_serializer_ = &response;
  }

  inline
  void therapy_measurement_set_sskel::
  serializer_maps (::xml_schema::serializer_map& therapy,
                   ::xml_schema::serializer_map& response)
  {
    this->therapy_serializer_map_ = &therapy;
    this->response_serializer_map_ = &response;
  }

  inline
  therapy_measurement_set_sskel::
  therapy_measurement_set_sskel ()
  : therapy_measurement_set_impl_ (0),
    ID_serializer_ (0),
    therapy_serializer_ (0),
    therapy_serializer_map_ (0),
    response_serializer_ (0),
    response_serializer_map_ (0)
  {
  }

  inline
  therapy_measurement_set_sskel::
  therapy_measurement_set_sskel (therapy_measurement_set_sskel* impl, void*)
  : ::xsde::cxx::serializer::validating::complex_content (impl, 0),
    therapy_measurement_set_impl_ (impl),
    ID_serializer_ (0),
    therapy_serializer_ (0),
    therapy_serializer_map_ (0),
    response_serializer_ (0),
    response_serializer_map_ (0)
  {
  }

  // PKPD_sskel
  //

  inline
  void PKPD_sskel::
  drug_serializer (::pkpd::drug_pk_sskel& s)
  {
    this->drug_serializer_ = &s;
  }

  inline
  void PKPD_sskel::
  drug_serializer (::xml_schema::serializer_map& m)
  {
    this->drug_serializer_map_ = &m;
  }

  inline
  void PKPD_sskel::
  pharmacodynamics_serializer (::pkpd::pharmacodynamics_sskel& s)
  {
    this->pharmacodynamics_serializer_ = &s;
  }

  inline
  void PKPD_sskel::
  pharmacodynamics_serializer (::xml_schema::serializer_map& m)
  {
    this->pharmacodynamics_serializer_map_ = &m;
  }

  inline
  void PKPD_sskel::
  serializers (::pkpd::drug_pk_sskel& drug,
               ::pkpd::pharmacodynamics_sskel& pharmacodynamics)
  {
    this->drug_serializer_ = &drug;
    this->pharmacodynamics_serializer_ = &pharmacodynamics;
  }

  inline
  void PKPD_sskel::
  serializer_maps (::xml_schema::serializer_map& drug,
                   ::xml_schema::serializer_map& pharmacodynamics)
  {
    this->drug_serializer_map_ = &drug;
    this->pharmacodynamics_serializer_map_ = &pharmacodynamics;
  }

  inline
  PKPD_sskel::
  PKPD_sskel ()
  : PKPD_impl_ (0),
    drug_serializer_ (0),
    drug_serializer_map_ (0),
    pharmacodynamics_serializer_ (0),
    pharmacodynamics_serializer_map_ (0)
  {
  }

  inline
  PKPD_sskel::
  PKPD_sskel (PKPD_sskel* impl, void*)
  : ::xsde::cxx::serializer::validating::complex_content (impl, 0),
    PKPD_impl_ (impl),
    drug_serializer_ (0),
    drug_serializer_map_ (0),
    pharmacodynamics_serializer_ (0),
    pharmacodynamics_serializer_map_ (0)
  {
  }
}

// Begin epilogue.
//
//
// End epilogue.

