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

#include "phenotype_base-sskel.hpp"

#include <assert.h>

#include <string.h>
#include <xsde/cxx/serializer/substitution-map.hxx>
#include <xsde/cxx/serializer/validating/inheritance-map.hxx>

static
const ::xsde::cxx::serializer::substitution_map_init
_xsde_substitution_map_init_;

static
const ::xsde::cxx::serializer::validating::inheritance_map_init
_xsde_inheritance_map_init_;

namespace phenotype_base
{
  // phenotype_type_sskel
  //

  const char* phenotype_type_sskel::
  _static_type ()
  {
    return "phenotype_type phenotype_base";
  }

  const char* phenotype_type_sskel::
  _dynamic_type () const
  {
    return _static_type ();
  }

  static
  const ::xsde::cxx::serializer::validating::inheritance_map_entry
  _xsde_phenotype_type_sskel_inheritance_map_entry_ (
    phenotype_type_sskel::_static_type (),
    ::xml_schema::string_sskel::_static_type ());

  void phenotype_type_sskel::
  pre (const ::std::string& x)
  {
    assert (this->string_impl_);
    this->string_impl_->pre (x);
  }

  const char* const phenotype_type_sskel::_xsde_phenotype_type_sskel_enums_[3UL] = 
  {
    "current",
    "expected",
    "target"
  };

  // phenotype_base_sskel
  //

  bool phenotype_base_sskel::
  type_present ()
  {
    return this->phenotype_base_impl_ ? this->phenotype_base_impl_->type_present () : false;
  }

  bool phenotype_base_sskel::
  adhesion_present ()
  {
    return this->phenotype_base_impl_ ? this->phenotype_base_impl_->adhesion_present () : false;
  }

  bool phenotype_base_sskel::
  geometrical_properties_present ()
  {
    return this->phenotype_base_impl_ ? this->phenotype_base_impl_->geometrical_properties_present () : false;
  }

  bool phenotype_base_sskel::
  mass_present ()
  {
    return this->phenotype_base_impl_ ? this->phenotype_base_impl_->mass_present () : false;
  }

  bool phenotype_base_sskel::
  mechanics_present ()
  {
    return this->phenotype_base_impl_ ? this->phenotype_base_impl_->mechanics_present () : false;
  }

  bool phenotype_base_sskel::
  motility_present ()
  {
    return this->phenotype_base_impl_ ? this->phenotype_base_impl_->motility_present () : false;
  }

  bool phenotype_base_sskel::
  PKPD_present ()
  {
    return this->phenotype_base_impl_ ? this->phenotype_base_impl_->PKPD_present () : false;
  }

  bool phenotype_base_sskel::
  timescale_present ()
  {
    return this->phenotype_base_impl_ ? this->phenotype_base_impl_->timescale_present () : false;
  }

  bool phenotype_base_sskel::
  transport_processes_present ()
  {
    return this->phenotype_base_impl_ ? this->phenotype_base_impl_->transport_processes_present () : false;
  }

  bool phenotype_base_sskel::
  custom_present ()
  {
    return this->phenotype_base_impl_ ? this->phenotype_base_impl_->custom_present () : false;
  }

  void phenotype_base_sskel::
  _reset ()
  {
    if (this->resetting_)
      return;

    typedef ::xsde::cxx::serializer::validating::complex_content base;
    base::_reset ();

    if (this->type_serializer_)
      this->type_serializer_->_reset ();

    this->resetting_ = true;

    if (this->adhesion_serializer_)
      this->adhesion_serializer_->_reset ();

    if (this->adhesion_serializer_map_)
      this->adhesion_serializer_map_->reset ();

    if (this->geometrical_properties_serializer_)
      this->geometrical_properties_serializer_->_reset ();

    if (this->geometrical_properties_serializer_map_)
      this->geometrical_properties_serializer_map_->reset ();

    if (this->mass_serializer_)
      this->mass_serializer_->_reset ();

    if (this->mass_serializer_map_)
      this->mass_serializer_map_->reset ();

    if (this->mechanics_serializer_)
      this->mechanics_serializer_->_reset ();

    if (this->mechanics_serializer_map_)
      this->mechanics_serializer_map_->reset ();

    if (this->motility_serializer_)
      this->motility_serializer_->_reset ();

    if (this->motility_serializer_map_)
      this->motility_serializer_map_->reset ();

    if (this->PKPD_serializer_)
      this->PKPD_serializer_->_reset ();

    if (this->PKPD_serializer_map_)
      this->PKPD_serializer_map_->reset ();

    if (this->timescale_serializer_)
      this->timescale_serializer_->_reset ();

    if (this->timescale_serializer_map_)
      this->timescale_serializer_map_->reset ();

    if (this->transport_processes_serializer_)
      this->transport_processes_serializer_->_reset ();

    if (this->transport_processes_serializer_map_)
      this->transport_processes_serializer_map_->reset ();

    if (this->custom_serializer_)
      this->custom_serializer_->_reset ();

    if (this->custom_serializer_map_)
      this->custom_serializer_map_->reset ();

    this->resetting_ = false;
  }

  const char* phenotype_base_sskel::
  _static_type ()
  {
    return "phenotype_base phenotype_base";
  }

  const char* phenotype_base_sskel::
  _dynamic_type () const
  {
    return _static_type ();
  }

  // expected_timescale_sskel
  //

  bool expected_timescale_sskel::
  cell_cycle_ID_present ()
  {
    return this->expected_timescale_impl_ ? this->expected_timescale_impl_->cell_cycle_ID_present () : false;
  }

  bool expected_timescale_sskel::
  cell_cycle_phase_ID_present ()
  {
    return this->expected_timescale_impl_ ? this->expected_timescale_impl_->cell_cycle_phase_ID_present () : false;
  }

  bool expected_timescale_sskel::
  cell_death_ID_present ()
  {
    return this->expected_timescale_impl_ ? this->expected_timescale_impl_->cell_death_ID_present () : false;
  }

  void expected_timescale_sskel::
  _reset ()
  {
    typedef ::common::units_decimal_nonnegative_sskel base;
    base::_reset ();

    if (this->cell_cycle_ID_serializer_)
      this->cell_cycle_ID_serializer_->_reset ();

    if (this->cell_cycle_phase_ID_serializer_)
      this->cell_cycle_phase_ID_serializer_->_reset ();

    if (this->cell_death_ID_serializer_)
      this->cell_death_ID_serializer_->_reset ();
  }

  const char* expected_timescale_sskel::
  _static_type ()
  {
    return "expected_timescale phenotype_base";
  }

  const char* expected_timescale_sskel::
  _dynamic_type () const
  {
    return _static_type ();
  }

  static
  const ::xsde::cxx::serializer::validating::inheritance_map_entry
  _xsde_expected_timescale_sskel_inheritance_map_entry_ (
    expected_timescale_sskel::_static_type (),
    ::common::units_decimal_nonnegative_sskel::_static_type ());

  void expected_timescale_sskel::
  pre (const ::common::units_decimal_nonnegative& x)
  {
    assert (this->units_decimal_nonnegative_impl_);
    this->units_decimal_nonnegative_impl_->pre (x);
  }

  // cell_parts_sskel
  //

  bool cell_parts_sskel::
  ID_present ()
  {
    return this->cell_parts_impl_ ? this->cell_parts_impl_->ID_present () : false;
  }

  bool cell_parts_sskel::
  phenotype_next ()
  {
    return this->cell_parts_impl_ ? this->cell_parts_impl_->phenotype_next () : false;
  }

  bool cell_parts_sskel::
  cell_part_next ()
  {
    return this->cell_parts_impl_ ? this->cell_parts_impl_->cell_part_next () : false;
  }

  bool cell_parts_sskel::
  custom_present ()
  {
    return this->cell_parts_impl_ ? this->cell_parts_impl_->custom_present () : false;
  }

  void cell_parts_sskel::
  _reset ()
  {
    if (this->resetting_)
      return;

    typedef ::xsde::cxx::serializer::validating::complex_content base;
    base::_reset ();

    if (this->name_serializer_)
      this->name_serializer_->_reset ();

    if (this->ID_serializer_)
      this->ID_serializer_->_reset ();

    this->resetting_ = true;

    if (this->phenotype_serializer_)
      this->phenotype_serializer_->_reset ();

    if (this->phenotype_serializer_map_)
      this->phenotype_serializer_map_->reset ();

    if (this->cell_part_serializer_)
      this->cell_part_serializer_->_reset ();

    if (this->cell_part_serializer_map_)
      this->cell_part_serializer_map_->reset ();

    if (this->custom_serializer_)
      this->custom_serializer_->_reset ();

    if (this->custom_serializer_map_)
      this->custom_serializer_map_->reset ();

    this->resetting_ = false;
  }

  const char* cell_parts_sskel::
  _static_type ()
  {
    return "cell_parts phenotype_base";
  }

  const char* cell_parts_sskel::
  _dynamic_type () const
  {
    return _static_type ();
  }
}

namespace phenotype_base
{
  // Element validation and serialization for phenotype_base_sskel.
  //
  void phenotype_base_sskel::
  _serialize_content ()
  {
    ::xsde::cxx::serializer::context& ctx = this->_context ();

    // adhesion
    //
    if (this->adhesion_present ())
    {
      ctx.type_id (0);
      const ::phenotype_common::adhesion& r = this->adhesion ();

      const void* t = ctx.type_id ();
      const char* dt = 0;
      ::phenotype_common::adhesion_sskel* s = 0;

      if (t == 0 && this->adhesion_serializer_ != 0)
        s = this->adhesion_serializer_;
      else if (this->adhesion_serializer_map_ != 0)
      {
        ::xml_schema::serializer_base* b = this->adhesion_serializer_map_->find (t);

        if (b != 0)
        {
          dt = b->_dynamic_type ();
          const char* st = ::phenotype_common::adhesion_sskel::_static_type ();

          if (strcmp (dt, st) == 0)
            dt = 0;

          if (dt != 0 && !::xsde::cxx::serializer::validating::inheritance_map_instance ().check (dt, st))
          {
            ctx.schema_error (::xsde::cxx::schema_error::not_derived);
            return;
          }

          s = static_cast< ::phenotype_common::adhesion_sskel* > (b);
        }
      }

      if (s)
      {
        s->pre (r);
        this->_start_element ("adhesion");
        if (dt != 0)
          this->_set_type (dt);

        s->_pre_impl (ctx);

        if (ctx.error_type ())
          return;

        s->_serialize_attributes ();

        if (ctx.error_type ())
          return;

        s->_serialize_content ();

        if (ctx.error_type ())
          return;

        s->_post_impl ();

        if (ctx.error_type ())
          return;

        this->_end_element ();
        s->post ();
      }
    }

    // geometrical_properties
    //
    if (this->geometrical_properties_present ())
    {
      ctx.type_id (0);
      const ::phenotype_common::geometrical_properties& r = this->geometrical_properties ();

      const void* t = ctx.type_id ();
      const char* dt = 0;
      ::phenotype_common::geometrical_properties_sskel* s = 0;

      if (t == 0 && this->geometrical_properties_serializer_ != 0)
        s = this->geometrical_properties_serializer_;
      else if (this->geometrical_properties_serializer_map_ != 0)
      {
        ::xml_schema::serializer_base* b = this->geometrical_properties_serializer_map_->find (t);

        if (b != 0)
        {
          dt = b->_dynamic_type ();
          const char* st = ::phenotype_common::geometrical_properties_sskel::_static_type ();

          if (strcmp (dt, st) == 0)
            dt = 0;

          if (dt != 0 && !::xsde::cxx::serializer::validating::inheritance_map_instance ().check (dt, st))
          {
            ctx.schema_error (::xsde::cxx::schema_error::not_derived);
            return;
          }

          s = static_cast< ::phenotype_common::geometrical_properties_sskel* > (b);
        }
      }

      if (s)
      {
        s->pre (r);
        this->_start_element ("geometrical_properties");
        if (dt != 0)
          this->_set_type (dt);

        s->_pre_impl (ctx);

        if (ctx.error_type ())
          return;

        s->_serialize_attributes ();

        if (ctx.error_type ())
          return;

        s->_serialize_content ();

        if (ctx.error_type ())
          return;

        s->_post_impl ();

        if (ctx.error_type ())
          return;

        this->_end_element ();
        s->post ();
      }
    }

    // mass
    //
    if (this->mass_present ())
    {
      ctx.type_id (0);
      const ::phenotype_common::mass& r = this->mass ();

      const void* t = ctx.type_id ();
      const char* dt = 0;
      ::phenotype_common::mass_sskel* s = 0;

      if (t == 0 && this->mass_serializer_ != 0)
        s = this->mass_serializer_;
      else if (this->mass_serializer_map_ != 0)
      {
        ::xml_schema::serializer_base* b = this->mass_serializer_map_->find (t);

        if (b != 0)
        {
          dt = b->_dynamic_type ();
          const char* st = ::phenotype_common::mass_sskel::_static_type ();

          if (strcmp (dt, st) == 0)
            dt = 0;

          if (dt != 0 && !::xsde::cxx::serializer::validating::inheritance_map_instance ().check (dt, st))
          {
            ctx.schema_error (::xsde::cxx::schema_error::not_derived);
            return;
          }

          s = static_cast< ::phenotype_common::mass_sskel* > (b);
        }
      }

      if (s)
      {
        s->pre (r);
        this->_start_element ("mass");
        if (dt != 0)
          this->_set_type (dt);

        s->_pre_impl (ctx);

        if (ctx.error_type ())
          return;

        s->_serialize_attributes ();

        if (ctx.error_type ())
          return;

        s->_serialize_content ();

        if (ctx.error_type ())
          return;

        s->_post_impl ();

        if (ctx.error_type ())
          return;

        this->_end_element ();
        s->post ();
      }
    }

    // mechanics
    //
    if (this->mechanics_present ())
    {
      ctx.type_id (0);
      const ::phenotype_common::mechanics& r = this->mechanics ();

      const void* t = ctx.type_id ();
      const char* dt = 0;
      ::phenotype_common::mechanics_sskel* s = 0;

      if (t == 0 && this->mechanics_serializer_ != 0)
        s = this->mechanics_serializer_;
      else if (this->mechanics_serializer_map_ != 0)
      {
        ::xml_schema::serializer_base* b = this->mechanics_serializer_map_->find (t);

        if (b != 0)
        {
          dt = b->_dynamic_type ();
          const char* st = ::phenotype_common::mechanics_sskel::_static_type ();

          if (strcmp (dt, st) == 0)
            dt = 0;

          if (dt != 0 && !::xsde::cxx::serializer::validating::inheritance_map_instance ().check (dt, st))
          {
            ctx.schema_error (::xsde::cxx::schema_error::not_derived);
            return;
          }

          s = static_cast< ::phenotype_common::mechanics_sskel* > (b);
        }
      }

      if (s)
      {
        s->pre (r);
        this->_start_element ("mechanics");
        if (dt != 0)
          this->_set_type (dt);

        s->_pre_impl (ctx);

        if (ctx.error_type ())
          return;

        s->_serialize_attributes ();

        if (ctx.error_type ())
          return;

        s->_serialize_content ();

        if (ctx.error_type ())
          return;

        s->_post_impl ();

        if (ctx.error_type ())
          return;

        this->_end_element ();
        s->post ();
      }
    }

    // motility
    //
    if (this->motility_present ())
    {
      ctx.type_id (0);
      const ::phenotype_common::motility& r = this->motility ();

      const void* t = ctx.type_id ();
      const char* dt = 0;
      ::phenotype_common::motility_sskel* s = 0;

      if (t == 0 && this->motility_serializer_ != 0)
        s = this->motility_serializer_;
      else if (this->motility_serializer_map_ != 0)
      {
        ::xml_schema::serializer_base* b = this->motility_serializer_map_->find (t);

        if (b != 0)
        {
          dt = b->_dynamic_type ();
          const char* st = ::phenotype_common::motility_sskel::_static_type ();

          if (strcmp (dt, st) == 0)
            dt = 0;

          if (dt != 0 && !::xsde::cxx::serializer::validating::inheritance_map_instance ().check (dt, st))
          {
            ctx.schema_error (::xsde::cxx::schema_error::not_derived);
            return;
          }

          s = static_cast< ::phenotype_common::motility_sskel* > (b);
        }
      }

      if (s)
      {
        s->pre (r);
        this->_start_element ("motility");
        if (dt != 0)
          this->_set_type (dt);

        s->_pre_impl (ctx);

        if (ctx.error_type ())
          return;

        s->_serialize_attributes ();

        if (ctx.error_type ())
          return;

        s->_serialize_content ();

        if (ctx.error_type ())
          return;

        s->_post_impl ();

        if (ctx.error_type ())
          return;

        this->_end_element ();
        s->post ();
      }
    }

    // PKPD
    //
    if (this->PKPD_present ())
    {
      ctx.type_id (0);
      const ::pkpd::PKPD& r = this->PKPD ();

      const void* t = ctx.type_id ();
      const char* dt = 0;
      ::pkpd::PKPD_sskel* s = 0;

      if (t == 0 && this->PKPD_serializer_ != 0)
        s = this->PKPD_serializer_;
      else if (this->PKPD_serializer_map_ != 0)
      {
        ::xml_schema::serializer_base* b = this->PKPD_serializer_map_->find (t);

        if (b != 0)
        {
          dt = b->_dynamic_type ();
          const char* st = ::pkpd::PKPD_sskel::_static_type ();

          if (strcmp (dt, st) == 0)
            dt = 0;

          if (dt != 0 && !::xsde::cxx::serializer::validating::inheritance_map_instance ().check (dt, st))
          {
            ctx.schema_error (::xsde::cxx::schema_error::not_derived);
            return;
          }

          s = static_cast< ::pkpd::PKPD_sskel* > (b);
        }
      }

      if (s)
      {
        s->pre (r);
        this->_start_element ("PKPD");
        if (dt != 0)
          this->_set_type (dt);

        s->_pre_impl (ctx);

        if (ctx.error_type ())
          return;

        s->_serialize_attributes ();

        if (ctx.error_type ())
          return;

        s->_serialize_content ();

        if (ctx.error_type ())
          return;

        s->_post_impl ();

        if (ctx.error_type ())
          return;

        this->_end_element ();
        s->post ();
      }
    }

    // timescale
    //
    if (this->timescale_present ())
    {
      ctx.type_id (0);
      const ::phenotype_base::expected_timescale& r = this->timescale ();

      const void* t = ctx.type_id ();
      const char* dt = 0;
      ::phenotype_base::expected_timescale_sskel* s = 0;

      if (t == 0 && this->timescale_serializer_ != 0)
        s = this->timescale_serializer_;
      else if (this->timescale_serializer_map_ != 0)
      {
        ::xml_schema::serializer_base* b = this->timescale_serializer_map_->find (t);

        if (b != 0)
        {
          dt = b->_dynamic_type ();
          const char* st = ::phenotype_base::expected_timescale_sskel::_static_type ();

          if (strcmp (dt, st) == 0)
            dt = 0;

          if (dt != 0 && !::xsde::cxx::serializer::validating::inheritance_map_instance ().check (dt, st))
          {
            ctx.schema_error (::xsde::cxx::schema_error::not_derived);
            return;
          }

          s = static_cast< ::phenotype_base::expected_timescale_sskel* > (b);
        }
      }

      if (s)
      {
        s->pre (r);
        this->_start_element ("timescale");
        if (dt != 0)
          this->_set_type (dt);

        s->_pre_impl (ctx);

        if (ctx.error_type ())
          return;

        s->_serialize_attributes ();

        if (ctx.error_type ())
          return;

        s->_serialize_content ();

        if (ctx.error_type ())
          return;

        s->_post_impl ();

        if (ctx.error_type ())
          return;

        this->_end_element ();
        s->post ();
      }
    }

    // transport_processes
    //
    if (this->transport_processes_present ())
    {
      ctx.type_id (0);
      const ::phenotype_common::transport_processes& r = this->transport_processes ();

      const void* t = ctx.type_id ();
      const char* dt = 0;
      ::phenotype_common::transport_processes_sskel* s = 0;

      if (t == 0 && this->transport_processes_serializer_ != 0)
        s = this->transport_processes_serializer_;
      else if (this->transport_processes_serializer_map_ != 0)
      {
        ::xml_schema::serializer_base* b = this->transport_processes_serializer_map_->find (t);

        if (b != 0)
        {
          dt = b->_dynamic_type ();
          const char* st = ::phenotype_common::transport_processes_sskel::_static_type ();

          if (strcmp (dt, st) == 0)
            dt = 0;

          if (dt != 0 && !::xsde::cxx::serializer::validating::inheritance_map_instance ().check (dt, st))
          {
            ctx.schema_error (::xsde::cxx::schema_error::not_derived);
            return;
          }

          s = static_cast< ::phenotype_common::transport_processes_sskel* > (b);
        }
      }

      if (s)
      {
        s->pre (r);
        this->_start_element ("transport_processes");
        if (dt != 0)
          this->_set_type (dt);

        s->_pre_impl (ctx);

        if (ctx.error_type ())
          return;

        s->_serialize_attributes ();

        if (ctx.error_type ())
          return;

        s->_serialize_content ();

        if (ctx.error_type ())
          return;

        s->_post_impl ();

        if (ctx.error_type ())
          return;

        this->_end_element ();
        s->post ();
      }
    }

    // custom
    //
    if (this->custom_present ())
    {
      ctx.type_id (0);
      const ::common::custom& r = this->custom ();

      const void* t = ctx.type_id ();
      const char* dt = 0;
      ::common::custom_sskel* s = 0;

      if (t == 0 && this->custom_serializer_ != 0)
        s = this->custom_serializer_;
      else if (this->custom_serializer_map_ != 0)
      {
        ::xml_schema::serializer_base* b = this->custom_serializer_map_->find (t);

        if (b != 0)
        {
          dt = b->_dynamic_type ();
          const char* st = ::common::custom_sskel::_static_type ();

          if (strcmp (dt, st) == 0)
            dt = 0;

          if (dt != 0 && !::xsde::cxx::serializer::validating::inheritance_map_instance ().check (dt, st))
          {
            ctx.schema_error (::xsde::cxx::schema_error::not_derived);
            return;
          }

          s = static_cast< ::common::custom_sskel* > (b);
        }
      }

      if (s)
      {
        s->pre (r);
        this->_start_element ("custom");
        if (dt != 0)
          this->_set_type (dt);

        s->_pre_impl (ctx);

        if (ctx.error_type ())
          return;

        s->_serialize_attributes ();

        if (ctx.error_type ())
          return;

        s->_serialize_content ();

        if (ctx.error_type ())
          return;

        s->_post_impl ();

        if (ctx.error_type ())
          return;

        this->_end_element ();
        s->post ();
      }
    }
  }

  // Element validation and serialization for cell_parts_sskel.
  //
  void cell_parts_sskel::
  _serialize_content ()
  {
    ::xsde::cxx::serializer::context& ctx = this->_context ();

    // phenotype
    //
    {
      size_t i = 0;
      for (; i < 3UL && this->phenotype_next (); ++i)
      {
        ctx.type_id (0);
        const ::phenotype_base::phenotype_base& r = this->phenotype ();

        const void* t = ctx.type_id ();
        const char* dt = 0;
        ::phenotype_base::phenotype_base_sskel* s = 0;

        if (t == 0 && this->phenotype_serializer_ != 0)
          s = this->phenotype_serializer_;
        else if (this->phenotype_serializer_map_ != 0)
        {
          ::xml_schema::serializer_base* b = this->phenotype_serializer_map_->find (t);

          if (b != 0)
          {
            dt = b->_dynamic_type ();
            const char* st = ::phenotype_base::phenotype_base_sskel::_static_type ();

            if (strcmp (dt, st) == 0)
              dt = 0;

            if (dt != 0 && !::xsde::cxx::serializer::validating::inheritance_map_instance ().check (dt, st))
            {
              ctx.schema_error (::xsde::cxx::schema_error::not_derived);
              return;
            }

            s = static_cast< ::phenotype_base::phenotype_base_sskel* > (b);
          }
        }

        if (s)
        {
          s->pre (r);
          this->_start_element ("phenotype");
          if (dt != 0)
            this->_set_type (dt);

          s->_pre_impl (ctx);

          if (ctx.error_type ())
            return;

          s->_serialize_attributes ();

          if (ctx.error_type ())
            return;

          s->_serialize_content ();

          if (ctx.error_type ())
            return;

          s->_post_impl ();

          if (ctx.error_type ())
            return;

          this->_end_element ();
          s->post ();
        }
      }
    }

    // cell_part
    //
    while (this->cell_part_next ())
    {
      ctx.type_id (0);
      const ::phenotype_base::cell_parts& r = this->cell_part ();

      const void* t = ctx.type_id ();
      const char* dt = 0;
      ::phenotype_base::cell_parts_sskel* s = 0;

      if (t == 0 && this->cell_part_serializer_ != 0)
        s = this->cell_part_serializer_;
      else if (this->cell_part_serializer_map_ != 0)
      {
        ::xml_schema::serializer_base* b = this->cell_part_serializer_map_->find (t);

        if (b != 0)
        {
          dt = b->_dynamic_type ();
          const char* st = ::phenotype_base::cell_parts_sskel::_static_type ();

          if (strcmp (dt, st) == 0)
            dt = 0;

          if (dt != 0 && !::xsde::cxx::serializer::validating::inheritance_map_instance ().check (dt, st))
          {
            ctx.schema_error (::xsde::cxx::schema_error::not_derived);
            return;
          }

          s = static_cast< ::phenotype_base::cell_parts_sskel* > (b);
        }
      }

      if (s)
      {
        s->pre (r);
        this->_start_element ("cell_part");
        if (dt != 0)
          this->_set_type (dt);

        s->_pre_impl (ctx);

        if (ctx.error_type ())
          return;

        s->_serialize_attributes ();

        if (ctx.error_type ())
          return;

        s->_serialize_content ();

        if (ctx.error_type ())
          return;

        s->_post_impl ();

        if (ctx.error_type ())
          return;

        this->_end_element ();
        s->post ();
      }
    }

    // custom
    //
    if (this->custom_present ())
    {
      ctx.type_id (0);
      const ::common::custom& r = this->custom ();

      const void* t = ctx.type_id ();
      const char* dt = 0;
      ::common::custom_sskel* s = 0;

      if (t == 0 && this->custom_serializer_ != 0)
        s = this->custom_serializer_;
      else if (this->custom_serializer_map_ != 0)
      {
        ::xml_schema::serializer_base* b = this->custom_serializer_map_->find (t);

        if (b != 0)
        {
          dt = b->_dynamic_type ();
          const char* st = ::common::custom_sskel::_static_type ();

          if (strcmp (dt, st) == 0)
            dt = 0;

          if (dt != 0 && !::xsde::cxx::serializer::validating::inheritance_map_instance ().check (dt, st))
          {
            ctx.schema_error (::xsde::cxx::schema_error::not_derived);
            return;
          }

          s = static_cast< ::common::custom_sskel* > (b);
        }
      }

      if (s)
      {
        s->pre (r);
        this->_start_element ("custom");
        if (dt != 0)
          this->_set_type (dt);

        s->_pre_impl (ctx);

        if (ctx.error_type ())
          return;

        s->_serialize_attributes ();

        if (ctx.error_type ())
          return;

        s->_serialize_content ();

        if (ctx.error_type ())
          return;

        s->_post_impl ();

        if (ctx.error_type ())
          return;

        this->_end_element ();
        s->post ();
      }
    }
  }
}

namespace phenotype_base
{
  // Attribute validation and serialization for phenotype_base_sskel.
  //
  void phenotype_base_sskel::
  _serialize_attributes ()
  {
    ::xsde::cxx::serializer::context& ctx = this->_context ();

    // type
    //
    if (this->type_present ())
    {
      const ::phenotype_base::phenotype_type& r = this->type ();

      if (this->type_serializer_)
      {
        this->type_serializer_->pre (r);
        this->_start_attribute ("type");
        this->type_serializer_->_pre_impl (ctx);

        if (ctx.error_type ())
          return;

        this->type_serializer_->_serialize_content ();

        if (ctx.error_type ())
          return;

        this->type_serializer_->_post_impl ();

        if (ctx.error_type ())
          return;

        this->_end_attribute ();
        this->type_serializer_->post ();
      }
    }
  }

  // Attribute validation and serialization for expected_timescale_sskel.
  //
  void expected_timescale_sskel::
  _serialize_attributes ()
  {
    ::xsde::cxx::serializer::context& ctx = this->_context ();

    typedef ::common::units_decimal_nonnegative_sskel base;
    base::_serialize_attributes ();

    if (ctx.error_type ())
      return;

    // cell_cycle_ID
    //
    if (this->cell_cycle_ID_present ())
    {
      unsigned int r = this->cell_cycle_ID ();

      if (this->cell_cycle_ID_serializer_)
      {
        this->cell_cycle_ID_serializer_->pre (r);
        this->_start_attribute ("cell_cycle_ID");
        this->cell_cycle_ID_serializer_->_pre_impl (ctx);

        if (ctx.error_type ())
          return;

        this->cell_cycle_ID_serializer_->_serialize_content ();

        if (ctx.error_type ())
          return;

        this->cell_cycle_ID_serializer_->_post_impl ();

        if (ctx.error_type ())
          return;

        this->_end_attribute ();
        this->cell_cycle_ID_serializer_->post ();
      }
    }

    // cell_cycle_phase_ID
    //
    if (this->cell_cycle_phase_ID_present ())
    {
      unsigned int r = this->cell_cycle_phase_ID ();

      if (this->cell_cycle_phase_ID_serializer_)
      {
        this->cell_cycle_phase_ID_serializer_->pre (r);
        this->_start_attribute ("cell_cycle_phase_ID");
        this->cell_cycle_phase_ID_serializer_->_pre_impl (ctx);

        if (ctx.error_type ())
          return;

        this->cell_cycle_phase_ID_serializer_->_serialize_content ();

        if (ctx.error_type ())
          return;

        this->cell_cycle_phase_ID_serializer_->_post_impl ();

        if (ctx.error_type ())
          return;

        this->_end_attribute ();
        this->cell_cycle_phase_ID_serializer_->post ();
      }
    }

    // cell_death_ID
    //
    if (this->cell_death_ID_present ())
    {
      unsigned int r = this->cell_death_ID ();

      if (this->cell_death_ID_serializer_)
      {
        this->cell_death_ID_serializer_->pre (r);
        this->_start_attribute ("cell_death_ID");
        this->cell_death_ID_serializer_->_pre_impl (ctx);

        if (ctx.error_type ())
          return;

        this->cell_death_ID_serializer_->_serialize_content ();

        if (ctx.error_type ())
          return;

        this->cell_death_ID_serializer_->_post_impl ();

        if (ctx.error_type ())
          return;

        this->_end_attribute ();
        this->cell_death_ID_serializer_->post ();
      }
    }
  }

  // Attribute validation and serialization for cell_parts_sskel.
  //
  void cell_parts_sskel::
  _serialize_attributes ()
  {
    ::xsde::cxx::serializer::context& ctx = this->_context ();

    // name
    //
    {
      const ::std::string& r = this->name ();

      if (this->name_serializer_)
      {
        this->name_serializer_->pre (r);
        this->_start_attribute ("name");
        this->name_serializer_->_pre_impl (ctx);

        if (ctx.error_type ())
          return;

        this->name_serializer_->_serialize_content ();

        if (ctx.error_type ())
          return;

        this->name_serializer_->_post_impl ();

        if (ctx.error_type ())
          return;

        this->_end_attribute ();
        this->name_serializer_->post ();
      }
      else
      {
        this->_schema_error (::xsde::cxx::schema_error::expected_attribute);
        return;
      }
    }

    // ID
    //
    if (this->ID_present ())
    {
      unsigned int r = this->ID ();

      if (this->ID_serializer_)
      {
        this->ID_serializer_->pre (r);
        this->_start_attribute ("ID");
        this->ID_serializer_->_pre_impl (ctx);

        if (ctx.error_type ())
          return;

        this->ID_serializer_->_serialize_content ();

        if (ctx.error_type ())
          return;

        this->ID_serializer_->_post_impl ();

        if (ctx.error_type ())
          return;

        this->_end_attribute ();
        this->ID_serializer_->post ();
      }
    }
  }
}

// Begin epilogue.
//
//
// End epilogue.

