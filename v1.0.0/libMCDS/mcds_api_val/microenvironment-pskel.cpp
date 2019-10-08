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

#include "microenvironment-pskel.hpp"

#include <assert.h>

#include <string.h>
#include <xsde/cxx/parser/substitution-map.hxx>
#include <xsde/cxx/parser/validating/inheritance-map.hxx>

static
const ::xsde::cxx::parser::substitution_map_init
_xsde_substitution_map_init_;

static
const ::xsde::cxx::parser::validating::inheritance_map_init
_xsde_inheritance_map_init_;

namespace microenvironment
{
  // domain_pskel
  //

  void domain_pskel::
  name (const ::std::string& x)
  {
    if (this->domain_impl_)
      this->domain_impl_->name (x);
  }

  void domain_pskel::
  variables (::variables::list_of_variables* x)
  {
    if (this->domain_impl_)
      this->domain_impl_->variables (x);
  }

  void domain_pskel::
  experimental_condition (::variables::experimental_conditions* x)
  {
    if (this->domain_impl_)
      this->domain_impl_->experimental_condition (x);
  }

  void domain_pskel::
  mesh (::mesh::mesh* x)
  {
    if (this->domain_impl_)
      this->domain_impl_->mesh (x);
  }

  void domain_pskel::
  data (::variables::data* x)
  {
    if (this->domain_impl_)
      this->domain_impl_->data (x);
  }

  void domain_pskel::
  custom (::common::custom* x)
  {
    if (this->domain_impl_)
      this->domain_impl_->custom (x);
  }

  void domain_pskel::
  _reset ()
  {
    if (this->resetting_)
      return;

    typedef ::xsde::cxx::parser::validating::complex_content base;
    base::_reset ();

    this->v_state_stack_.clear ();

    this->v_all_count_.clear ();

    if (this->name_parser_)
      this->name_parser_->_reset ();

    this->resetting_ = true;

    if (this->variables_parser_)
      this->variables_parser_->_reset ();

    if (this->variables_parser_map_)
      this->variables_parser_map_->reset ();

    if (this->experimental_condition_parser_)
      this->experimental_condition_parser_->_reset ();

    if (this->experimental_condition_parser_map_)
      this->experimental_condition_parser_map_->reset ();

    if (this->mesh_parser_)
      this->mesh_parser_->_reset ();

    if (this->mesh_parser_map_)
      this->mesh_parser_map_->reset ();

    if (this->data_parser_)
      this->data_parser_->_reset ();

    if (this->data_parser_map_)
      this->data_parser_map_->reset ();

    if (this->custom_parser_)
      this->custom_parser_->_reset ();

    if (this->custom_parser_map_)
      this->custom_parser_map_->reset ();

    this->resetting_ = false;
  }

  const char* domain_pskel::
  _static_type ()
  {
    return "domain microenvironment";
  }

  const char* domain_pskel::
  _dynamic_type () const
  {
    return _static_type ();
  }

  // microenvironment_pskel
  //

  void microenvironment_pskel::
  domain (::microenvironment::domain* x)
  {
    if (this->microenvironment_impl_)
      this->microenvironment_impl_->domain (x);
  }

  void microenvironment_pskel::
  vascular_network (::vascular::vascular_network* x)
  {
    if (this->microenvironment_impl_)
      this->microenvironment_impl_->vascular_network (x);
  }

  void microenvironment_pskel::
  basement_membrane (::basement::basement_membrane* x)
  {
    if (this->microenvironment_impl_)
      this->microenvironment_impl_->basement_membrane (x);
  }

  void microenvironment_pskel::
  custom (::common::custom* x)
  {
    if (this->microenvironment_impl_)
      this->microenvironment_impl_->custom (x);
  }

  void microenvironment_pskel::
  _reset ()
  {
    if (this->resetting_)
      return;

    typedef ::xsde::cxx::parser::validating::complex_content base;
    base::_reset ();

    this->v_state_stack_.clear ();

    this->resetting_ = true;

    if (this->domain_parser_)
      this->domain_parser_->_reset ();

    if (this->domain_parser_map_)
      this->domain_parser_map_->reset ();

    if (this->vascular_network_parser_)
      this->vascular_network_parser_->_reset ();

    if (this->vascular_network_parser_map_)
      this->vascular_network_parser_map_->reset ();

    if (this->basement_membrane_parser_)
      this->basement_membrane_parser_->_reset ();

    if (this->basement_membrane_parser_map_)
      this->basement_membrane_parser_map_->reset ();

    if (this->custom_parser_)
      this->custom_parser_->_reset ();

    if (this->custom_parser_map_)
      this->custom_parser_map_->reset ();

    this->resetting_ = false;
  }

  const char* microenvironment_pskel::
  _static_type ()
  {
    return "microenvironment microenvironment";
  }

  const char* microenvironment_pskel::
  _dynamic_type () const
  {
    return _static_type ();
  }
}

#include <assert.h>

namespace microenvironment
{
  // Element validation and dispatch functions for domain_pskel.
  //
  bool domain_pskel::
  _start_element_impl (const ::xsde::cxx::ro_string& ns,
                       const ::xsde::cxx::ro_string& n,
                       const char* t)
  {
    XSDE_UNUSED (t);

    ::xsde::cxx::parser::context& ctx = this->_context ();

    v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
    v_state_descr_* vd = vs.data + (vs.size - 1);

    all_0 (vd->state, v_all_count_.top (), ns, n, t, true);

    if (vd->state != ~0UL || ctx.error_type ())
      vd->count++;
    else
      return false;

    return true;
  }

  bool domain_pskel::
  _end_element_impl (const ::xsde::cxx::ro_string& ns,
                     const ::xsde::cxx::ro_string& n)
  {
    v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
    v_state_descr_& vd = vs.data[vs.size - 1];

    all_0 (vd.state, v_all_count_.top (), ns, n, 0, false);

    return true;
  }

  void domain_pskel::
  _pre_e_validate ()
  {
    this->v_state_stack_.push ();
    static_cast< v_state_* > (this->v_state_stack_.top ())->size = 0;

    v_all_count_.push ();

    v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
    v_state_descr_& vd = vs.data[vs.size++];

    vd.func = 0;
    vd.state = 0;
    vd.count = 0;
  }

  void domain_pskel::
  _post_e_validate ()
  {
    v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
    v_state_descr_& vd = vs.data[vs.size - 1];

    if (vd.count != 0)
    {
      ::xsde::cxx::ro_string empty;
      all_0 (vd.state, v_all_count_.top (), empty, empty, 0, true);
    }


    vs.size--;
    v_all_count_.pop ();

    this->v_state_stack_.pop ();
  }

  void domain_pskel::
  all_0 (unsigned long& state,
         unsigned char* count,
         const ::xsde::cxx::ro_string& ns,
         const ::xsde::cxx::ro_string& n,
         const char* t,
         bool start)
  {
    XSDE_UNUSED (t);

    ::xsde::cxx::parser::context& ctx = this->_context ();

    if (n == "variables" && ns.empty ())
    {
      if (count[0UL] == 0)
      {
        if (start)
        {
          ::variables::list_of_variables_pskel* p = 0;

          if (t == 0 && this->variables_parser_ != 0)
            p = this->variables_parser_;
          else
          {
            const char* ts = ::variables::list_of_variables_pskel::_static_type ();

            if (t == 0)
              t = ts;

            if (this->variables_parser_ != 0 && strcmp (t, ts) == 0)
              p = this->variables_parser_;
            else
            {
              if (t != ts &&
                  !::xsde::cxx::parser::validating::inheritance_map_instance ().check (t, ts))
              {
                ctx.schema_error (::xsde::cxx::schema_error::not_derived);
                return;
              }

              if (this->variables_parser_map_ != 0)
                p = static_cast< ::variables::list_of_variables_pskel* > (
                  this->variables_parser_map_->find (t));
            }
          }

          if (p)
          {
            p->pre ();
            ctx.nested_parser (p);
          }
        }
        else
        {
          ::variables::list_of_variables_pskel* p =
          static_cast< ::variables::list_of_variables_pskel* > (ctx.nested_parser ());

          if (p != 0)
          {
            ::variables::list_of_variables* tmp = p->post_list_of_variables ();
            this->variables (tmp);
          }

          count[0UL] = 1;
        }
      }
      else
      {
        assert (start);
        state = ~0UL;
      }
    }
    else if (n == "experimental_condition" && ns.empty ())
    {
      if (count[1UL] == 0)
      {
        if (start)
        {
          ::variables::experimental_conditions_pskel* p = 0;

          if (t == 0 && this->experimental_condition_parser_ != 0)
            p = this->experimental_condition_parser_;
          else
          {
            const char* ts = ::variables::experimental_conditions_pskel::_static_type ();

            if (t == 0)
              t = ts;

            if (this->experimental_condition_parser_ != 0 && strcmp (t, ts) == 0)
              p = this->experimental_condition_parser_;
            else
            {
              if (t != ts &&
                  !::xsde::cxx::parser::validating::inheritance_map_instance ().check (t, ts))
              {
                ctx.schema_error (::xsde::cxx::schema_error::not_derived);
                return;
              }

              if (this->experimental_condition_parser_map_ != 0)
                p = static_cast< ::variables::experimental_conditions_pskel* > (
                  this->experimental_condition_parser_map_->find (t));
            }
          }

          if (p)
          {
            p->pre ();
            ctx.nested_parser (p);
          }
        }
        else
        {
          ::variables::experimental_conditions_pskel* p =
          static_cast< ::variables::experimental_conditions_pskel* > (ctx.nested_parser ());

          if (p != 0)
          {
            ::variables::experimental_conditions* tmp = p->post_experimental_conditions ();
            this->experimental_condition (tmp);
          }

          count[1UL] = 1;
        }
      }
      else
      {
        assert (start);
        state = ~0UL;
      }
    }
    else if (n == "mesh" && ns.empty ())
    {
      if (count[2UL] == 0)
      {
        if (start)
        {
          ::mesh::mesh_pskel* p = 0;

          if (t == 0 && this->mesh_parser_ != 0)
            p = this->mesh_parser_;
          else
          {
            const char* ts = ::mesh::mesh_pskel::_static_type ();

            if (t == 0)
              t = ts;

            if (this->mesh_parser_ != 0 && strcmp (t, ts) == 0)
              p = this->mesh_parser_;
            else
            {
              if (t != ts &&
                  !::xsde::cxx::parser::validating::inheritance_map_instance ().check (t, ts))
              {
                ctx.schema_error (::xsde::cxx::schema_error::not_derived);
                return;
              }

              if (this->mesh_parser_map_ != 0)
                p = static_cast< ::mesh::mesh_pskel* > (
                  this->mesh_parser_map_->find (t));
            }
          }

          if (p)
          {
            p->pre ();
            ctx.nested_parser (p);
          }
        }
        else
        {
          ::mesh::mesh_pskel* p =
          static_cast< ::mesh::mesh_pskel* > (ctx.nested_parser ());

          if (p != 0)
          {
            ::mesh::mesh* tmp = p->post_mesh ();
            this->mesh (tmp);
          }

          count[2UL] = 1;
        }
      }
      else
      {
        assert (start);
        state = ~0UL;
      }
    }
    else if (n == "data" && ns.empty ())
    {
      if (count[3UL] == 0)
      {
        if (start)
        {
          ::variables::data_pskel* p = 0;

          if (t == 0 && this->data_parser_ != 0)
            p = this->data_parser_;
          else
          {
            const char* ts = ::variables::data_pskel::_static_type ();

            if (t == 0)
              t = ts;

            if (this->data_parser_ != 0 && strcmp (t, ts) == 0)
              p = this->data_parser_;
            else
            {
              if (t != ts &&
                  !::xsde::cxx::parser::validating::inheritance_map_instance ().check (t, ts))
              {
                ctx.schema_error (::xsde::cxx::schema_error::not_derived);
                return;
              }

              if (this->data_parser_map_ != 0)
                p = static_cast< ::variables::data_pskel* > (
                  this->data_parser_map_->find (t));
            }
          }

          if (p)
          {
            p->pre ();
            ctx.nested_parser (p);
          }
        }
        else
        {
          ::variables::data_pskel* p =
          static_cast< ::variables::data_pskel* > (ctx.nested_parser ());

          if (p != 0)
          {
            ::variables::data* tmp = p->post_data ();
            this->data (tmp);
          }

          count[3UL] = 1;
        }
      }
      else
      {
        assert (start);
        state = ~0UL;
      }
    }
    else if (n == "custom" && ns.empty ())
    {
      if (count[4UL] == 0)
      {
        if (start)
        {
          ::common::custom_pskel* p = 0;

          if (t == 0 && this->custom_parser_ != 0)
            p = this->custom_parser_;
          else
          {
            const char* ts = ::common::custom_pskel::_static_type ();

            if (t == 0)
              t = ts;

            if (this->custom_parser_ != 0 && strcmp (t, ts) == 0)
              p = this->custom_parser_;
            else
            {
              if (t != ts &&
                  !::xsde::cxx::parser::validating::inheritance_map_instance ().check (t, ts))
              {
                ctx.schema_error (::xsde::cxx::schema_error::not_derived);
                return;
              }

              if (this->custom_parser_map_ != 0)
                p = static_cast< ::common::custom_pskel* > (
                  this->custom_parser_map_->find (t));
            }
          }

          if (p)
          {
            p->pre ();
            ctx.nested_parser (p);
          }
        }
        else
        {
          ::common::custom_pskel* p =
          static_cast< ::common::custom_pskel* > (ctx.nested_parser ());

          if (p != 0)
          {
            ::common::custom* tmp = p->post_custom ();
            this->custom (tmp);
          }

          count[4UL] = 1;
        }
      }
      else
      {
        assert (start);
        state = ~0UL;
      }
    }
    else if (n.empty () && ns.empty ())
    {
      state = ~0UL;
    }
    else
      state = ~0UL;
  }

  // Element validation and dispatch functions for microenvironment_pskel.
  //
  bool microenvironment_pskel::
  _start_element_impl (const ::xsde::cxx::ro_string& ns,
                       const ::xsde::cxx::ro_string& n,
                       const char* t)
  {
    XSDE_UNUSED (t);

    ::xsde::cxx::parser::context& ctx = this->_context ();

    v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
    v_state_descr_* vd = vs.data + (vs.size - 1);

    if (vd->func == 0 && vd->state == 0)
    {
      typedef ::xsde::cxx::parser::validating::complex_content base;
      if (base::_start_element_impl (ns, n, t))
        return true;
      else
        vd->state = 1;
    }

    while (vd->func != 0)
    {
      (this->*vd->func) (vd->state, vd->count, ns, n, t, true);

      vd = vs.data + (vs.size - 1);

      if (vd->state == ~0UL && !ctx.error_type ())
        vd = vs.data + (--vs.size - 1);
      else
        break;
    }

    if (vd->func == 0)
    {
      if (vd->state != ~0UL)
      {
        unsigned long s = ~0UL;

        if (n == "domain" && ns.empty ())
          s = 0UL;
        else if (n == "vascular_network" && ns.empty ())
          s = 1UL;
        else if (n == "basement_membrane" && ns.empty ())
          s = 2UL;
        else if (n == "custom" && ns.empty ())
          s = 3UL;

        if (s != ~0UL)
        {
          vd->count++;
          vd->state = ~0UL;

          vd = vs.data + vs.size++;
          vd->func = &microenvironment_pskel::sequence_0;
          vd->state = s;
          vd->count = 0;

          this->sequence_0 (vd->state, vd->count, ns, n, t, true);
        }
        else
        {
          return false;
        }
      }
      else
        return false;
    }

    return true;
  }

  bool microenvironment_pskel::
  _end_element_impl (const ::xsde::cxx::ro_string& ns,
                     const ::xsde::cxx::ro_string& n)
  {
    v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
    v_state_descr_& vd = vs.data[vs.size - 1];

    if (vd.func == 0 && vd.state == 0)
    {
      typedef ::xsde::cxx::parser::validating::complex_content base;
      if (!base::_end_element_impl (ns, n))
        assert (false);
      return true;
    }

    assert (vd.func != 0);
    (this->*vd.func) (vd.state, vd.count, ns, n, 0, false);

    if (vd.state == ~0UL)
      vs.size--;

    return true;
  }

  void microenvironment_pskel::
  _pre_e_validate ()
  {
    this->v_state_stack_.push ();
    static_cast< v_state_* > (this->v_state_stack_.top ())->size = 0;

    v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
    v_state_descr_& vd = vs.data[vs.size++];

    vd.func = 0;
    vd.state = 0;
    vd.count = 0;
  }

  void microenvironment_pskel::
  _post_e_validate ()
  {
    ::xsde::cxx::parser::context& ctx = this->_context ();

    v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
    v_state_descr_* vd = vs.data + (vs.size - 1);

    ::xsde::cxx::ro_string empty;
    while (vd->func != 0)
    {
      (this->*vd->func) (vd->state, vd->count, empty, empty, 0, true);

      if (ctx.error_type ())
        return;

      assert (vd->state == ~0UL);
      vd = vs.data + (--vs.size - 1);
    }


    this->v_state_stack_.pop ();
  }

  void microenvironment_pskel::
  sequence_0 (unsigned long& state,
              unsigned long& count,
              const ::xsde::cxx::ro_string& ns,
              const ::xsde::cxx::ro_string& n,
              const char* t,
              bool start)
  {
    ::xsde::cxx::parser::context& ctx = this->_context ();

    XSDE_UNUSED (t);
    XSDE_UNUSED (ctx);

    switch (state)
    {
      case 0UL:
      {
        if (n == "domain" && ns.empty ())
        {
          if (start)
          {
            ::microenvironment::domain_pskel* p = 0;

            if (t == 0 && this->domain_parser_ != 0)
              p = this->domain_parser_;
            else
            {
              const char* ts = ::microenvironment::domain_pskel::_static_type ();

              if (t == 0)
                t = ts;

              if (this->domain_parser_ != 0 && strcmp (t, ts) == 0)
                p = this->domain_parser_;
              else
              {
                if (t != ts &&
                    !::xsde::cxx::parser::validating::inheritance_map_instance ().check (t, ts))
                {
                  ctx.schema_error (::xsde::cxx::schema_error::not_derived);
                  return;
                }

                if (this->domain_parser_map_ != 0)
                  p = static_cast< ::microenvironment::domain_pskel* > (
                    this->domain_parser_map_->find (t));
              }
            }

            if (p)
            {
              p->pre ();
              ctx.nested_parser (p);
            }
          }
          else
          {
            ::microenvironment::domain_pskel* p =
            static_cast< ::microenvironment::domain_pskel* > (ctx.nested_parser ());

            if (p != 0)
            {
              ::microenvironment::domain* tmp = p->post_domain ();
              this->domain (tmp);
            }

            count++;
          }

          break;
        }
        else
        {
          assert (start);
          count = 0;
          state = 1UL;
          // Fall through.
        }
      }
      case 1UL:
      {
        if (n == "vascular_network" && ns.empty ())
        {
          if (start)
          {
            ::vascular::vascular_network_pskel* p = 0;

            if (t == 0 && this->vascular_network_parser_ != 0)
              p = this->vascular_network_parser_;
            else
            {
              const char* ts = ::vascular::vascular_network_pskel::_static_type ();

              if (t == 0)
                t = ts;

              if (this->vascular_network_parser_ != 0 && strcmp (t, ts) == 0)
                p = this->vascular_network_parser_;
              else
              {
                if (t != ts &&
                    !::xsde::cxx::parser::validating::inheritance_map_instance ().check (t, ts))
                {
                  ctx.schema_error (::xsde::cxx::schema_error::not_derived);
                  return;
                }

                if (this->vascular_network_parser_map_ != 0)
                  p = static_cast< ::vascular::vascular_network_pskel* > (
                    this->vascular_network_parser_map_->find (t));
              }
            }

            if (p)
            {
              p->pre ();
              ctx.nested_parser (p);
            }
          }
          else
          {
            ::vascular::vascular_network_pskel* p =
            static_cast< ::vascular::vascular_network_pskel* > (ctx.nested_parser ());

            if (p != 0)
            {
              ::vascular::vascular_network* tmp = p->post_vascular_network ();
              this->vascular_network (tmp);
            }

            count++;
          }

          break;
        }
        else
        {
          assert (start);
          count = 0;
          state = 2UL;
          // Fall through.
        }
      }
      case 2UL:
      {
        if (n == "basement_membrane" && ns.empty ())
        {
          if (start)
          {
            ::basement::basement_membrane_pskel* p = 0;

            if (t == 0 && this->basement_membrane_parser_ != 0)
              p = this->basement_membrane_parser_;
            else
            {
              const char* ts = ::basement::basement_membrane_pskel::_static_type ();

              if (t == 0)
                t = ts;

              if (this->basement_membrane_parser_ != 0 && strcmp (t, ts) == 0)
                p = this->basement_membrane_parser_;
              else
              {
                if (t != ts &&
                    !::xsde::cxx::parser::validating::inheritance_map_instance ().check (t, ts))
                {
                  ctx.schema_error (::xsde::cxx::schema_error::not_derived);
                  return;
                }

                if (this->basement_membrane_parser_map_ != 0)
                  p = static_cast< ::basement::basement_membrane_pskel* > (
                    this->basement_membrane_parser_map_->find (t));
              }
            }

            if (p)
            {
              p->pre ();
              ctx.nested_parser (p);
            }
          }
          else
          {
            ::basement::basement_membrane_pskel* p =
            static_cast< ::basement::basement_membrane_pskel* > (ctx.nested_parser ());

            if (p != 0)
            {
              ::basement::basement_membrane* tmp = p->post_basement_membrane ();
              this->basement_membrane (tmp);
            }

            count++;
          }

          break;
        }
        else
        {
          assert (start);
          count = 0;
          state = 3UL;
          // Fall through.
        }
      }
      case 3UL:
      {
        if (n == "custom" && ns.empty ())
        {
          if (start)
          {
            ::common::custom_pskel* p = 0;

            if (t == 0 && this->custom_parser_ != 0)
              p = this->custom_parser_;
            else
            {
              const char* ts = ::common::custom_pskel::_static_type ();

              if (t == 0)
                t = ts;

              if (this->custom_parser_ != 0 && strcmp (t, ts) == 0)
                p = this->custom_parser_;
              else
              {
                if (t != ts &&
                    !::xsde::cxx::parser::validating::inheritance_map_instance ().check (t, ts))
                {
                  ctx.schema_error (::xsde::cxx::schema_error::not_derived);
                  return;
                }

                if (this->custom_parser_map_ != 0)
                  p = static_cast< ::common::custom_pskel* > (
                    this->custom_parser_map_->find (t));
              }
            }

            if (p)
            {
              p->pre ();
              ctx.nested_parser (p);
            }
          }
          else
          {
            ::common::custom_pskel* p =
            static_cast< ::common::custom_pskel* > (ctx.nested_parser ());

            if (p != 0)
            {
              ::common::custom* tmp = p->post_custom ();
              this->custom (tmp);
            }

            count = 0;
            state = ~0UL;
          }

          break;
        }
        else
        {
          assert (start);
          count = 0;
          state = ~0UL;
          // Fall through.
        }
      }
      case ~0UL:
        break;
    }
  }
}

namespace microenvironment
{
  // Attribute validation and dispatch functions for domain_pskel.
  //
  bool domain_pskel::
  _attribute_impl_phase_one (const ::xsde::cxx::ro_string& ns,
                             const ::xsde::cxx::ro_string& n,
                             const ::xsde::cxx::ro_string& s)
  {
    ::xsde::cxx::parser::context& ctx = this->_context ();

    if (n == "name" && ns.empty ())
    {
      if (this->name_parser_)
      {
        this->name_parser_->pre ();

        this->name_parser_->_pre_impl (ctx);

        if (!ctx.error_type ())
          this->name_parser_->_characters (s);

        if (!ctx.error_type ())
          this->name_parser_->_post_impl ();

        if (!ctx.error_type ())
        {
          const ::std::string& tmp = this->name_parser_->post_string ();

          this->name (tmp);
        }
      }

      return true;
    }

    return false;
  }
}

namespace microenvironment
{
}

// Begin epilogue.
//
//
// End epilogue.

