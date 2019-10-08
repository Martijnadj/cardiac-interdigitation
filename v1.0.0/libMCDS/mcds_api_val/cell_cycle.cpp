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

#include <xsde/cxx/pre.hxx>

// Begin prologue.
//
//
// End prologue.

#include "cell_cycle.hpp"

#include <stdlib.h>
#include <new>

#include <xsde/cxx/guard.hxx>

namespace cell_cycle
{
  // death_type
  //

  static const char* _xsde_death_type_enumerators_[] = 
  {
    "apoptosis",
    "necrosis",
    "autophagy"
  };

  const char* death_type::
  string () const
  {
    return _xsde_death_type_enumerators_[value_];
  }

  // death_rate_type
  //

  death_rate_type::
  death_rate_type ()
  {
  }

  death_rate_type::
  ~death_rate_type ()
  {
  }

  void death_rate_type::
  _copy (death_rate_type& c) const
  {
    XSDE_UNUSED (c);

    const ::common::units_decimal_nonnegative& b = *this;
    b._copy (c);
    c.type (this->type ());
  }

  death_rate_type* death_rate_type::
  _clone () const
  {
    death_rate_type* c = new death_rate_type;
    ::xsde::cxx::guard< death_rate_type > g (c);
    this->_copy (*c);
    g.release ();
    return c;
  }

  // cell_cycle_arrest
  //

  cell_cycle_arrest::
  cell_cycle_arrest ()
  {
    this->condition_ = 0;
  }

  cell_cycle_arrest::
  ~cell_cycle_arrest ()
  {
    delete this->condition_;
  }

  void cell_cycle_arrest::
  _copy (cell_cycle_arrest& c) const
  {
    XSDE_UNUSED (c);

    if (this->condition_present ())
    {
      ::cell_cycle::arrest_condition* m = this->condition ()._clone ();
      c.condition (m);
    }
  }

  cell_cycle_arrest* cell_cycle_arrest::
  _clone () const
  {
    cell_cycle_arrest* c = new cell_cycle_arrest;
    ::xsde::cxx::guard< cell_cycle_arrest > g (c);
    this->_copy (*c);
    g.release ();
    return c;
  }

  // transition
  //

  transition::
  transition ()
  {
    this->checkpoint_failure_probability_ = 0;
    this->subsequent_phase_present_ = false;
    this->transition_rate_ = 0;
  }

  transition::
  ~transition ()
  {
    delete this->checkpoint_failure_probability_;
    delete this->transition_rate_;
  }

  void transition::
  _copy (transition& c) const
  {
    XSDE_UNUSED (c);

    if (this->checkpoint_failure_probability_present ())
    {
      ::common::units_decimal* m = this->checkpoint_failure_probability ()._clone ();
      c.checkpoint_failure_probability (m);
    }

    if (this->subsequent_phase_present ())
      c.subsequent_phase (this->subsequent_phase ());

    this->threshold ().copy (c.threshold ());

    if (this->transition_rate_present ())
    {
      ::common::units_decimal* m = this->transition_rate ()._clone ();
      c.transition_rate (m);
    }
  }

  transition* transition::
  _clone () const
  {
    transition* c = new transition;
    ::xsde::cxx::guard< transition > g (c);
    this->_copy (*c);
    g.release ();
    return c;
  }

  // cell_cycle_phase
  //

  cell_cycle_phase::
  cell_cycle_phase ()
  {
    this->ID_present_ = false;
    this->birth_rate_ = 0;
    this->duration_ = 0;
    this->net_birth_rate_ = 0;
    this->population_doubling_time_ = 0;
    this->cell_cycle_arrest_ = 0;
    this->custom_ = 0;
  }

  cell_cycle_phase::
  ~cell_cycle_phase ()
  {
    delete this->birth_rate_;
    delete this->duration_;
    delete this->net_birth_rate_;
    delete this->population_doubling_time_;
    delete this->cell_cycle_arrest_;
    delete this->custom_;
  }

  void cell_cycle_phase::
  _copy (cell_cycle_phase& c) const
  {
    XSDE_UNUSED (c);

    c.name (this->name ());

    if (this->ID_present ())
      c.ID (this->ID ());

    if (this->birth_rate_present ())
    {
      ::common::units_decimal_nonnegative* m = this->birth_rate ()._clone ();
      c.birth_rate (m);
    }

    if (this->duration_present ())
    {
      ::common::units_decimal_nonnegative* m = this->duration ()._clone ();
      c.duration (m);
    }

    this->death_rate ().copy (c.death_rate ());

    if (this->net_birth_rate_present ())
    {
      ::common::units_decimal* m = this->net_birth_rate ()._clone ();
      c.net_birth_rate (m);
    }

    if (this->population_doubling_time_present ())
    {
      ::common::units_decimal_nonnegative* m = this->population_doubling_time ()._clone ();
      c.population_doubling_time (m);
    }

    if (this->cell_cycle_arrest_present ())
    {
      ::cell_cycle::cell_cycle_arrest* m = this->cell_cycle_arrest ()._clone ();
      c.cell_cycle_arrest (m);
    }

    this->transition ().copy (c.transition ());

    this->cell_part ().copy (c.cell_part ());

    if (this->custom_present ())
    {
      ::common::custom* m = this->custom ()._clone ();
      c.custom (m);
    }
  }

  cell_cycle_phase* cell_cycle_phase::
  _clone () const
  {
    cell_cycle_phase* c = new cell_cycle_phase;
    ::xsde::cxx::guard< cell_cycle_phase > g (c);
    this->_copy (*c);
    g.release ();
    return c;
  }

  // summary_elements
  //

  summary_elements::
  summary_elements ()
  {
    this->birth_rate_ = 0;
    this->duration_ = 0;
    this->net_birth_rate_ = 0;
    this->population_doubling_time_ = 0;
  }

  summary_elements::
  ~summary_elements ()
  {
    delete this->birth_rate_;
    delete this->duration_;
    delete this->net_birth_rate_;
    delete this->population_doubling_time_;
  }

  void summary_elements::
  _copy (summary_elements& c) const
  {
    XSDE_UNUSED (c);

    if (this->birth_rate_present ())
    {
      ::common::units_decimal_nonnegative* m = this->birth_rate ()._clone ();
      c.birth_rate (m);
    }

    if (this->duration_present ())
    {
      ::common::units_decimal_nonnegative* m = this->duration ()._clone ();
      c.duration (m);
    }

    this->death_rate ().copy (c.death_rate ());

    if (this->net_birth_rate_present ())
    {
      ::common::units_decimal* m = this->net_birth_rate ()._clone ();
      c.net_birth_rate (m);
    }

    if (this->population_doubling_time_present ())
    {
      ::common::units_decimal_nonnegative* m = this->population_doubling_time ()._clone ();
      c.population_doubling_time (m);
    }
  }

  summary_elements* summary_elements::
  _clone () const
  {
    summary_elements* c = new summary_elements;
    ::xsde::cxx::guard< summary_elements > g (c);
    this->_copy (*c);
    g.release ();
    return c;
  }

  // cell_cycle
  //

  cell_cycle::
  cell_cycle ()
  {
    this->ID_present_ = false;
    this->summary_elements_ = 0;
    this->custom_ = 0;
  }

  cell_cycle::
  ~cell_cycle ()
  {
    delete this->summary_elements_;
    delete this->custom_;
  }

  void cell_cycle::
  _copy (cell_cycle& c) const
  {
    XSDE_UNUSED (c);

    c.model (this->model ());

    if (this->ID_present ())
      c.ID (this->ID ());

    this->cell_cycle_phase ().copy (c.cell_cycle_phase ());

    this->cell_death ().copy (c.cell_death ());

    if (this->summary_elements_present ())
    {
      ::cell_cycle::summary_elements* m = this->summary_elements ()._clone ();
      c.summary_elements (m);
    }

    if (this->custom_present ())
    {
      ::common::custom* m = this->custom ()._clone ();
      c.custom (m);
    }
  }

  cell_cycle* cell_cycle::
  _clone () const
  {
    cell_cycle* c = new cell_cycle;
    ::xsde::cxx::guard< cell_cycle > g (c);
    this->_copy (*c);
    g.release ();
    return c;
  }

  // cell_death
  //

  cell_death::
  cell_death ()
  {
    this->ID_present_ = false;
    this->duration_ = 0;
    this->custom_ = 0;
  }

  cell_death::
  ~cell_death ()
  {
    delete this->duration_;
    delete this->custom_;
  }

  void cell_death::
  _copy (cell_death& c) const
  {
    XSDE_UNUSED (c);

    c.type (this->type ());

    if (this->ID_present ())
      c.ID (this->ID ());

    {
      ::common::units_decimal* m = this->duration ()._clone ();
      c.duration (m);
    }

    this->cell_part ().copy (c.cell_part ());

    if (this->custom_present ())
    {
      ::common::custom* m = this->custom ()._clone ();
      c.custom (m);
    }
  }

  cell_death* cell_death::
  _clone () const
  {
    cell_death* c = new cell_death;
    ::xsde::cxx::guard< cell_death > g (c);
    this->_copy (*c);
    g.release ();
    return c;
  }

  // arrest_type
  //

  static const char* _xsde_arrest_type_enumerators_[] = 
  {
    "maximum_cell_density",
    "maximum_cell_surface_density",
    "maximum_cell_volume_density",
    "maximum_cell_number",
    "maximum_volume_fraction",
    "maximum_area_fraction"
  };

  const char* arrest_type::
  string () const
  {
    return _xsde_arrest_type_enumerators_[value_];
  }

  // arrest_condition
  //

  arrest_condition::
  arrest_condition ()
  {
    this->type_present_ = false;
  }

  arrest_condition::
  ~arrest_condition ()
  {
  }

  void arrest_condition::
  _copy (arrest_condition& c) const
  {
    XSDE_UNUSED (c);

    const ::common::units_decimal& b = *this;
    b._copy (c);
    if (this->type_present ())
      c.type (this->type ());
  }

  arrest_condition* arrest_condition::
  _clone () const
  {
    arrest_condition* c = new arrest_condition;
    ::xsde::cxx::guard< arrest_condition > g (c);
    this->_copy (*c);
    g.release ();
    return c;
  }

  // cycles_and_deaths
  //

  cycles_and_deaths::
  cycles_and_deaths ()
  {
  }

  cycles_and_deaths::
  ~cycles_and_deaths ()
  {
  }

  void cycles_and_deaths::
  _copy (cycles_and_deaths& c) const
  {
    XSDE_UNUSED (c);

    this->cell_cycle ().copy (c.cell_cycle ());

    this->cell_death ().copy (c.cell_death ());
  }

  cycles_and_deaths* cycles_and_deaths::
  _clone () const
  {
    cycles_and_deaths* c = new cycles_and_deaths;
    ::xsde::cxx::guard< cycles_and_deaths > g (c);
    this->_copy (*c);
    g.release ();
    return c;
  }
}

// Begin epilogue.
//
//
// End epilogue.

#include <xsde/cxx/post.hxx>

