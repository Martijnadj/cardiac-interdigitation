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

#ifndef BASEMENT_MEMBRANE_SIMPL_HPP
#define BASEMENT_MEMBRANE_SIMPL_HPP

#include <xsde/cxx/pre.hxx>

// Begin prologue.
//
//
// End prologue.

#ifndef XSDE_OMIT_SAGGR
#  define XSDE_OMIT_SAGGR
#  define BASEMENT_MEMBRANE_SIMPL_HPP_CLEAR_OMIT_SAGGR
#endif

#include "basement_membrane-sskel.hpp"

#include <xsde/cxx/stack.hxx>

#include "common-simpl.hpp"

#include "mesh-simpl.hpp"

namespace basement
{
  class basement_edge_simpl: public basement_edge_sskel
  {
    public:
    basement_edge_simpl ();

    virtual void
    pre (const ::basement::basement_edge&);

    // Elements.
    //
    virtual const ::common::units_decimal_nonnegative&
    tensile_strength ();

    virtual const ::common::custom&
    custom ();

    public:
    ::mesh::edge_simpl base_impl_;

    public:
    struct basement_edge_simpl_state
    {
      const ::basement::basement_edge* basement_edge_;
    };

    basement_edge_simpl_state basement_edge_simpl_state_;
  };

  class basement_face_simpl: public basement_face_sskel
  {
    public:
    basement_face_simpl ();

    virtual void
    pre (const ::basement::basement_face&);

    // Elements.
    //
    virtual const ::common::units_decimal_nonnegative&
    thickness ();

    virtual const ::common::custom&
    custom ();

    public:
    ::mesh::face_simpl base_impl_;

    public:
    struct basement_face_simpl_state
    {
      const ::basement::basement_face* basement_face_;
    };

    basement_face_simpl_state basement_face_simpl_state_;
  };

  class nodes_simpl: public nodes_sskel
  {
    public:
    nodes_simpl ();

    virtual void
    pre (const ::basement::nodes&);

    // Elements.
    //
    virtual const ::mesh::node&
    node ();

    virtual const ::common::custom&
    custom ();

    virtual void
    post ();

    virtual void
    _reset ();

    public:
    struct nodes_simpl_state
    {
      const ::basement::nodes* nodes_;
    };

    nodes_simpl_state nodes_simpl_state_first_;
    ::xsde::cxx::stack nodes_simpl_state_;
  };

  class egdes_simpl: public egdes_sskel
  {
    public:
    virtual void
    pre (const ::basement::egdes&);

    // Elements.
    //
    virtual const ::basement::basement_edge&
    edge ();

    virtual const ::common::custom&
    custom ();

    public:
    struct egdes_simpl_state
    {
      const ::basement::egdes* egdes_;
    };

    egdes_simpl_state egdes_simpl_state_;
  };

  class faces_simpl: public faces_sskel
  {
    public:
    virtual void
    pre (const ::basement::faces&);

    // Elements.
    //
    virtual const ::basement::basement_face&
    face ();

    virtual const ::common::custom&
    custom ();

    public:
    struct faces_simpl_state
    {
      const ::basement::faces* faces_;
    };

    faces_simpl_state faces_simpl_state_;
  };

  class basement_membrane_simpl: public basement_membrane_sskel
  {
    public:
    virtual void
    pre (const ::basement::basement_membrane&);

    // Attributes.
    //
    virtual bool
    ID_present ();

    virtual unsigned int
    ID ();

    // Elements.
    //
    virtual const ::basement::nodes&
    nodes ();

    virtual const ::basement::egdes&
    edges ();

    virtual const ::basement::faces&
    faces ();

    virtual const ::common::custom&
    custom ();

    public:
    struct basement_membrane_simpl_state
    {
      const ::basement::basement_membrane* basement_membrane_;
    };

    basement_membrane_simpl_state basement_membrane_simpl_state_;
  };
}

#ifdef BASEMENT_MEMBRANE_SIMPL_HPP_CLEAR_OMIT_SAGGR
#  undef XSDE_OMIT_SAGGR
#endif

#ifndef XSDE_OMIT_SAGGR

#endif // XSDE_OMIT_SAGGR

// Begin epilogue.
//
//
// End epilogue.

#include <xsde/cxx/post.hxx>

#endif // BASEMENT_MEMBRANE_SIMPL_HPP
