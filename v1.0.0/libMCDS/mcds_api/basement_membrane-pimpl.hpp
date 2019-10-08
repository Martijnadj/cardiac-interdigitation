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

#ifndef BASEMENT_MEMBRANE_PIMPL_HPP
#define BASEMENT_MEMBRANE_PIMPL_HPP

#include <xsde/cxx/pre.hxx>

// Begin prologue.
//
//
// End prologue.

#ifndef XSDE_OMIT_PAGGR
#  define XSDE_OMIT_PAGGR
#  define BASEMENT_MEMBRANE_PIMPL_HPP_CLEAR_OMIT_PAGGR
#endif

#include "basement_membrane-pskel.hpp"

#include <xsde/cxx/stack.hxx>

#include "common-pimpl.hpp"

#include "mesh-pimpl.hpp"

namespace basement
{
  class basement_edge_pimpl: public basement_edge_pskel
  {
    public:
    basement_edge_pimpl (bool = false);

    ~basement_edge_pimpl ();

    virtual void
    _reset ();

    virtual void
    pre ();

    // Elements.
    //
    virtual void
    tensile_strength (::common::units_decimal_nonnegative*);

    virtual void
    custom (::common::custom*);

    virtual ::basement::basement_edge*
    post_basement_edge ();

    public:
    void
    pre_impl (::basement::basement_edge*);

    public:
    ::mesh::edge_pimpl base_impl_;

    public:
    struct basement_edge_pimpl_state
    {
      ::basement::basement_edge* basement_edge_;
    };

    basement_edge_pimpl_state basement_edge_pimpl_state_;
    bool basement_edge_pimpl_base_;
  };

  class basement_face_pimpl: public basement_face_pskel
  {
    public:
    basement_face_pimpl (bool = false);

    ~basement_face_pimpl ();

    virtual void
    _reset ();

    virtual void
    pre ();

    // Elements.
    //
    virtual void
    thickness (::common::units_decimal_nonnegative*);

    virtual void
    custom (::common::custom*);

    virtual ::basement::basement_face*
    post_basement_face ();

    public:
    void
    pre_impl (::basement::basement_face*);

    public:
    ::mesh::face_pimpl base_impl_;

    public:
    struct basement_face_pimpl_state
    {
      ::basement::basement_face* basement_face_;
    };

    basement_face_pimpl_state basement_face_pimpl_state_;
    bool basement_face_pimpl_base_;
  };

  class nodes_pimpl: public nodes_pskel
  {
    public:
    nodes_pimpl (bool = false);

    ~nodes_pimpl ();

    virtual void
    _reset ();

    virtual void
    pre ();

    // Elements.
    //
    virtual void
    node (::mesh::node*);

    virtual void
    custom (::common::custom*);

    virtual ::basement::nodes*
    post_nodes ();

    public:
    void
    pre_impl (::basement::nodes*);

    public:
    struct nodes_pimpl_state
    {
      ::basement::nodes* nodes_;
    };

    nodes_pimpl_state nodes_pimpl_state_first_;
    ::xsde::cxx::stack nodes_pimpl_state_;
    bool nodes_pimpl_base_;
  };

  class egdes_pimpl: public egdes_pskel
  {
    public:
    egdes_pimpl (bool = false);

    ~egdes_pimpl ();

    virtual void
    _reset ();

    virtual void
    pre ();

    // Elements.
    //
    virtual void
    edge (::basement::basement_edge*);

    virtual void
    custom (::common::custom*);

    virtual ::basement::egdes*
    post_egdes ();

    public:
    void
    pre_impl (::basement::egdes*);

    public:
    struct egdes_pimpl_state
    {
      ::basement::egdes* egdes_;
    };

    egdes_pimpl_state egdes_pimpl_state_;
    bool egdes_pimpl_base_;
  };

  class faces_pimpl: public faces_pskel
  {
    public:
    faces_pimpl (bool = false);

    ~faces_pimpl ();

    virtual void
    _reset ();

    virtual void
    pre ();

    // Elements.
    //
    virtual void
    face (::basement::basement_face*);

    virtual void
    custom (::common::custom*);

    virtual ::basement::faces*
    post_faces ();

    public:
    void
    pre_impl (::basement::faces*);

    public:
    struct faces_pimpl_state
    {
      ::basement::faces* faces_;
    };

    faces_pimpl_state faces_pimpl_state_;
    bool faces_pimpl_base_;
  };

  class basement_membrane_pimpl: public basement_membrane_pskel
  {
    public:
    basement_membrane_pimpl (bool = false);

    ~basement_membrane_pimpl ();

    virtual void
    _reset ();

    virtual void
    pre ();

    // Attributes.
    //
    virtual void
    ID (unsigned int);

    // Elements.
    //
    virtual void
    nodes (::basement::nodes*);

    virtual void
    edges (::basement::egdes*);

    virtual void
    faces (::basement::faces*);

    virtual void
    custom (::common::custom*);

    virtual ::basement::basement_membrane*
    post_basement_membrane ();

    public:
    void
    pre_impl (::basement::basement_membrane*);

    public:
    struct basement_membrane_pimpl_state
    {
      ::basement::basement_membrane* basement_membrane_;
    };

    basement_membrane_pimpl_state basement_membrane_pimpl_state_;
    bool basement_membrane_pimpl_base_;
  };
}

#ifdef BASEMENT_MEMBRANE_PIMPL_HPP_CLEAR_OMIT_PAGGR
#  undef XSDE_OMIT_PAGGR
#endif

#ifndef XSDE_OMIT_PAGGR

#endif // XSDE_OMIT_PAGGR

// Begin epilogue.
//
//
// End epilogue.

#include <xsde/cxx/post.hxx>

#endif // BASEMENT_MEMBRANE_PIMPL_HPP
