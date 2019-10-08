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

#ifndef PHENOTYPE_DATASET_SIMPL_HPP
#define PHENOTYPE_DATASET_SIMPL_HPP

#include <xsde/cxx/pre.hxx>

// Begin prologue.
//
//
// End prologue.

#ifndef XSDE_OMIT_SAGGR
#  define XSDE_OMIT_SAGGR
#  define PHENOTYPE_DATASET_SIMPL_HPP_CLEAR_OMIT_SAGGR
#endif

#include "phenotype_dataset-sskel.hpp"

#include <xsde/cxx/stack.hxx>

#include "common-simpl.hpp"

#include "microenvironment-simpl.hpp"

#include "phenotype-simpl.hpp"

#include "phenotype_base-simpl.hpp"

namespace phenotype_dataset
{
  class phenotype_dataset_simpl: public phenotype_dataset_sskel
  {
    public:
    virtual void
    pre (const ::phenotype_dataset::phenotype_dataset&);

    // Attributes.
    //
    virtual bool
    keywords_present ();

    virtual ::std::string
    keywords ();

    virtual bool
    ID_present ();

    virtual unsigned long long
    ID ();

    // Elements.
    //
    virtual bool
    microenvironment_present ();

    virtual const ::microenvironment::microenvironment&
    microenvironment ();

    virtual bool
    phenotype_next ();

    virtual const ::phenotype::phenotype&
    phenotype ();

    virtual bool
    cell_part_next ();

    virtual const ::phenotype_base::cell_parts&
    cell_part ();

    virtual bool
    custom_present ();

    virtual const ::common::custom&
    custom ();

    public:
    struct phenotype_dataset_simpl_state
    {
      const ::phenotype_dataset::phenotype_dataset* phenotype_dataset_;
      ::phenotype_dataset::phenotype_dataset::phenotype_const_iterator phenotype_;
      ::phenotype_dataset::phenotype_dataset::phenotype_const_iterator phenotype_end_;
      ::phenotype_dataset::phenotype_dataset::cell_part_const_iterator cell_part_;
      ::phenotype_dataset::phenotype_dataset::cell_part_const_iterator cell_part_end_;
    };

    phenotype_dataset_simpl_state phenotype_dataset_simpl_state_;
  };
}

#ifdef PHENOTYPE_DATASET_SIMPL_HPP_CLEAR_OMIT_SAGGR
#  undef XSDE_OMIT_SAGGR
#endif

#ifndef XSDE_OMIT_SAGGR

#endif // XSDE_OMIT_SAGGR

// Begin epilogue.
//
//
// End epilogue.

#include <xsde/cxx/post.hxx>

#endif // PHENOTYPE_DATASET_SIMPL_HPP
