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

#include "cell_line-simpl.hpp"

#include <xsde/cxx/serializer/validating/string-common.hxx>

namespace cell_line
{
  // cell_line_simpl
  //

  void cell_line_simpl::
  pre (const ::cell_line::cell_line& x)
  {
    this->cell_line_simpl_state_.cell_line_ = &x;
    this->cell_line_simpl_state_.phenotype_dataset_ = 
    this->cell_line_simpl_state_.cell_line_->phenotype_dataset ().begin ();
    this->cell_line_simpl_state_.phenotype_dataset_end_ = 
    this->cell_line_simpl_state_.cell_line_->phenotype_dataset ().end ();
  }

  bool cell_line_simpl::
  ID_present ()
  {
    return this->cell_line_simpl_state_.cell_line_->ID_present ();
  }

  ::std::string cell_line_simpl::
  ID ()
  {
    return this->cell_line_simpl_state_.cell_line_->ID ();
  }

  bool cell_line_simpl::
  label_present ()
  {
    return this->cell_line_simpl_state_.cell_line_->label_present ();
  }

  ::std::string cell_line_simpl::
  label ()
  {
    return this->cell_line_simpl_state_.cell_line_->label ();
  }

  bool cell_line_simpl::
  curated_present ()
  {
    return this->cell_line_simpl_state_.cell_line_->curated_present ();
  }

  bool cell_line_simpl::
  curated ()
  {
    return this->cell_line_simpl_state_.cell_line_->curated ();
  }

  bool cell_line_simpl::
  metadata_present ()
  {
    return this->cell_line_simpl_state_.cell_line_->metadata_present ();
  }

  const ::metadata::metadata& cell_line_simpl::
  metadata ()
  {
    return this->cell_line_simpl_state_.cell_line_->metadata ();
  }

  bool cell_line_simpl::
  phenotype_dataset_next ()
  {
    return this->cell_line_simpl_state_.phenotype_dataset_ != 
    this->cell_line_simpl_state_.phenotype_dataset_end_;
  }

  const ::phenotype_dataset::phenotype_dataset& cell_line_simpl::
  phenotype_dataset ()
  {
    return *this->cell_line_simpl_state_.phenotype_dataset_++;
  }

  bool cell_line_simpl::
  custom_present ()
  {
    return this->cell_line_simpl_state_.cell_line_->custom_present ();
  }

  const ::common::custom& cell_line_simpl::
  custom ()
  {
    return this->cell_line_simpl_state_.cell_line_->custom ();
  }

  // DCLs_simpl
  //

  void DCLs_simpl::
  pre (const ::cell_line::DCLs& x)
  {
    this->DCLs_simpl_state_.DCLs_ = &x;
    this->DCLs_simpl_state_.cell_line_ = 
    this->DCLs_simpl_state_.DCLs_->cell_line ().begin ();
    this->DCLs_simpl_state_.cell_line_end_ = 
    this->DCLs_simpl_state_.DCLs_->cell_line ().end ();
  }

  bool DCLs_simpl::
  cell_line_next ()
  {
    return this->DCLs_simpl_state_.cell_line_ != 
    this->DCLs_simpl_state_.cell_line_end_;
  }

  const ::cell_line::cell_line& DCLs_simpl::
  cell_line ()
  {
    return *this->DCLs_simpl_state_.cell_line_++;
  }
}

// Begin epilogue.
//
//
// End epilogue.

