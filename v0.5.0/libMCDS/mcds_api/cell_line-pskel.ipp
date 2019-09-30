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

namespace cell_line
{
  // cell_line_pskel
  //

  inline
  void cell_line_pskel::
  ID_parser (::xml_schema::string_pskel& p)
  {
    this->ID_parser_ = &p;
  }

  inline
  void cell_line_pskel::
  label_parser (::xml_schema::string_pskel& p)
  {
    this->label_parser_ = &p;
  }

  inline
  void cell_line_pskel::
  curated_parser (::xml_schema::boolean_pskel& p)
  {
    this->curated_parser_ = &p;
  }

  inline
  void cell_line_pskel::
  metadata_parser (::metadata::metadata_pskel& p)
  {
    this->metadata_parser_ = &p;
  }

  inline
  void cell_line_pskel::
  metadata_parser (::xml_schema::parser_map& m)
  {
    this->metadata_parser_map_ = &m;
  }

  inline
  void cell_line_pskel::
  phenotype_dataset_parser (::phenotype_dataset::phenotype_dataset_pskel& p)
  {
    this->phenotype_dataset_parser_ = &p;
  }

  inline
  void cell_line_pskel::
  phenotype_dataset_parser (::xml_schema::parser_map& m)
  {
    this->phenotype_dataset_parser_map_ = &m;
  }

  inline
  void cell_line_pskel::
  custom_parser (::common::custom_pskel& p)
  {
    this->custom_parser_ = &p;
  }

  inline
  void cell_line_pskel::
  custom_parser (::xml_schema::parser_map& m)
  {
    this->custom_parser_map_ = &m;
  }

  inline
  void cell_line_pskel::
  parsers (::xml_schema::string_pskel& ID,
           ::xml_schema::string_pskel& label,
           ::xml_schema::boolean_pskel& curated,
           ::metadata::metadata_pskel& metadata,
           ::phenotype_dataset::phenotype_dataset_pskel& phenotype_dataset,
           ::common::custom_pskel& custom)
  {
    this->ID_parser_ = &ID;
    this->label_parser_ = &label;
    this->curated_parser_ = &curated;
    this->metadata_parser_ = &metadata;
    this->phenotype_dataset_parser_ = &phenotype_dataset;
    this->custom_parser_ = &custom;
  }

  inline
  void cell_line_pskel::
  parser_maps (::xml_schema::parser_map& metadata,
               ::xml_schema::parser_map& phenotype_dataset,
               ::xml_schema::parser_map& custom)
  {
    this->metadata_parser_map_ = &metadata;
    this->phenotype_dataset_parser_map_ = &phenotype_dataset;
    this->custom_parser_map_ = &custom;
  }

  inline
  cell_line_pskel::
  cell_line_pskel ()
  : cell_line_impl_ (0),
    ID_parser_ (0),
    label_parser_ (0),
    curated_parser_ (0),
    metadata_parser_ (0),
    metadata_parser_map_ (0),
    phenotype_dataset_parser_ (0),
    phenotype_dataset_parser_map_ (0),
    custom_parser_ (0),
    custom_parser_map_ (0),
    v_state_stack_ (sizeof (v_state_), &v_state_first_)
  {
  }

  inline
  cell_line_pskel::
  cell_line_pskel (cell_line_pskel* impl, void*)
  : ::xsde::cxx::parser::validating::complex_content (impl, 0),
    cell_line_impl_ (impl),
    ID_parser_ (0),
    label_parser_ (0),
    curated_parser_ (0),
    metadata_parser_ (0),
    metadata_parser_map_ (0),
    phenotype_dataset_parser_ (0),
    phenotype_dataset_parser_map_ (0),
    custom_parser_ (0),
    custom_parser_map_ (0),
    v_state_stack_ (sizeof (v_state_), &v_state_first_)
  {
  }

  // DCLs_pskel
  //

  inline
  void DCLs_pskel::
  cell_line_parser (::cell_line::cell_line_pskel& p)
  {
    this->cell_line_parser_ = &p;
  }

  inline
  void DCLs_pskel::
  cell_line_parser (::xml_schema::parser_map& m)
  {
    this->cell_line_parser_map_ = &m;
  }

  inline
  void DCLs_pskel::
  parsers (::cell_line::cell_line_pskel& cell_line)
  {
    this->cell_line_parser_ = &cell_line;
  }

  inline
  void DCLs_pskel::
  parser_maps (::xml_schema::parser_map& cell_line)
  {
    this->cell_line_parser_map_ = &cell_line;
  }

  inline
  DCLs_pskel::
  DCLs_pskel ()
  : DCLs_impl_ (0),
    cell_line_parser_ (0),
    cell_line_parser_map_ (0),
    v_state_stack_ (sizeof (v_state_), &v_state_first_)
  {
  }

  inline
  DCLs_pskel::
  DCLs_pskel (DCLs_pskel* impl, void*)
  : ::xsde::cxx::parser::validating::complex_content (impl, 0),
    DCLs_impl_ (impl),
    cell_line_parser_ (0),
    cell_line_parser_map_ (0),
    v_state_stack_ (sizeof (v_state_), &v_state_first_)
  {
  }
}

// Begin epilogue.
//
//
// End epilogue.

