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

namespace phenotype_dataset
{
  // phenotype_dataset_sskel
  //

  inline
  void phenotype_dataset_sskel::
  keywords_serializer (::xml_schema::string_sskel& keywords)
  {
    this->keywords_serializer_ = &keywords;
  }

  inline
  void phenotype_dataset_sskel::
  ID_serializer (::xml_schema::unsigned_long_sskel& ID)
  {
    this->ID_serializer_ = &ID;
  }

  inline
  void phenotype_dataset_sskel::
  microenvironment_serializer (::microenvironment::microenvironment_sskel& s)
  {
    this->microenvironment_serializer_ = &s;
  }

  inline
  void phenotype_dataset_sskel::
  microenvironment_serializer (::xml_schema::serializer_map& m)
  {
    this->microenvironment_serializer_map_ = &m;
  }

  inline
  void phenotype_dataset_sskel::
  phenotype_serializer (::phenotype::phenotype_sskel& s)
  {
    this->phenotype_serializer_ = &s;
  }

  inline
  void phenotype_dataset_sskel::
  phenotype_serializer (::xml_schema::serializer_map& m)
  {
    this->phenotype_serializer_map_ = &m;
  }

  inline
  void phenotype_dataset_sskel::
  cell_part_serializer (::phenotype_base::cell_parts_sskel& s)
  {
    this->cell_part_serializer_ = &s;
  }

  inline
  void phenotype_dataset_sskel::
  cell_part_serializer (::xml_schema::serializer_map& m)
  {
    this->cell_part_serializer_map_ = &m;
  }

  inline
  void phenotype_dataset_sskel::
  custom_serializer (::common::custom_sskel& s)
  {
    this->custom_serializer_ = &s;
  }

  inline
  void phenotype_dataset_sskel::
  custom_serializer (::xml_schema::serializer_map& m)
  {
    this->custom_serializer_map_ = &m;
  }

  inline
  void phenotype_dataset_sskel::
  serializers (::xml_schema::string_sskel& keywords,
               ::xml_schema::unsigned_long_sskel& ID,
               ::microenvironment::microenvironment_sskel& microenvironment,
               ::phenotype::phenotype_sskel& phenotype,
               ::phenotype_base::cell_parts_sskel& cell_part,
               ::common::custom_sskel& custom)
  {
    this->keywords_serializer_ = &keywords;
    this->ID_serializer_ = &ID;
    this->microenvironment_serializer_ = &microenvironment;
    this->phenotype_serializer_ = &phenotype;
    this->cell_part_serializer_ = &cell_part;
    this->custom_serializer_ = &custom;
  }

  inline
  void phenotype_dataset_sskel::
  serializer_maps (::xml_schema::serializer_map& microenvironment,
                   ::xml_schema::serializer_map& phenotype,
                   ::xml_schema::serializer_map& cell_part,
                   ::xml_schema::serializer_map& custom)
  {
    this->microenvironment_serializer_map_ = &microenvironment;
    this->phenotype_serializer_map_ = &phenotype;
    this->cell_part_serializer_map_ = &cell_part;
    this->custom_serializer_map_ = &custom;
  }

  inline
  phenotype_dataset_sskel::
  phenotype_dataset_sskel ()
  : phenotype_dataset_impl_ (0),
    keywords_serializer_ (0),
    ID_serializer_ (0),
    microenvironment_serializer_ (0),
    microenvironment_serializer_map_ (0),
    phenotype_serializer_ (0),
    phenotype_serializer_map_ (0),
    cell_part_serializer_ (0),
    cell_part_serializer_map_ (0),
    custom_serializer_ (0),
    custom_serializer_map_ (0)
  {
  }

  inline
  phenotype_dataset_sskel::
  phenotype_dataset_sskel (phenotype_dataset_sskel* impl, void*)
  : ::xsde::cxx::serializer::non_validating::complex_content (impl, 0),
    phenotype_dataset_impl_ (impl),
    keywords_serializer_ (0),
    ID_serializer_ (0),
    microenvironment_serializer_ (0),
    microenvironment_serializer_map_ (0),
    phenotype_serializer_ (0),
    phenotype_serializer_map_ (0),
    cell_part_serializer_ (0),
    cell_part_serializer_map_ (0),
    custom_serializer_ (0),
    custom_serializer_map_ (0)
  {
  }
}

// Begin epilogue.
//
//
// End epilogue.

