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

#ifndef PHENOTYPE_DATASET_FWD_HPP
#define PHENOTYPE_DATASET_FWD_HPP

#include <xsde/cxx/version.hxx>

#if (XSDE_INT_VERSION != 3020000L)
#error XSD/e runtime version mismatch
#endif

#include <xsde/cxx/pre.hxx>

// Begin prologue.
//
//
// End prologue.

#include <string>
#include <xsde/cxx/hybrid/xml-schema.hxx>
#include <xsde/cxx/hybrid/sequence.hxx>

namespace xml_schema
{
  using ::xsde::cxx::hybrid::any_type;
  typedef ::std::string any_simple_type;

  typedef signed char byte;
  using ::xsde::cxx::hybrid::byte_base;

  typedef unsigned char unsigned_byte;
  using ::xsde::cxx::hybrid::unsigned_byte_base;

  typedef short short_;
  using ::xsde::cxx::hybrid::short_base;

  typedef unsigned short unsigned_short;
  using ::xsde::cxx::hybrid::unsigned_short_base;

  typedef int int_;
  using ::xsde::cxx::hybrid::int_base;

  typedef unsigned int unsigned_int;
  using ::xsde::cxx::hybrid::unsigned_int_base;

  typedef long long long_;
  using ::xsde::cxx::hybrid::long_base;

  typedef unsigned long long unsigned_long;
  using ::xsde::cxx::hybrid::unsigned_long_base;

  typedef long integer;
  using ::xsde::cxx::hybrid::integer_base;

  typedef long non_positive_integer;
  using ::xsde::cxx::hybrid::non_positive_integer_base;

  typedef unsigned long non_negative_integer;
  using ::xsde::cxx::hybrid::non_negative_integer_base;

  typedef unsigned long positive_integer;
  using ::xsde::cxx::hybrid::positive_integer_base;

  typedef long negative_integer;
  using ::xsde::cxx::hybrid::negative_integer_base;

  typedef bool boolean;
  using ::xsde::cxx::hybrid::boolean_base;

  typedef float float_;
  using ::xsde::cxx::hybrid::float_base;

  typedef double double_;
  using ::xsde::cxx::hybrid::double_base;

  typedef double decimal;
  using ::xsde::cxx::hybrid::decimal_base;

  typedef ::std::string string;

  typedef ::std::string normalized_string;

  typedef ::std::string token;

  typedef ::std::string name;

  typedef ::std::string nmtoken;

  typedef ::xsde::cxx::string_sequence nmtokens;

  typedef ::std::string ncname;

  typedef ::std::string language;

  typedef ::std::string id;

  typedef ::std::string idref;

  typedef ::xsde::cxx::string_sequence idrefs;

  typedef ::std::string uri;

  using ::xsde::cxx::qname;

  using ::xsde::cxx::buffer;
  typedef ::xsde::cxx::buffer base64_binary;
  typedef ::xsde::cxx::buffer hex_binary;

  using ::xsde::cxx::time_zone;
  using ::xsde::cxx::date;
  using ::xsde::cxx::date_time;
  using ::xsde::cxx::duration;
  using ::xsde::cxx::gday;
  using ::xsde::cxx::gmonth;
  using ::xsde::cxx::gmonth_day;
  using ::xsde::cxx::gyear;
  using ::xsde::cxx::gyear_month;
  using ::xsde::cxx::time;

  using ::xsde::cxx::hybrid::pod_sequence;
  using ::xsde::cxx::hybrid::fix_sequence;
  using ::xsde::cxx::hybrid::var_sequence;
  using ::xsde::cxx::string_sequence;
  using ::xsde::cxx::hybrid::data_sequence;
}

namespace phenotype_dataset
{
  class phenotype_dataset;
}


// Begin epilogue.
//
//
// End epilogue.

#include <xsde/cxx/post.hxx>

#endif // PHENOTYPE_DATASET_FWD_HPP
