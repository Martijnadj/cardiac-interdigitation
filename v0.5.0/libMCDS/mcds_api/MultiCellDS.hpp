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

#ifndef MULTI_CELL_DS_HPP
#define MULTI_CELL_DS_HPP

#include <xsde/cxx/version.hxx>

#if (XSDE_INT_VERSION != 3020000L)
#error XSD/e runtime version mismatch
#endif

#include <xsde/cxx/config.hxx>

#ifndef XSDE_ENCODING_UTF8
#error the generated code uses the UTF-8 encodingwhile the XSD/e runtime does not (reconfigure the runtime or change the --char-encoding value)
#endif

#ifndef XSDE_STL
#error the generated code uses STL while the XSD/e runtime does not (reconfigure the runtime or add --no-stl)
#endif

#ifndef XSDE_EXCEPTIONS
#error the generated code uses exceptions while the XSD/e runtime does not (reconfigure the runtime or add --no-exceptions)
#endif

#ifndef XSDE_LONGLONG
#error the generated code uses long long while the XSD/e runtime does not (reconfigure the runtime or add --no-long-long)
#endif

#ifdef XSDE_CUSTOM_ALLOCATOR
#error the XSD/e runtime uses custom allocator while the generated code does not (reconfigure the runtime or add --custom-allocator)
#endif

#include <xsde/cxx/pre.hxx>

// Begin prologue.
//
//
// End prologue.

#include "MultiCellDS-fwd.hpp"

#ifndef XSDE_DONT_INCLUDE_INLINE
#define XSDE_DONT_INCLUDE_INLINE

#include "common.hpp"

#include "metadata.hpp"

#include "cell_line.hpp"

#include "cell.hpp"

#include "microenvironment.hpp"

#undef XSDE_DONT_INCLUDE_INLINE
#else

#include "common.hpp"

#include "metadata.hpp"

#include "cell_line.hpp"

#include "cell.hpp"

#include "microenvironment.hpp"

#endif // XSDE_DONT_INCLUDE_INLINE

// MCDS_type (fixed-length)
//
class MCDS_type
{
  public:
  enum value_type
  {
    cell_line,
    snapshot_simulation,
    snapshot_experiment,
    snapshot_clinical
  };

  MCDS_type ();
  MCDS_type (value_type);

  void
  value (value_type);

  operator value_type () const
  {
    return value_;
  }

  const char*
  string () const;

  private:
  value_type value_;
};

// MultiCellDS (variable-length)
//
class MultiCellDS
{
  private:
  MultiCellDS (const MultiCellDS&);
  MultiCellDS& operator= (const MultiCellDS&);

  public:
  MultiCellDS ();

  MultiCellDS*
  _clone () const;

  ~MultiCellDS ();

  // version
  //
  bool
  version_present () const;

  void
  version_present (bool);

  const ::std::string&
  version () const;

  ::std::string&
  version ();

  void
  version (const ::std::string&);

  // type
  //
  bool
  type_present () const;

  void
  type_present (bool);

  const ::MCDS_type&
  type () const;

  ::MCDS_type&
  type ();

  void
  type (const ::MCDS_type&);

  // cell_line
  //
  typedef ::xsde::cxx::hybrid::var_sequence< ::cell_line::cell_line > cell_line_sequence;
  typedef cell_line_sequence::iterator cell_line_iterator;
  typedef cell_line_sequence::const_iterator cell_line_const_iterator;

  const cell_line_sequence&
  cell_line () const;

  cell_line_sequence&
  cell_line ();

  // metadata
  //
  bool
  metadata_present () const;

  const ::metadata::metadata&
  metadata () const;

  ::metadata::metadata&
  metadata ();

  void
  metadata (::metadata::metadata*);

  ::metadata::metadata*
  metadata_detach ();

  // microenvironment
  //
  bool
  microenvironment_present () const;

  const ::microenvironment::microenvironment&
  microenvironment () const;

  ::microenvironment::microenvironment&
  microenvironment ();

  void
  microenvironment (::microenvironment::microenvironment*);

  ::microenvironment::microenvironment*
  microenvironment_detach ();

  // cellular_information
  //
  bool
  cellular_information_present () const;

  const ::cell::cellular_information&
  cellular_information () const;

  ::cell::cellular_information&
  cellular_information ();

  void
  cellular_information (::cell::cellular_information*);

  ::cell::cellular_information*
  cellular_information_detach ();

  void
  _copy (MultiCellDS&) const;

  private:
  ::std::string version_;
  unsigned char version_present_;
  ::MCDS_type type_;
  unsigned char type_present_;
  cell_line_sequence cell_line_;
  ::metadata::metadata* metadata_;
  ::microenvironment::microenvironment* microenvironment_;
  ::cell::cellular_information* cellular_information_;
};

#ifndef XSDE_DONT_INCLUDE_INLINE

#include "common.ipp"

#include "metadata.ipp"

#include "cell_line.ipp"

#include "cell.ipp"

#include "microenvironment.ipp"

#endif // XSDE_DONT_INCLUDE_INLINE

#ifndef XSDE_DONT_INCLUDE_INLINE
#include "MultiCellDS.ipp"
#endif // XSDE_DONT_INCLUDE_INLINE

// Begin epilogue.
//
//
// End epilogue.

#include <xsde/cxx/post.hxx>

#endif // MULTI_CELL_DS_HPP
