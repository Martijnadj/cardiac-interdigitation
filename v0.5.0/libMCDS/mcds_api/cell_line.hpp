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

#ifndef CELL_LINE_HPP
#define CELL_LINE_HPP

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

#include "cell_line-fwd.hpp"

#ifndef XSDE_DONT_INCLUDE_INLINE
#define XSDE_DONT_INCLUDE_INLINE

#include "metadata.hpp"

#include "phenotype_dataset.hpp"

#include "common.hpp"

#undef XSDE_DONT_INCLUDE_INLINE
#else

#include "metadata.hpp"

#include "phenotype_dataset.hpp"

#include "common.hpp"

#endif // XSDE_DONT_INCLUDE_INLINE

namespace cell_line
{
  // cell_line (variable-length)
  //
  class cell_line
  {
    private:
    cell_line (const cell_line&);
    cell_line& operator= (const cell_line&);

    public:
    cell_line ();

    cell_line*
    _clone () const;

    ~cell_line ();

    // ID
    //
    bool
    ID_present () const;

    void
    ID_present (bool);

    const ::std::string&
    ID () const;

    ::std::string&
    ID ();

    void
    ID (const ::std::string&);

    // label
    //
    bool
    label_present () const;

    void
    label_present (bool);

    const ::std::string&
    label () const;

    ::std::string&
    label ();

    void
    label (const ::std::string&);

    // curated
    //
    bool
    curated_present () const;

    void
    curated_present (bool);

    bool
    curated () const;

    bool&
    curated ();

    void
    curated (bool);

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

    // phenotype_dataset
    //
    typedef ::xsde::cxx::hybrid::var_sequence< ::phenotype_dataset::phenotype_dataset > phenotype_dataset_sequence;
    typedef phenotype_dataset_sequence::iterator phenotype_dataset_iterator;
    typedef phenotype_dataset_sequence::const_iterator phenotype_dataset_const_iterator;

    const phenotype_dataset_sequence&
    phenotype_dataset () const;

    phenotype_dataset_sequence&
    phenotype_dataset ();

    // custom
    //
    bool
    custom_present () const;

    const ::common::custom&
    custom () const;

    ::common::custom&
    custom ();

    void
    custom (::common::custom*);

    ::common::custom*
    custom_detach ();

    void
    _copy (cell_line&) const;

    private:
    ::std::string ID_;
    unsigned char ID_present_;
    ::std::string label_;
    unsigned char label_present_;
    bool curated_;
    unsigned char curated_present_;
    ::metadata::metadata* metadata_;
    phenotype_dataset_sequence phenotype_dataset_;
    ::common::custom* custom_;
  };

  // DCLs (variable-length)
  //
  class DCLs
  {
    private:
    DCLs (const DCLs&);
    DCLs& operator= (const DCLs&);

    public:
    DCLs ();

    DCLs*
    _clone () const;

    ~DCLs ();

    // cell_line
    //
    typedef ::xsde::cxx::hybrid::var_sequence< ::cell_line::cell_line > cell_line_sequence;
    typedef cell_line_sequence::iterator cell_line_iterator;
    typedef cell_line_sequence::const_iterator cell_line_const_iterator;

    const cell_line_sequence&
    cell_line () const;

    cell_line_sequence&
    cell_line ();

    void
    _copy (DCLs&) const;

    private:
    cell_line_sequence cell_line_;
  };
}

#ifndef XSDE_DONT_INCLUDE_INLINE

#include "metadata.ipp"

#include "phenotype_dataset.ipp"

#include "common.ipp"

#endif // XSDE_DONT_INCLUDE_INLINE

#ifndef XSDE_DONT_INCLUDE_INLINE
#include "cell_line.ipp"
#endif // XSDE_DONT_INCLUDE_INLINE

// Begin epilogue.
//
//
// End epilogue.

#include <xsde/cxx/post.hxx>

#endif // CELL_LINE_HPP
