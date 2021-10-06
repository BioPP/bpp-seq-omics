//
// File: BedGraphFeatureReader.cpp
// Created by: Julien Dutheil
// Created on: Tue Feb 2 2016
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

   This software is a computer program whose purpose is to provide classes
   for sequences analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#include "BedGraphFeatureReader.h"

// From bpp-core:
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Numeric/VectorTools.h>

// From the STL:
#include <string>
#include <iostream>

using namespace bpp;
using namespace std;

const std::string BedGraphFeatureReader::BED_VALUE = "BedValue";


void BedGraphFeatureReader::getNextLine_()
{
  nextLine_ = "";
  while (TextTools::isEmpty(nextLine_) || nextLine_.size() < 2 || nextLine_[0] == '#')
  {
    if (input_.eof())
    {
      nextLine_ = "";
      return;
    }
    getline(input_, nextLine_);
  }
}

const BasicSequenceFeature BedGraphFeatureReader::nextFeature()
{
  if (!hasMoreFeature())
    throw Exception("BedGraphFeatureReader::nextFeature(). No more feature in file.");

  // Parse current line:
  StringTokenizer st(nextLine_, "\t");
  if (st.numberOfRemainingTokens() != 4)
    throw Exception("BedGraphFeatureReader::nextFeature(). Wrong BedGraph file format: should have 4 tab delimited columns.");

  // if ok, we can parse each column:
  string seqId       = st.nextToken();
  unsigned int start = TextTools::to<unsigned int>(st.nextToken());
  unsigned int end   = TextTools::to<unsigned int>(st.nextToken());
  string value       = st.nextToken();
  string id          = "bed" + TextTools::toString(++id_);
  BasicSequenceFeature feature(id, seqId, "bed_graph", "", start, end, '.', -1);

  // Set value attributes:
  if (value != ".")
    feature.setAttribute(BED_VALUE, value);

  // Read the next line:
  getNextLine_();

  return feature;
}

std::string BedGraphFeatureReader::toString(const bpp::SequenceFeature& f)
{
  std::vector< std::string > v;
  v.push_back(f.getSequenceId());
  v.push_back(bpp::TextTools::toString(f.getStart()));
  v.push_back(bpp::TextTools::toString(f.getEnd()));
  string value = f.getAttribute(BED_VALUE);
  v.push_back(bpp::TextTools::toString(value == "" ? "." : value));
  return bpp::VectorTools::paste(v, "\t");
}
