// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
