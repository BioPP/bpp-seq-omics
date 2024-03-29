// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "GffFeatureReader.h"

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

const std::string GffFeatureReader::GFF_PHASE = "GFF_PHASE";
const std::string GffFeatureReader::GFF_NAME = "Name";
const std::string GffFeatureReader::GFF_ALIAS = "GFF_ALIAS";
const std::string GffFeatureReader::GFF_PARENT = "Parent";
const std::string GffFeatureReader::GFF_TARGET = "Target";
const std::string GffFeatureReader::GFF_GAP = "Gap";
const std::string GffFeatureReader::GFF_DERIVES_FROM = "GFF_DERIVES_FROM";
const std::string GffFeatureReader::GFF_NOTE = "Note";
const std::string GffFeatureReader::GFF_DBXREF = "Dbxref";
const std::string GffFeatureReader::GFF_ONTOLOGY_TERM = "Ontology_term";
const std::string GffFeatureReader::GFF_IS_CIRCULAR = "Is_circular";


void GffFeatureReader::getNextLine_()
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

const BasicSequenceFeature GffFeatureReader::nextFeature()
{
  if (!hasMoreFeature())
    throw Exception("GffFeatureReader::nextFeature(). No more feature in file.");

  // Parse current line:
  StringTokenizer st(nextLine_, "\t");
  if (st.numberOfRemainingTokens() != 9)
    throw Exception("GffFeatureReader::nextFeature(). Wrong GFF3 file format: should have 9 tab delimited columns.");

  // if ok, we can parse each column:
  string seqId       = st.nextToken();
  string source      = st.nextToken();
  string type        = st.nextToken();
  unsigned int start = TextTools::to<unsigned int>(st.nextToken()) - 1;
  unsigned int end   = TextTools::to<unsigned int>(st.nextToken());
  double score       = TextTools::to<double>(st.nextToken());
  string strand      = st.nextToken();
  string phase       = st.nextToken();
  string attrDesc    = st.nextToken();
  map<string, string> attributes;
  KeyvalTools::multipleKeyvals(attrDesc, attributes, ";", false);
  string id = attributes["ID"];
  BasicSequenceFeature feature(id, seqId, source, type, start, end, strand[0], score);

  // Set phase attributes:
  if (phase != ".")
    feature.setAttribute(GFF_PHASE, phase);

  // now check additional attributes:
  for (map<string, string>::iterator it = attributes.begin(); it != attributes.end(); ++it)
  {
    if (it->first != "ID")
      feature.setAttribute(it->first, it->second); // We accept all attributes, even if they are not standard.
  }

  // Read the next line:
  getNextLine_();

  return feature;
}

std::string GffFeatureReader::toString(const bpp::SequenceFeature& f)
{
  std::vector< std::string > v;
  std::vector< std::string > attr;
  std::set< std::string > attrNames = f.getAttributeList();
  v.push_back(f.getSequenceId());
  v.push_back(f.getSource());
  v.push_back(f.getType());
  v.push_back(bpp::TextTools::toString(f.getStart() + 1));
  v.push_back(bpp::TextTools::toString(f.getEnd()));
  v.push_back(bpp::TextTools::toString(f.getScore()));
  if (f.isStranded())
  {
    if (f.isNegativeStrand())
    {
      v.push_back("-");
    }
    else
    {
      v.push_back("+");
    }
  }
  else
  {
    v.push_back(".");
  }
  if (f.getAttribute(GFF_PHASE) == "")
  {
    v.push_back(".");
  }
  else
  {
    v.push_back(f.getAttribute(GFF_PHASE));
  }

  if (f.getId() != "")
  {
    attr.push_back("ID=" + f.getId());
  }
  for (std::set< std::string >::iterator it = attrNames.begin(); it != attrNames.end(); it++)
  {
    attr.push_back(*it + "=" + f.getAttribute(*it));
  }
  v.push_back(bpp::VectorTools::paste(attr, ";"));
  return bpp::VectorTools::paste(v, "\t");
}
