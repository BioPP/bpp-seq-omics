// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "GtfFeatureReader.h"

// From bpp-core:
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>

// From the STL:
#include <string>
#include <iostream>

using namespace bpp;
using namespace std;

const string GtfFeatureReader::GTF_PHASE = "GTF_PHASE";
const string GtfFeatureReader::GTF_GENE_ID = "gene_id";
const string GtfFeatureReader::GTF_TRANSCRIPT_ID = "transcript_id";

void GtfFeatureReader::getNextLine_()
{
  nextLine_ = "";
  // Remove comments
  /*
     todo: remove characters between # and end-of-line taking care of double quotes
   */
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

const BasicSequenceFeature GtfFeatureReader::nextFeature()
{
  if (!hasMoreFeature())
    throw Exception("GtfFeatureReader::nextFeature(). No more feature in file.");

  // Parse current line:
  StringTokenizer st(nextLine_, "\t");
  if (st.numberOfRemainingTokens() != 9)
    throw Exception("GtfFeatureReader::nextFeature(). Wrong GTF file format: should have 9 tab delimited columns.");

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
  StringTokenizer st1(attrDesc, ";");
  while (st1.hasMoreToken())
  {
    string item(st1.nextToken());
    if (TextTools::isEmpty(item))
      continue;
    item = TextTools::removeSurroundingWhiteSpaces(item);
    auto idx = item.find_first_of(' ');
    if (idx != string::npos)
    {
      string key(item.substr(0, idx));
      string value(item.substr(idx));
      // remove first "
      while (
        value.size() > 0
        && (
          value[0] == '"'
          || TextTools::isWhiteSpaceCharacter(value[0])
          )
        )
      {
        value.erase(value.begin());
      }
      // remove last "
      while (
        value.size() > 0
        && (
          value[value.size() - 1] == '"'
          || TextTools::isWhiteSpaceCharacter(value[value.size() - 1])
          )
        )
      {
        value.erase(value.end() - 1);
      }
      attributes[key] = value;
    // std::cout << "[" << key << "] = [" << value << "]" << std::endl;
    } else {
      attributes[item] = ""; //Empty value.
    }
  }
  // KeyvalTools::multipleKeyvals(attrDesc, attributes, ";", false);
  // std::string id = attributes["ID"];
  string id = "";
  BasicSequenceFeature feature(id, seqId, source, type, start, end, strand[0], score);

  // Set phase attributes:
  phase = TextTools::removeSurroundingWhiteSpaces(phase);
  if (phase != ".")
    feature.setAttribute(GTF_PHASE, phase);

  // now check additional attributes:
  for (map<std::string, std::string>::iterator it = attributes.begin(); it != attributes.end(); ++it)
  {
    feature.setAttribute(it->first, it->second); // We accept all attributes, even if they are not standard.
    // std::cout << "[" << it->first << "] = [" << it->second << "]" << std::endl;
    // std::cout << "phase: " << phase << std::endl;
  }

  // Read the next line:
  getNextLine_();

  return feature;
}
