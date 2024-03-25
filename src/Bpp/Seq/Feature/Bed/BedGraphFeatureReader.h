// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _BEDGRAPHFEATUREREADER_H_
#define _BEDGRAPHFEATUREREADER_H_

#include "../SequenceFeature.h"
#include "../FeatureReader.h"

// From bpp-core:
#include <Bpp/Exceptions.h>

// From the STL:
#include <string>
#include <vector>

namespace bpp
{
/**
 * @brief A simple reader for features in the BedGraph format.
 *
 * Format description at UCSC: https://genome.ucsc.edu/goldenpath/help/bedgraph.html
 *
 * Note: The value associated to each feature is stored as a string attribute, using
 * tag BegGraphFeatureReader::BED_VALUE. No check is performed regarding its value.
 * An automatic id is generated, and the source tag is set to "beg_graph".
 *
 * @author Julien Dutheil
 */
class BedGraphFeatureReader :
  public virtual FeatureReader
{
public:
  static const std::string BED_VALUE;

private:
  std::istream& input_;
  std::string nextLine_;
  unsigned int id_;

public:
  BedGraphFeatureReader(std::istream& input) :
    input_(input), nextLine_(), id_(0)
  {
    bool start = false;
    do
    {
      getNextLine_();
      if (nextLine_.size() >= 5 && nextLine_.substr(0, 5) == "track")
      {
        start = true;
      }
    }
    while (!start && !input_.eof());
    if (input_.eof())
      throw Exception("BedGraphFeatureReader::constructor: Invalid BedGraph file, missing proper header.");
    getNextLine_();
  }

public:
  bool hasMoreFeature() const { return nextLine_ != ""; }
  const BasicSequenceFeature nextFeature();

  void getAllFeatures(SequenceFeatureSet& features)
  {
    while (hasMoreFeature())
    {
      features.addFeature(nextFeature());
    }
  }
  void getFeaturesOfType(const std::string& type, SequenceFeatureSet& features)
  {
    while (hasMoreFeature())
    {
      BasicSequenceFeature feature = nextFeature();
      if (feature.getType() == type)
        features.addFeature(feature);
    }
  }
  void getFeaturesOfSequence(const std::string& seqId, SequenceFeatureSet& features)
  {
    while (hasMoreFeature())
    {
      BasicSequenceFeature feature = nextFeature();
      if (feature.getSequenceId() == seqId)
        features.addFeature(feature);
    }
  }

  /**
   * @param f A sequence feature.
   * @return A string describing the feature, in GFF format.
   */
  static std::string toString(const bpp::SequenceFeature& f);

  /**
   * @brief Out put a string description of a feature to a stream.
   *
   * A end of line character will be appended after the description.
   *
   * @param f A sequence feature.
   * @param out An output stream.
   */
  static void toString(const bpp::SequenceFeature& f, std::ostream& out)
  {
    out << toString(f) << std::endl;
  }

  static void toString(const bpp::SequenceFeatureSet& fs, std::ostream& out)
  {
    for (unsigned int i = 0; i < fs.getNumberOfFeatures(); ++i)
    {
      toString(fs[i], out);
    }
  }

private:
  void getNextLine_();
};
} // end of namespace bpp

#endif // _BEDGRAPHFEATUREREADER_H_
