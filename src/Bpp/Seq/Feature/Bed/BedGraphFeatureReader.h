//
// File: BedGraphFeatureReader.h
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
 * Format desciption at UCSC: https://genome.ucsc.edu/goldenpath/help/bedgraph.html
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
    while(!start && !input_.eof());
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

#endif//_BEDGRAPHFEATUREREADER_H_
