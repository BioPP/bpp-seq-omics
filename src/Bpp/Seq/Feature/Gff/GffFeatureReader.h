// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _GFFFEATUREREADER_H_
#define _GFFFEATUREREADER_H_

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
 * @brief A simple reader implementing the Gene Finding Feature format.
 *
 * The reference norm in use is the one of GFF3 http://www.sequenceontology.org/gff3.shtml .
 * This class is a "beta" class, and may undeavour interface changes in the future.
 *
 * Note that in GFF, coordinates are [a, b] 1-based. They will therefore be converted to [a, b[ 0-based,
 * as specified for the SequenceFeature object.
 *
 * @author Julien Dutheil, Sylvain Gaillard
 */
class GffFeatureReader :
  public virtual FeatureReader
{
public:
  static const std::string GFF_STRAND;
  static const std::string GFF_PHASE;
  static const std::string GFF_NAME;
  static const std::string GFF_ALIAS;
  static const std::string GFF_PARENT;
  static const std::string GFF_TARGET;
  static const std::string GFF_GAP;
  static const std::string GFF_DERIVES_FROM;
  static const std::string GFF_NOTE;
  static const std::string GFF_DBXREF;
  static const std::string GFF_ONTOLOGY_TERM;
  static const std::string GFF_IS_CIRCULAR;

private:
  std::istream& input_;
  std::string nextLine_;

public:
  GffFeatureReader(std::istream& input) :
    input_(input), nextLine_()
  {
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

#endif // _GFFFEATUREREADER_H_
