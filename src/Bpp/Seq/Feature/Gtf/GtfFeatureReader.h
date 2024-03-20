// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _GTFFEATUREREADER_H_
#define _GTFFEATUREREADER_H_

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
 * @brief A simple reader implementing the Gene Transfer Format.
 *
 * The reference norm in use is the one of GTF2.2 http://mblab.wustl.edu/GTF22.html .
 * This class is a "beta" class, and may undeavour interface changes in the future.
 *
 * Note that in GTF, coordinates are [a, b] 1-based. They will therefore be converted to [a, b[ 0-based,
 * as specified for the SequenceFeature object.
 *
 * @author Sylvain Gaillard
 */
class GtfFeatureReader :
  public virtual FeatureReader
{
public:
  static const std::string GTF_PHASE;
  static const std::string GTF_GENE_ID;
  static const std::string GTF_TRANSCRIPT_ID;

private:
  std::istream& input_;
  std::string nextLine_;

public:
  GtfFeatureReader(std::istream& input) :
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

private:
  void getNextLine_();
};
} // end of namespace bpp

#endif//_GTFFEATUREREADER_H_
