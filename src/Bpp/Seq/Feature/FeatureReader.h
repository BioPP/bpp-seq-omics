// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _FEATUREREADER_H_
#define _FEATUREREADER_H_

#include "SequenceFeature.h"

// From bpp-core:
#include <Bpp/Exceptions.h>

// From the STL:
#include <string>
#include <vector>

namespace bpp
{
/**
 * @brief Interface for feature readers.
 *
 * @author Julien Dutheil, Sylvain Gaillard
 */
class FeatureReader
{
public:
  FeatureReader() {}
  virtual ~FeatureReader() {}

public:
  virtual bool hasMoreFeature() const = 0;
  virtual const BasicSequenceFeature nextFeature() = 0;
  virtual void getAllFeatures(SequenceFeatureSet& features) = 0;
  virtual void getFeaturesOfType(const std::string& type, SequenceFeatureSet& features) = 0;
  virtual void getFeaturesOfSequence(const std::string& seqId, SequenceFeatureSet& features) = 0;
};
} // end of namespace bpp

#endif // _FEATUREREADER_H_
