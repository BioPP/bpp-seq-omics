// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _SEQUENCESTREAMTOMAFITERATOR_H_
#define _SEQUENCESTREAMTOMAFITERATOR_H_

#include "AbstractMafIterator.h"
#include <Bpp/Seq/Alphabet/CaseMaskedAlphabet.h>
#include <Bpp/Seq/Io/ISequenceStream.h>

// From the STL:
#include <iostream>
#include <memory>

namespace bpp
{
/**
 * @brief A MafIterator built from a sequence stream.
 *
 * Each block will contain one sequence from the original file.
 *
 * @author Julien Dutheil
 */
class SequenceStreamToMafIterator :
  public AbstractMafIterator
{
private:
  std::shared_ptr<ISequenceStream> seqStream_;
  std::shared_ptr<std::istream> stream_;
  bool zeroBasedCoords_;
  bool firstBlock_;

public:
  SequenceStreamToMafIterator(
      std::shared_ptr<ISequenceStream> seqStream,
      std::shared_ptr<std::istream> stream,
      bool parseMask = false,
      bool zeroBasedCoordinates = true) :
    seqStream_(seqStream),
    stream_(stream),
    zeroBasedCoords_(zeroBasedCoordinates),
    firstBlock_(true)
  {}

private:
  // Recopy is forbidden!
  SequenceStreamToMafIterator(const SequenceStreamToMafIterator& ss2mi) :
    seqStream_(), 
    stream_(nullptr),
    zeroBasedCoords_(ss2mi.zeroBasedCoords_),
    firstBlock_(ss2mi.firstBlock_)
  {}
  
  SequenceStreamToMafIterator& operator=(const SequenceStreamToMafIterator& ss2mi)
  {
    seqStream_.reset();
    stream_ = 0; 
    zeroBasedCoords_ = ss2mi.zeroBasedCoords_; 
    firstBlock_ = ss2mi.firstBlock_;
    return *this;
  }

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_();
};
} // end of namespace bpp.

#endif//_SEQUENCESTREAMTOMAFITERATOR_H_
