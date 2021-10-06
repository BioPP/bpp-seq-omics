//
// File: SequenceStreamToMafIterator.h
// Authors: Julien Dutheil
// Created: Wed Oct 31 2012
//

/*
   Copyright or © or Copr. Bio++ Development Team, (2012)

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

#ifndef _SEQUENCESTREAMTOMAFITERATOR_H_
#define _SEQUENCESTREAMTOMAFITERATOR_H_

#include "MafIterator.h"
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
  std::unique_ptr<ISequenceStream> seqStream_;
  std::istream* stream_;
  bool zeroBasedCoords_;
  bool firstBlock_;

public:
  SequenceStreamToMafIterator(ISequenceStream* seqStream, std::istream* stream, bool parseMask = false, bool zeroBasedCoordinates = true) :
    seqStream_(seqStream), stream_(stream), zeroBasedCoords_(zeroBasedCoordinates), firstBlock_(true) {}

private:
  // Recopy is forbidden!
  SequenceStreamToMafIterator(const SequenceStreamToMafIterator& ss2mi) :
    seqStream_(), stream_(0), zeroBasedCoords_(ss2mi.zeroBasedCoords_), firstBlock_(ss2mi.firstBlock_) {}
  SequenceStreamToMafIterator& operator=(const SequenceStreamToMafIterator& ss2mi)
  {
    seqStream_.reset(); stream_ = 0; zeroBasedCoords_ = ss2mi.zeroBasedCoords_; firstBlock_ = ss2mi.firstBlock_;
    return *this;
  }

private:
  MafBlock* analyseCurrentBlock_();
};
} // end of namespace bpp.

#endif//_SEQUENCESTREAMTOMAFITERATOR_H_
