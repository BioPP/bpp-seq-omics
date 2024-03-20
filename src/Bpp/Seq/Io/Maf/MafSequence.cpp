// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "MafSequence.h"

// From the STL:
#include <string>

using namespace bpp;
using namespace std;

unique_ptr<MafSequence> MafSequence::subSequence(size_t startAt, size_t length) const
{
  string subseq = toString().substr(startAt, length);
  size_t begin = begin_;
  if (hasCoordinates_)
  {
    for (size_t i = 0; i < startAt; ++i)
    {
      if (!getAlphabet()->isGap(operator[](i)))
        begin++;
    }
  }
  auto newSeq = make_unique<MafSequence>(getName(), subseq, begin, strand_, srcSize_);
  if (!hasCoordinates_)
    newSeq->removeCoordinates();
  vector<string> anno = getAnnotationTypes();
  for (size_t i = 0; i < anno.size(); ++i)
  {
    newSeq->addAnnotation(annotation(anno[i]).getPartAnnotation(startAt, length));
  }
  return newSeq;
}
