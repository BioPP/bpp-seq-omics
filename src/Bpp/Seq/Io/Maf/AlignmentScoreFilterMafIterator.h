// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _ALIGNMENTSCOREFILTERMAFITERATOR_H_
#define _ALIGNMENTSCOREFILTERMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>

namespace bpp
{
/**
 * @brief Discard block with low alignment scores.
 */
class AlignmentScoreFilterMafIterator :
  public AbstractFilterMafIterator
{
private:
  double minScore_;

public:
  /**
   * @brief Creates a new AlignmentScoreFilterMafIterator object.
   *
   * @param iterator The input iterator.
   * @param minimumScore Minimum score.
   */
  AlignmentScoreFilterMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      double minimumScore) :
    AbstractFilterMafIterator(iterator),
    minScore_(minimumScore)
  {}

private:
  AlignmentScoreFilterMafIterator(const AlignmentScoreFilterMafIterator& iterator) = delete;

  AlignmentScoreFilterMafIterator& operator=(const AlignmentScoreFilterMafIterator& iterator) = delete;

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_() override
  {
    bool test;
    do
    {
      currentBlock_ = iterator_->nextBlock();
      if (!currentBlock_) break;
      test = (currentBlock_->getScore() < minScore_);
      if (test)
      {
        if (logstream_)
        {
          (*logstream_ << "ALIGNMENT SCORE FILTER: block " << currentBlock_->getDescription() << " with score " << currentBlock_->getScore() << " was discarded.").endLine();
        }
        currentBlock_ = 0;
      }
    }
    while (test);
    return std::move(currentBlock_);
  }
};
} // end of namespace bpp.

#endif // _ALIGNMENTSCOREFILTERMAFITERATOR_H_
