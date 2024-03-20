// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _SEQUENCESTATISTICSMAFITERATOR_H_
#define _SEQUENCESTATISTICSMAFITERATOR_H_

#include "AbstractMafIterator.h"
#include "MafStatistics.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief Compute a series of sequence statistics for each block.
 *
 * Computed statistics are stored into a vector of double, which can be retrieved as well as statistics names.
 * Listeners can be set up to automatically analyse or write the output after iterations are over.
 *
 * The current implementation focuses on speed and memory efificiency, as it only stores in memory the current results of the statistics.
 * The only drawback of this, is that disk access might be high when writing the results,
 * although appropriate buffering should most likely circumvent the issue.
 * The code is easily extensible, however, to enable storage of all results into a matrix,
 * with writing only once at the end of iterations.
 */
class SequenceStatisticsMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::vector<std::shared_ptr<MafStatisticsInterface>> statistics_;
  std::vector<std::unique_ptr<BppNumberI>> results_;
  std::vector<std::string> names_;

public:
  /**
   * @param iterator The input iterator.
   * @param statistics A vector of pointers toward MafStatistics, to be computed simultaneously for each maf block.
   */
  SequenceStatisticsMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      const std::vector<std::shared_ptr<MafStatisticsInterface>> statistics);

private:
  SequenceStatisticsMafIterator(const SequenceStatisticsMafIterator& iterator) :
    AbstractFilterMafIterator(0),
    statistics_(iterator.statistics_),
    results_(),
    names_(iterator.names_)
  {}

  SequenceStatisticsMafIterator& operator=(const SequenceStatisticsMafIterator& iterator)
  {
    statistics_ = iterator.statistics_;
    results_.clear();
    names_ = iterator.names_;
    return *this;
  }

public:
  const std::vector<std::unique_ptr<BppNumberI>>& getResults() const { return results_; }
  const std::vector<std::string>& getResultsColumnNames() const { return names_; }

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_();
};
} // end of namespace bpp.

#endif // _SEQUENCESTATISTICSMAFITERATOR_H_
