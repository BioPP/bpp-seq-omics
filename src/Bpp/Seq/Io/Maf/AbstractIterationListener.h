// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _ABSTRACTITERATIONLISTENER_H_
#define _ABSTRACTITERATIONLISTENER_H_

#include "MafIterator.h"
#include "SequenceStatisticsMafIterator.h"

namespace bpp
{
/**
 * @brief Iteration listener that works with a SequenceStatisticsMafIterator,
 * enabling output of results in a file (partial implementation, format-independent)
 */
class AbstractStatisticsOutputIterationListener :
  public virtual IterationListenerInterface
{
protected:
  std::shared_ptr<SequenceStatisticsMafIterator> statsIterator_;

public:
  AbstractStatisticsOutputIterationListener(
      std::shared_ptr<SequenceStatisticsMafIterator> iterator) :
    statsIterator_(iterator) {}

  AbstractStatisticsOutputIterationListener(const AbstractStatisticsOutputIterationListener& listener) :
    statsIterator_(listener.statsIterator_) {}

  AbstractStatisticsOutputIterationListener& operator=(const AbstractStatisticsOutputIterationListener& listener)
  {
    statsIterator_ = listener.statsIterator_;
    return *this;
  }

  virtual ~AbstractStatisticsOutputIterationListener() {}
};

/**
 * @brief Iteration listener that works with a SequenceStatisticsMafIterator,
 * enabling output of results in a file in CSV format
 */
class CsvStatisticsOutputIterationListener :
  public AbstractStatisticsOutputIterationListener
{
private:
  std::shared_ptr<OutputStream> output_;
  std::string sep_;
  std::string refSpecies_;

public:
  CsvStatisticsOutputIterationListener(
      std::shared_ptr<SequenceStatisticsMafIterator> iterator,
      const std::string& refSpecies,
      std::shared_ptr<OutputStream> output,
      const std::string& sep = "\t") :
    AbstractStatisticsOutputIterationListener(iterator), output_(output), sep_(sep), refSpecies_(refSpecies) {}

  CsvStatisticsOutputIterationListener(const CsvStatisticsOutputIterationListener& listener) :
    AbstractStatisticsOutputIterationListener(listener), output_(listener.output_), sep_(listener.sep_), refSpecies_(listener.refSpecies_) {}

  CsvStatisticsOutputIterationListener& operator=(const CsvStatisticsOutputIterationListener& listener)
  {
    AbstractStatisticsOutputIterationListener::operator=(listener);
    output_ = listener.output_;
    sep_ = listener.sep_;
    refSpecies_ = listener.refSpecies_;
    return *this;
  }

  virtual ~CsvStatisticsOutputIterationListener() {}

public:
  virtual void iterationStarts();
  virtual void iterationMoves(const MafBlock& currentBlock);
  virtual void iterationStops() {}
};
} // end of namespace bpp.

#endif // _ABSTRACTITERATIONLISTENER_H_
