// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "AbstractIterationListener.h"

// From the STL:
#include <vector>

using namespace std;
using namespace bpp;

void CsvStatisticsOutputIterationListener::iterationStarts()
{
  const vector<string>& header = statsIterator_->getResultsColumnNames();
  *output_ << "Chr" << sep_ << "Start" << sep_ << "Stop";
  for (size_t i = 0; i < header.size(); ++i)
  {
    *output_ << sep_ << header[i];
  }
  output_->endLine();
}

void CsvStatisticsOutputIterationListener::iterationMoves(const MafBlock& currentBlock)
{
  auto& values = statsIterator_->getResults();
  if (currentBlock.hasSequenceForSpecies(refSpecies_))
  {
    const auto& refSeq = currentBlock.sequenceForSpecies(refSpecies_);
    if (refSeq.hasCoordinates())
      *output_ << refSeq.getChromosome() << sep_ << refSeq.start() << sep_ << refSeq.stop();
    else
      *output_ << "NA" << sep_ << "NA" << sep_ << "NA";
  }
  else
  {
    *output_ << "NA" << sep_ << "NA" << sep_ << "NA";
  }
  for (size_t i = 0; i < values.size(); ++i)
  {
    *output_ << sep_ << (values[i] ? values[i]->toString() : "NA");
  }
  output_->endLine();
}
