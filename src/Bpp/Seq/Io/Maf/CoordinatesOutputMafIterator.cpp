// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "CoordinatesOutputMafIterator.h"

using namespace bpp;
using namespace std;

void CoordinatesOutputMafIterator::writeHeader_(ostream& out) const
{
  for (size_t i = 0; i < species_.size(); ++i)
  {
    if (i > 0)
      out << "\t";
    string sp = species_[i];
    out << sp << ".chr\t" << sp << ".strand\t" << sp << ".start\t" << sp << ".stop";
    if (includeSrcSize_)
      out << "\t" << sp << ".src";
  }
  out << endl;
}

std::unique_ptr<MafBlock> CoordinatesOutputMafIterator::analyseCurrentBlock_()
{
  currentBlock_ = iterator_->nextBlock();
  if (currentBlock_)
  {
    for (size_t i = 0; i < species_.size(); ++i)
    {
      if (i > 0)
        *output_ << "\t";
      vector<const MafSequence*> seqs = currentBlock_->getSequencesForSpecies(species_[i]);
      if (seqs.size() > 1)
        throw Exception("CoordinatesOutputMafIterator::analyseCurrentBlock_(). There is more than one sequence for species '" + species_[i] + "' in current block.");
      else if (seqs.size() == 0)
      {
        *output_ << "NA\tNA\tNA\tNA";
        if (includeSrcSize_)
          *output_ << "\tNA";
      }
      else
      {
        *output_ << seqs[0]->getChromosome() << "\t" << seqs[0]->getStrand() << "\t" << seqs[0]->start() << "\t" << seqs[0]->stop();
        if (includeSrcSize_)
          *output_ << "\t" << seqs[0]->getSrcSize();
      }
    }
    *output_ << endl;
  }
  return move(currentBlock_);
}
