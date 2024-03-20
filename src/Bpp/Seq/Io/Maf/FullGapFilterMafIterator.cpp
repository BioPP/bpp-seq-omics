// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "FullGapFilterMafIterator.h"

// From bpp-seq
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/SiteTools.h>

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>

using namespace std;

unique_ptr<MafBlock> FullGapFilterMafIterator::analyseCurrentBlock_()
{
  auto block = iterator_->nextBlock();
  if (!block)
    return nullptr;

  // We create a copy of the ingroup alignement for better efficiency:
  VectorSiteContainer vsc(AlphabetTools::DNA_ALPHABET);
  for (size_t i = 0; i < species_.size(); ++i)
  {
    if (block->hasSequenceForSpecies(species_[i]))
    {
      auto tmpSeq = make_unique<Sequence>(block->sequenceForSpecies(species_[i]));
      vsc.addSequence(tmpSeq->getName(), tmpSeq);
    }
  }
  if (vsc.getNumberOfSequences() == 0)
    return block; // Block ignored as it does not contain any of the focus species.

  // Now check the positions that are only made of gaps:
  if (verbose_)
  {
    ApplicationTools::message->endLine();
    ApplicationTools::displayTask("Cleaning block for gap sites", true);
  }
  size_t n = block->getNumberOfSites();
  vector<size_t> start;
  vector<unsigned int> count;
  bool test = false;
  for (size_t i = 0; i < n; ++i)
  {
    const Site& site = vsc.site(i);
    if (SiteTools::isGapOnly(site))
    {
      if (test)
      {
        count[count.size() - 1]++;
      }
      else
      {
        start.push_back(i);
        count.push_back(1);
        test = true;
      }
    }
    else
    {
      test = false;
    }
  }
  // Now remove blocks:
  size_t totalRemoved = 0;
  for (size_t i = start.size(); i > 0; --i)
  {
    if (verbose_)
      ApplicationTools::displayGauge(start.size() - i, start.size() - 1, '=');
    block->deleteSites(start[i - 1], count[i - 1]);
    totalRemoved += count[i - 1];
  }
  if (verbose_)
    ApplicationTools::displayTaskDone();

  // Correct coordinates:
  if (totalRemoved > 0)
  {
    for (size_t i = 0; i < block->getNumberOfSequences(); ++i)
    {
      if (!VectorTools::contains(species_, block->sequence(i).getSpecies()))
      {
        block->removeCoordinatesFromSequence(i);
      }
    }
  }
  if (logstream_)
  {
    (*logstream_ << "FULL GAP CLEANER: " << totalRemoved << " positions have been removed.").endLine();
  }
  return block;
}
