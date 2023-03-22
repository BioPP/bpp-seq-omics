//
// File: FullGapFilterMafIterator.cpp
// Authors: Julien Dutheil
// Created: Tue Sep 07 2010
//

/*
   Copyright or © or Copr. Bio++ Development Team, (2010)

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
