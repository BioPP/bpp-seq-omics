//
// File: BlockMergerMafIterator.h
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

#ifndef _BLOCKMERGERMAFITERATOR_H_
#define _BLOCKMERGERMAFITERATOR_H_

#include "MafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief Merge blocks if some of their sequences are contiguous.
 *
 * The user specifies the focus species. Sequences that are not in this set will
 * be automatically merged and their coordinates removed.
 * The scores, if any, will be averaged for the block, weighted by the corresponding block sizes.
 * The pass value will be removed if it is different for the two blocks.
 * It is possible to define a maximum distance for the merging. Setting a distance of zero implies that the blocks
 * have to be exactly contiguous. Alternatively, the appropriate number of 'N' will be inserted in all species.
 * All species however have to be distant of the exact same amount.
 */
class BlockMergerMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::vector<std::string> species_;
  MafBlock* incomingBlock_;
  std::vector<std::string> ignoreChrs_; // These chromosomes will never be merged (ex: 'Un').
  unsigned int maxDist_;
  bool renameChimericChromosomes_;
  std::map<std::string, unsigned int> chimericChromosomeCounts_;

public:
  BlockMergerMafIterator(MafIterator* iterator, const std::vector<std::string>& species, unsigned int maxDist = 0, bool renameChimericChromosomes = false) :
    AbstractFilterMafIterator(iterator),
    species_(species),
    incomingBlock_(0),
    ignoreChrs_(),
    maxDist_(maxDist),
    renameChimericChromosomes_(renameChimericChromosomes),
    chimericChromosomeCounts_()
  {
    incomingBlock_ = iterator->nextBlock();
  }

private:
  BlockMergerMafIterator(const BlockMergerMafIterator& iterator) :
    AbstractFilterMafIterator(0),
    species_(iterator.species_),
    incomingBlock_(iterator.incomingBlock_),
    ignoreChrs_(iterator.ignoreChrs_),
    maxDist_(iterator.maxDist_),
    renameChimericChromosomes_(iterator.renameChimericChromosomes_),
    chimericChromosomeCounts_(iterator.chimericChromosomeCounts_)
  {}

  BlockMergerMafIterator& operator=(const BlockMergerMafIterator& iterator)
  {
    species_                   = iterator.species_;
    incomingBlock_             = iterator.incomingBlock_;
    ignoreChrs_                = iterator.ignoreChrs_;
    maxDist_                   = iterator.maxDist_;
    renameChimericChromosomes_ = iterator.renameChimericChromosomes_;
    chimericChromosomeCounts_  = iterator.chimericChromosomeCounts_;
    return *this;
  }

public:
  /**
   * brief Add a chromosome that should be ignored to the list.
   * @param chr The name of the chromosome to be ignored.
   */
  void ignoreChromosome(const std::string& chr)
  {
    ignoreChrs_.push_back(chr);
  }

private:
  MafBlock* analyseCurrentBlock_();
};
} // end of namespace bpp.

#endif//_BLOCKMERGERMAFITERATOR_H_
