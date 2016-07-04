//
// File: AlignmentFilterMafIterator.h
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

#ifndef _ALIGNMENTFILTERMAFITERATOR_H_
#define _ALIGNMENTFILTERMAFITERATOR_H_

#include "MafIterator.h"

//From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp {

/**
 * @brief Filter maf blocks to remove ambiguously aligned or non-informative regions.
 *
 * Regions with a too high proportion of gaps, unknown character or high entropy in a set of species will be removed,
 * and blocks adjusted accordingly. 
 *
 * The total entropy of a window is defined as @f$ \frac{\sum_i {E(5)}_i}{w} @f$ where @f$E(5)_i@f$ is the shannon entropy of site i, computed using a logarithm of base 5, and w is the window size.
 * As a result, the total entropy is normalized so that it falls between 0 and 1. A logarithm of base 5 is used to account for the fact that gaps are considered as a state.
 *
 * In case a sequence from the list is missing, it can be either ignored or counted as a full sequence of gaps.
 */
class AlignmentFilterMafIterator:
  public AbstractFilterMafIterator,
  public virtual MafTrashIterator
{
  private:
    std::vector<std::string> species_;
    unsigned int windowSize_;
    unsigned int step_;
    unsigned int maxGap_;
    double maxPropGap_;
    double maxEnt_;
    std::deque<MafBlock*> blockBuffer_;
    std::deque<MafBlock*> trashBuffer_;
    std::deque< std::vector<int> > window_;
    bool keepTrashedBlocks_;
    bool missingAsGap_;
    bool relative_;

  public:

    /**
     * @brief Create a new AlignmentFilterMafIterator with absolute thresholds.
     *
     * @param iterator Input iterator
     * @param species Selection of species on which filtering criteria are applied.
     * Results of filtering will be applied to all species.
     * @param windowSize Size of the sliding window (nt).
     * @param step Step by which windows are moved (nt).
     * @param maxGap Maximum number of gaps allowed in the window.
     * @param maxEnt Maximum entropy allowed in the window.
     * @param keepTrashedBlocks Removed windows are kept as separate blocks.
     * @param missingAsGap Add missing species as gap sequences where needed.
     */
    AlignmentFilterMafIterator(MafIterator* iterator,
        const std::vector<std::string>& species,
        unsigned int windowSize,
        unsigned int step,
        unsigned int maxGap,
        double maxEnt,
        bool keepTrashedBlocks,
        bool missingAsGap) :
      AbstractFilterMafIterator(iterator),
      species_(species),
      windowSize_(windowSize),
      step_(step),
      maxGap_(maxGap),
      maxPropGap_(),
      maxEnt_(maxEnt),
      blockBuffer_(),
      trashBuffer_(),
      window_(species.size()),
      keepTrashedBlocks_(keepTrashedBlocks),
      missingAsGap_(missingAsGap),
      relative_(false)
    {}

    /**
     * @brief Create a new AlignmentFilterMafIterator with relative thresholds.
     *
     * @param iterator Input iterator
     * @param species Selection of species on which filtering criteria are applied.
     * Results of filtering will be applied to all species.
     * @param windowSize Size of the sliding window (nt).
     * @param step Step by which windows are moved (nt).
     * @param maxPropGap Maximum proportion of gaps allowed in the window.
     * @param maxEnt Maximum entropy allowed in the window.
     * @param keepTrashedBlocks Removed windows are kept as separate blocks.
     * @param missingAsGap Add missing species as gap sequences where needed.
     */
    AlignmentFilterMafIterator(MafIterator* iterator,
        const std::vector<std::string>& species,
        unsigned int windowSize,
        unsigned int step,
        double maxPropGap,
        double maxEnt,
        bool keepTrashedBlocks,
        bool missingAsGap) :
      AbstractFilterMafIterator(iterator),
      species_(species),
      windowSize_(windowSize),
      step_(step),
      maxGap_(),
      maxPropGap_(maxPropGap),
      maxEnt_(maxEnt),
      blockBuffer_(),
      trashBuffer_(),
      window_(species.size()),
      keepTrashedBlocks_(keepTrashedBlocks),
      missingAsGap_(missingAsGap),
      relative_(true)
    {}


  public:
    MafBlock* nextRemovedBlock() throw (Exception) {
      if (trashBuffer_.size() == 0) return 0;
      MafBlock* block = trashBuffer_.front();
      trashBuffer_.pop_front();
      return block;
    }

  private:
    MafBlock* analyseCurrentBlock_() throw (Exception);
};

/**
 * @brief Filter maf blocks to remove ambiguously aligned or non-informative regions.
 *
 * This iterators offers a different algorithm than AlignmentFilterMafIterator.
 * It takes two parameters: g=maxGap and n=maxPos. Windows with more than n positions containing each of them more than g=maxPos gaps will be discarded.
 * In addition, consecutives patterns are only counted once.
 * In case a sequence from the list is missing, it can be either ignored or counted as a full sequence of gaps.
 */
class AlignmentFilter2MafIterator:
  public AbstractFilterMafIterator,
  public virtual MafTrashIterator
{
  private:
    std::vector<std::string> species_;
    unsigned int windowSize_;
    unsigned int step_;
    unsigned int maxGap_;
    double maxPropGap_;
    unsigned int maxPos_;
    std::deque<MafBlock*> blockBuffer_;
    std::deque<MafBlock*> trashBuffer_;
    std::deque< std::vector<bool> > window_;
    bool keepTrashedBlocks_;
    bool missingAsGap_;
    bool relative_;

  public:
    /**
     * @brief Create a new AlignmentFilter2MafIterator with absolute thresholds.
     *
     * @param iterator Input iterator
     * @param species Selection of species on which filtering criteria are applied.
     * Results of filtering will be applied to all species.
     * @param windowSize Size of the sliding window (nt).
     * @param step Step by which windows are moved (nt).
     * @param maxGap Maximum number of gaps allowed in the window.
     * @param maxPos Maximum number of gaps "events" allowed.
     * @param keepTrashedBlocks Removed windows are kept as separate blocks.
     * @param missingAsGap Add missing species as gap sequences where needed.
     */
    AlignmentFilter2MafIterator(MafIterator* iterator, const std::vector<std::string>& species, unsigned int windowSize, unsigned int step, unsigned int maxGap, unsigned int maxPos, bool keepTrashedBlocks, bool missingAsGap) :
      AbstractFilterMafIterator(iterator),
      species_(species),
      windowSize_(windowSize),
      step_(step),
      maxGap_(maxGap),
      maxPropGap_(),
      maxPos_(maxPos),
      blockBuffer_(),
      trashBuffer_(),
      window_(species.size()),
      keepTrashedBlocks_(keepTrashedBlocks),
      missingAsGap_(missingAsGap),
      relative_(false)
    {}

    /**
     * @brief Create a new AlignmentFilterMafIterator with relative thresholds.
     *
     * @param iterator Input iterator
     * @param species Selection of species on which filtering criteria are applied.
     * Results of filtering will be applied to all species.
     * @param windowSize Size of the sliding window (nt).
     * @param step Step by which windows are moved (nt).
     * @param maxPropGap Maximum proportion of gaps allowed in the window.
     * @param maxGap Maximum number of gaps allowed in the window.
     * @param keepTrashedBlocks Removed windows are kept as separate blocks.
     * @param missingAsGap Add missing species as gap sequences where needed.
     */
    AlignmentFilter2MafIterator(MafIterator* iterator, const std::vector<std::string>& species, unsigned int windowSize, unsigned int step, double maxPropGap, unsigned int maxPos, bool keepTrashedBlocks, bool missingAsGap) :
      AbstractFilterMafIterator(iterator),
      species_(species),
      windowSize_(windowSize),
      step_(step),
      maxGap_(),
      maxPropGap_(maxPropGap),
      maxPos_(maxPos),
      blockBuffer_(),
      trashBuffer_(),
      window_(species.size()),
      keepTrashedBlocks_(keepTrashedBlocks),
      missingAsGap_(missingAsGap),
      relative_(true)
    {}

  public:
    MafBlock* nextRemovedBlock() throw (Exception) {
      if (trashBuffer_.size() == 0) return 0;
      MafBlock* block = trashBuffer_.front();
      trashBuffer_.pop_front();
      return block;
    }

  private:
    MafBlock* analyseCurrentBlock_() throw (Exception);

};

} // end of namespace bpp.

#endif //_ALIGNMENTFILTERMAFITERATOR_H_
