//
// File: PlinkOutputMafIterator.h
// Authors: Julien Dutheil
// Created: Tue Dec 10 2015
//

/*
   Copyright or © or Copr. Bio++ Development Team, (2015)

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

#ifndef _PLINKOUTPUTMAFITERATOR_H_
#define _PLINKOUTPUTMAFITERATOR_H_

#include "MafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief This iterator outputs all biallelic SNPs in the PLINK format (ped and map files).
 */
class PlinkOutputMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::ostream* outputPed_;
  std::ostream* outputMap_;
  std::vector<std::string> species_;
  std::string refSpecies_;
  bool map3_;
  std::vector<std::string> ped_;
  std::string currentChr_;
  size_t lastPosition_;
  bool recodeChr_;
  std::map<std::string, unsigned int> chrCodes_;
  unsigned int currentCode_;
  int phenotype_;
  std::string colSeparator_; //Stored as a string to facilitate concatenation.

public:
  /**
   * @brief Build a new PlinkOutputMafIterator object.
   *
   * @warning In the current implementation, there is no way to deal with diploid individuals.
   * Each sequence is therefore considered as a haploid genome and will be written has a homozygous SNP
   * in the Ped file.
   *
   * @param iterator The input iterator.
   * @param outPed The output stream where to write the Ped file.
   * @param outMap The output stream where to write the Map file.
   * @param species A list of at least two species to compute SNPs.
   * Only blocks containing at least these two species will be used.
   * In case one species is duplicated in a block, the first sequence will be used.
   * @param reference The species to use as a reference for coordinates.
   * It does not have to be one of the selected species on which SNPs are computed.
   * @param map3 Tell if genetic distance column should be ommited in the map file. Otherwise set to 0.
   * @param recodeChr Tell if chromosomes should be recoded to numbers.
   * @param phenotype Phenotype value to set in the map file.
   * @param columnSeparator Character used to separate columns (PLINK officially supports space and tab).
   */
  PlinkOutputMafIterator(MafIterator* iterator,
                         std::ostream* outPed,
                         std::ostream* outMap,
                         const std::vector<std::string>& species,
                         const std::string& reference,
                         bool map3 = false,
                         bool recodeChr = false,
			 int phenotype = 0,
			 char columnSeparator = '\t') :
    AbstractFilterMafIterator(iterator),
    outputPed_(outPed), outputMap_(outMap), species_(species), refSpecies_(reference), map3_(map3),
    ped_(species.size()), currentChr_(""), lastPosition_(0), recodeChr_(recodeChr), chrCodes_(), currentCode_(1), phenotype_(0), colSeparator_(TextTools::toString(columnSeparator))
  {
    init_();
  }

private:
  PlinkOutputMafIterator(const PlinkOutputMafIterator& iterator) :
    AbstractFilterMafIterator(0),
    outputPed_(iterator.outputPed_),
    outputMap_(iterator.outputMap_),
    species_(iterator.species_),
    refSpecies_(iterator.refSpecies_),
    map3_(iterator.map3_),
    ped_(iterator.ped_),
    currentChr_(iterator.currentChr_),
    lastPosition_(iterator.lastPosition_),
    recodeChr_(iterator.recodeChr_),
    chrCodes_(iterator.chrCodes_),
    currentCode_(iterator.currentCode_),
    phenotype_(iterator.phenotype_),
    colSeparator_(iterator.colSeparator_)
  {}

  PlinkOutputMafIterator& operator=(const PlinkOutputMafIterator& iterator)
  {
    outputPed_       = iterator.outputPed_;
    outputMap_       = iterator.outputMap_;
    species_         = iterator.species_;
    refSpecies_      = iterator.refSpecies_;
    map3_            = iterator.map3_;
    ped_             = iterator.ped_;
    currentChr_      = iterator.currentChr_;
    lastPosition_    = iterator.lastPosition_;
    recodeChr_       = iterator.recodeChr_;
    chrCodes_        = iterator.chrCodes_;
    currentCode_     = iterator.currentCode_;
    phenotype_       = iterator.phenotype_;
    colSeparator_    = iterator.colSeparator_;
    return *this;
  }

public:
  MafBlock* analyseCurrentBlock_()
  {
    currentBlock_ = iterator_->nextBlock();
    if (outputMap_ && currentBlock_)
      parseBlock_(*outputMap_, *currentBlock_);
    if (outputMap_ && outputPed_ && !currentBlock_)
      writePedToFile_(*outputPed_); // Note we currently can output Map and no Ped, but not Ped without Map.
    return currentBlock_;
  }

private:
  void init_();
  void parseBlock_(std::ostream& out, const MafBlock& block);
  void writePedToFile_(std::ostream& out);
};
} // end of namespace bpp.

#endif//_PLINKOUTPUTMAFITERATOR_H_
