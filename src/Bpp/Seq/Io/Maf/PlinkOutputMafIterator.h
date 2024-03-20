// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _PLINKOUTPUTMAFITERATOR_H_
#define _PLINKOUTPUTMAFITERATOR_H_

#include "AbstractMafIterator.h"

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
  std::shared_ptr<std::ostream> outputPed_;
  std::shared_ptr<std::ostream> outputMap_;
  std::vector<std::string> species_;
  std::string refSpecies_;
  bool map3_;
  std::vector<std::string> ped_;
  std::string currentChr_;
  size_t lastPosition_;
  bool recodeChr_;
  std::map<std::string, unsigned int> chrCodes_;
  unsigned int currentCode_;
  bool makeDiploids_;
  int phenotype_;
  std::string colSeparator_; // Stored as a string to facilitate concatenation.
  size_t nbIndividuals_;

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
   * @param makeDiploids If true, combines genomes into diploids. In case of odd numbers, the last genome will be ignored. If false, each genome will be output as a homozygous diploid.
   * @param phenotype Phenotype value to set in the map file.
   * @param columnSeparator Character used to separate columns (PLINK officially supports space and tab).
   */
  PlinkOutputMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      std::shared_ptr<std::ostream> outPed,
      std::shared_ptr<std::ostream> outMap,
      const std::vector<std::string>& species,
      const std::string& reference,
      bool map3 = false,
      bool recodeChr = false,
      bool makeDiploids = false,
      int phenotype = 0,
      char columnSeparator = '\t') :
    AbstractFilterMafIterator(iterator),
    outputPed_(outPed),
    outputMap_(outMap),
    species_(species),
    refSpecies_(reference),
    map3_(map3),
    ped_(),
    currentChr_(""),
    lastPosition_(0),
    recodeChr_(recodeChr),
    chrCodes_(),
    currentCode_(1),
    makeDiploids_(makeDiploids),
    phenotype_(phenotype),
    colSeparator_(TextTools::toString(columnSeparator)),
    nbIndividuals_(0)
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
    makeDiploids_(iterator.makeDiploids_),
    phenotype_(iterator.phenotype_),
    colSeparator_(iterator.colSeparator_),
    nbIndividuals_(iterator.nbIndividuals_)
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
    makeDiploids_    = iterator.makeDiploids_;
    phenotype_       = iterator.phenotype_;
    colSeparator_    = iterator.colSeparator_;
    nbIndividuals_   = iterator.nbIndividuals_;
    return *this;
  }

public:
  std::unique_ptr<MafBlock> analyseCurrentBlock_()
  {
    currentBlock_ = iterator_->nextBlock();
    if (outputMap_ && currentBlock_)
      parseBlock_(*outputMap_, *currentBlock_);
    if (outputMap_ && outputPed_ && !currentBlock_)
      writePedToFile_(*outputPed_); // Note we currently can output Map and no Ped, but not Ped without Map.
    return std::move(currentBlock_);
  }

private:
  void init_();
  void parseBlock_(std::ostream& out, const MafBlock& block);
  void writePedToFile_(std::ostream& out);
};
} // end of namespace bpp.

#endif // _PLINKOUTPUTMAFITERATOR_H_
