// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _VCFOUTPUTMAFITERATOR_H_
#define _VCFOUTPUTMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief This iterator performs a simple SNP call from the MAF blocks, and outputs the results in the Variant Call Format (VCF).
 *
 * Only SNPs are supported for now.
 */
class VcfOutputMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::shared_ptr<std::ostream> output_;
  std::string refSpecies_;
  std::vector<std::vector<std::string>> genotypes_;
  bool outputAll_;
  bool generateDiploids_;

public:
  /**
   * @brief Build a new VcfOutputMafIterator object.
   *
   * @param iterator The input iterator.
   * @param out The output stream where to write the VCF file.
   * @param reference The species to use as a reference.
   * @param genotypes A list of species for which genotype information should be written in the VCF file. There will be one extra column per genotype, +1 format column.
   * @param outputAll If true, also output non-variable positions.
   * @param generateDiploids If true, output artificial "homozygous" diploids.
   */
  VcfOutputMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      std::shared_ptr<std::ostream> out,
      const std::string& reference,
      const std::vector<std::string>& genotypes,
      bool outputAll = false,
      bool generateDiploids = false) :
    AbstractFilterMafIterator(iterator),
    output_(out),
    refSpecies_(reference),
    genotypes_(),
    outputAll_(outputAll),
    generateDiploids_(generateDiploids)
  {
    for (auto g : genotypes)
    {
      std::vector<std::string> tmp;
      tmp.push_back(g);
      genotypes_.push_back(tmp);
    }
    if (output_)
      writeHeader_(*output_);
  }

  /**
   * @brief Build a new VcfOutputMafIterator object.
   *
   * @param iterator The input iterator.
   * @param out The output stream where to write the VCF file.
   * @param reference The species to use as a reference.
   * @param genotypes A list of species combinations for which genotype information should be written in the VCF file. There will be one extra column per genotype, +1 format column. When more than one sequence is specified in a combination, a (phased) polyploid genotype will be created.
   * @param outputAll If true, also output non-variable positions.
   */
  VcfOutputMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      std::shared_ptr<std::ostream> out,
      const std::string& reference,
      const std::vector< std::vector<std::string> >& genotypes,
      bool outputAll = false) :
    AbstractFilterMafIterator(iterator),
    output_(out),
    refSpecies_(reference),
    genotypes_(genotypes),
    outputAll_(outputAll),
    generateDiploids_(false)
  {
    if (output_)
      writeHeader_(*output_);
  }

private:
  VcfOutputMafIterator(const VcfOutputMafIterator& iterator) :
    AbstractFilterMafIterator(0),
    output_(iterator.output_),
    refSpecies_(iterator.refSpecies_),
    genotypes_(iterator.genotypes_),
    outputAll_(iterator.outputAll_),
    generateDiploids_(iterator.generateDiploids_)
  {}

  VcfOutputMafIterator& operator=(const VcfOutputMafIterator& iterator)
  {
    output_ = iterator.output_;
    refSpecies_ = iterator.refSpecies_;
    genotypes_ = iterator.genotypes_;
    outputAll_ = iterator.outputAll_;
    generateDiploids_ = iterator.generateDiploids_;
    return *this;
  }

public:
  std::unique_ptr<MafBlock> analyseCurrentBlock_()
  {
    currentBlock_ = iterator_->nextBlock();
    if (output_ && currentBlock_)
      writeBlock_(*output_, *currentBlock_);
    return move(currentBlock_);
  }

private:
  void writeHeader_(std::ostream& out) const;
  void writeBlock_(std::ostream& out, const MafBlock& block) const;
};
} // end of namespace bpp.

#endif//_VCFOUTPUTMAFITERATOR_H_
