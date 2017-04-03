//
// File: VcfOutputMafIterator.h
// Authors: Julien Dutheil
// Created: Tue Jan 05 2013
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

#ifndef _VCFOUTPUTMAFITERATOR_H_
#define _VCFOUTPUTMAFITERATOR_H_

#include "MafIterator.h"

//From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp {

/**
 * @brief This iterator performs a simple SNP call from the MAF blocks, and outputs the results in the Variant Call Format (VCF).
 *
 * Only substitutions are supported for now.
 */
class VcfOutputMafIterator:
  public AbstractFilterMafIterator
{
  private:
    std::ostream* output_;
    std::string refSpecies_;
    std::vector<std::string> genotypes_;
    bool outputAll_;

  public:
    /**
     * @brief Build a new VcfOutputMafIterator object.
     *
     * @param iterator The input iterator.
     * @param out The output stream where to write the VCF file.
     * @param reference The species to use as a reference.
     * @param genotypes A list of species for which genotype information should be written in the VCF file. There will be one extra column per genotype, +1 format column.
     * @param outputAll If true, also output non-variable positions.
     */
    VcfOutputMafIterator(MafIterator* iterator, std::ostream* out, const std::string& reference, const std::vector<std::string>& genotypes, bool outputAll = false) :
      AbstractFilterMafIterator(iterator), output_(out), refSpecies_(reference), genotypes_(genotypes), outputAll_(outputAll)
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
      outputAll_(iterator.outputAll_)
    {}
    
    VcfOutputMafIterator& operator=(const VcfOutputMafIterator& iterator)
    {
      output_ = iterator.output_;
      refSpecies_ = iterator.refSpecies_;
      genotypes_ = iterator.genotypes_;
      outputAll_ = iterator.outputAll_;
      return *this;
    }


  public:
    MafBlock* analyseCurrentBlock_() throw (Exception) {
      currentBlock_ = iterator_->nextBlock();
      if (output_ && currentBlock_)
        writeBlock_(*output_, *currentBlock_);
      return currentBlock_;
    }

  private:
    void writeHeader_(std::ostream& out) const;
    void writeBlock_(std::ostream& out, const MafBlock& block) const;
};

} // end of namespace bpp.

#endif //_VCFOUTPUTMAFITERATOR_H_
