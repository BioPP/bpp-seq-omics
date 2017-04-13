//
// File: MafParser.h
// Authors: Julien Dutheil
// Created: Tue Apr 27 2010
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

#ifndef _MAFPARSER_H_
#define _MAFPARSER_H_

#include "MafIterator.h"
#include <Bpp/Seq/Alphabet/CaseMaskedAlphabet.h>

//From the STL:
#include <iostream>

namespace bpp {

/**
 * @brief MAF file parser.
 * 
 * This class parses synteny blocks from Maf file.
 *
 * The MAF format is documented on the UCSC Genome Browser website:
 * <a href="http://genome.ucsc.edu/FAQ/FAQformat.html#format5">http://genome.ucsc.edu/FAQ/FAQformat.html#format5</a>
 *
 * @author Julien Dutheil
 */
class MafParser:
  public AbstractMafIterator
{
  private:
    std::istream* stream_;
    bool mask_;
    CaseMaskedAlphabet cmAlphabet_;
    bool firstBlock_;
    short dotOption_;

  public:
    /**
     * @brief Create a new instance of MafParser
     *
     * @param stream The input stream to read text from
     * @param parseMask Tell is masking (lower case) should be kept
     * @param dotOption (one of DOT_ERROR, DOT_ASGAP or DOT_RESOLVE)
     *        tells how dot should be treated. DOT_ERROR, the default,
     *        will return an exception. DOT_ASGAP will convert all dots
     *        to gaps, which will increase parsing time.
     */
    MafParser(std::istream* stream, bool parseMask = false, short dotOption = DOT_ERROR) :
      stream_(stream), mask_(parseMask), cmAlphabet_(&AlphabetTools::DNA_ALPHABET), firstBlock_(true), dotOption_(dotOption) {}

  private:
    //Recopy is forbidden!
    MafParser(const MafParser& maf):
      stream_(0), mask_(maf.mask_), cmAlphabet_(&AlphabetTools::DNA_ALPHABET),
      firstBlock_(maf.firstBlock_), dotOption_(maf.dotOption_) {}

    MafParser& operator=(const MafParser& maf) {
      stream_ = 0;
      mask_ = maf.mask_;
      firstBlock_ = maf.firstBlock_;
      dotOption_ = maf.dotOption_;
      return *this;
    }

  private:
    MafBlock* analyseCurrentBlock_() throw (Exception);

  public:
    static constexpr short DOT_ERROR = 0;
    static constexpr short DOT_ASGAP = 1;
    //static constexpr short DOT_RESOLVE = 2; // not yet supported

};

} // end of namespace bpp.

#endif //_MAFPARSER_H_

