// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _MAFPARSER_H_
#define _MAFPARSER_H_

#include "AbstractMafIterator.h"
#include <Bpp/Seq/Alphabet/CaseMaskedAlphabet.h>

// From the STL:
#include <iostream>

namespace bpp
{
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
class MafParser :
  public AbstractMafIterator
{
private:
  std::shared_ptr<std::istream> stream_;
  bool mask_;
  bool checkSequenceSize_;
  CaseMaskedAlphabet cmAlphabet_;
  bool firstBlock_;
  short dotOption_;

public:
  /**
   * @brief Create a new instance of MafParser
   *
   * @param stream The input stream to read text from
   * @param parseMask Tell is masking (lower case) should be kept
   * @param checkSize Tell if the size of sequence found should be
   *        compared to the specified one. An exception is thrown
   *        in case of mismatch (default). If set to no, a warning
   *        will be displayed if verbose is set to true.
   * @param dotOption (one of DOT_ERROR, DOT_ASGAP or DOT_ASUNRES)
   *        tells how dot should be treated. DOT_ERROR, the default,
   *        will return an exception. DOT_ASGAP will convert all dots
   *        to gaps and DOT_ASUNRES will convert them to 'N', which
   *        will increase parsing time.
   */
  MafParser(
      std::shared_ptr<std::istream> stream, 
      bool parseMask = false,
      bool checkSize = true,
      short dotOption = DOT_ERROR) :
    stream_(stream),
    mask_(parseMask),
    checkSequenceSize_(checkSize),
    cmAlphabet_(AlphabetTools::DNA_ALPHABET),
    firstBlock_(true),
    dotOption_(dotOption)
  {}

private:
  // Recopy is forbidden!
  MafParser(const MafParser& maf) :
    stream_(nullptr), mask_(maf.mask_), checkSequenceSize_(maf.checkSequenceSize_),
    cmAlphabet_(AlphabetTools::DNA_ALPHABET), firstBlock_(maf.firstBlock_),
    dotOption_(maf.dotOption_) {}

  MafParser& operator=(const MafParser& maf)
  {
    stream_ = nullptr;
    mask_ = maf.mask_;
    checkSequenceSize_ = maf.checkSequenceSize_;
    firstBlock_ = maf.firstBlock_;
    dotOption_ = maf.dotOption_;
    return *this;
  }

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_();

public:
  static constexpr short DOT_ERROR = 0;
  static constexpr short DOT_ASGAP = 1;
  static constexpr short DOT_ASUNRES = 2;
  // static constexpr short DOT_RESOLVE = 3; // not yet supported
};
} // end of namespace bpp.

#endif//_MAFPARSER_H_
