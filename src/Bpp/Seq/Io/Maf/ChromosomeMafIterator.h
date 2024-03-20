// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _CHROMOSOMEMAFITERATOR_H_
#define _CHROMOSOMEMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <set>

namespace bpp
{
/**
 * @brief Filter maf blocks to keep only blocks corresponding to a selection of chromosomes (of a reference sequence).
 */
class ChromosomeMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::string ref_;
  std::set<std::string> chr_;

public:
  /**
   * @param iterator The input iterator.
   * @param reference The reference species name.
   * @param chr the set of chromosome names to filter.
   */
  ChromosomeMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      const std::string& reference,
      const std::set<std::string>& chr) :
    AbstractFilterMafIterator(iterator),
    ref_(reference),
    chr_(chr)
  {}

  /**
   * @param iterator The input iterator.
   * @param reference The reference species name.
   * @param chr the chromosome name to filter.
   */
  ChromosomeMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      const std::string& reference,
      const std::string& chr) :
    AbstractFilterMafIterator(iterator),
    ref_(reference),
    chr_()
  {
    chr_.insert(chr);
  }

private:
  ChromosomeMafIterator(const ChromosomeMafIterator& iterator) :
    AbstractFilterMafIterator(0),
    ref_(iterator.ref_),
    chr_(iterator.chr_)
  {}

  ChromosomeMafIterator& operator=(const ChromosomeMafIterator& iterator)
  {
    ref_ = iterator.ref_;
    chr_ = iterator.chr_;
    return *this;
  }

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_();
};
} // end of namespace bpp.

#endif // _CHROMOSOMEMAFITERATOR_H_
