// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _CHROMOSOMERENAMINGMAFITERATOR_H_
#define _CHROMOSOMERENAMINGMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <map>

namespace bpp
{
/**
 * @brief Rename chromosomes according to a translation table.
 */
class ChromosomeRenamingMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::map<std::string, std::string> chrTranslation_;

public:
  /**
   * @param iterator The input iterator.
   * @param chrTranslation a map with original chromosome names as keys and translations as values. Only chromosomes matching one of the key will be translated, without further checking of the new name.
   */
  ChromosomeRenamingMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      const std::map<std::string, std::string>& chrTranslation) :
    AbstractFilterMafIterator(iterator),
    chrTranslation_(chrTranslation)
  {}

private:
  ChromosomeRenamingMafIterator(const ChromosomeRenamingMafIterator& iterator) :
    AbstractFilterMafIterator(0),
    chrTranslation_(iterator.chrTranslation_)
  {}

  ChromosomeRenamingMafIterator& operator=(const ChromosomeRenamingMafIterator& iterator)
  {
    chrTranslation_ = iterator.chrTranslation_;
    return *this;
  }

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_();
};
} // end of namespace bpp.

#endif // _CHROMOSOMERENAMINGMAFITERATOR_H_
