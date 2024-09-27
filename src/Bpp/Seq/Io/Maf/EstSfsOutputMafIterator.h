// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _ESTSFSOUTPUTMAFITERATOR_H_
#define _ESTSFSOUTPUTMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief Output data in the format read by Peter D. Keightley's EST-SFS
 *
 * One to three outgroups are supported. EST-SFS estimate ancestral alleles and compute unfolded site frequency spectra (uSFS).
 *
 * @see Keightley and Jackson, Genetics 209: 897-906 (2018).
 */
class EstSfsOutputMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::shared_ptr<std::ostream> output_;
  std::vector<std::string> ingroup_;
  std::vector<std::string> outgroup1_;
  std::vector<std::string> outgroup2_;
  std::vector<std::string> outgroup3_;

public:
  /**
   * @brief Build a new EstSfsOutputMafIterator object.
   *
   * Output data in the Input format of the EST-SFS program from Peter D. Keightley.
   *
   * @param iterator The input iterator.
   * @param out The output stream where to write the EST-SFS file.
   * @param ingroup The list of genome ids to use as a part of the focal species.
   * @param outgroup1 The list of genome ids to use as the first outgroup.
   * @param outgroup2 The list of genome ids to use as the second outgroup.
   * @param outgroup3 The list of genome ids to use as the third outgroup.
   */
  EstSfsOutputMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      std::shared_ptr<std::ostream> out,
      const std::vector<std::string>& ingroup,
      const std::vector<std::string>& outgroup1,
      const std::vector<std::string>& outgroup2,
      const std::vector<std::string>& outgroup3) :
    AbstractFilterMafIterator(iterator),
    output_(out),
    ingroup_(ingroup),
    outgroup1_(outgroup1),
    outgroup2_(outgroup2),
    outgroup3_(outgroup3)
  {
    if (outgroup1.size() == 0) {
      throw Exception("EstSfsOutputMafIterator::constructor. At least one outgroup species is required.");
    }
  }

  EstSfsOutputMafIterator(const EstSfsOutputMafIterator& iterator) = delete;
    
  EstSfsOutputMafIterator& operator=(const EstSfsOutputMafIterator& iterator) = delete;
 
public:
  std::unique_ptr<MafBlock> analyseCurrentBlock_()
  {
    currentBlock_ = iterator_->nextBlock();
    if (output_ && currentBlock_)
      writeBlock_(*output_, *currentBlock_);
    return std::move(currentBlock_);
  }

private:
  void writeBlock_(std::ostream& out, const MafBlock& block) const;
};
} // end of namespace bpp.

#endif // _ESTSFSOUTPUTMAFITERATOR_H_
