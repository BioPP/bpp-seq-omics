// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _WINDOWSPLITMAFITERATOR_H_
#define _WINDOWSPLITMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief Splits block into windows of given sizes.
 */
class WindowSplitMafIterator :
  public AbstractFilterMafIterator
{
private:
  size_t windowSize_;
  size_t windowStep_;
  short align_;
  std::deque<std::unique_ptr<MafBlock>> blockBuffer_;
  bool keepSmallBlocks_;

public:
  static const short RAGGED_LEFT;
  static const short RAGGED_RIGHT;
  static const short CENTER;
  static const short ADJUST;

public:
  WindowSplitMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      size_t windowSize,
      size_t windowStep,
      short splitOption = CENTER,
      bool keepSmallBlocks = false) :
    AbstractFilterMafIterator(iterator),
    windowSize_(windowSize), 
    windowStep_(windowStep),
    align_(splitOption),
    blockBuffer_(),
    keepSmallBlocks_(keepSmallBlocks)
  {
    if (splitOption != RAGGED_LEFT && splitOption != RAGGED_RIGHT
        && splitOption != CENTER && splitOption != ADJUST)
      throw Exception("WindowSplitMafIterator: unvalid split option: " + TextTools::toString(splitOption));
    if (splitOption != RAGGED_LEFT && windowStep != windowSize)
      throw Exception("WindowSplitMafIterator: overlapping windows are only supported together with the RAGGED_LEFT option.");
  }

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_();
};
} // end of namespace bpp.

#endif//_WINDOWSPLITMAFITERATOR_H_
