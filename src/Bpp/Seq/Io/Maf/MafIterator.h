// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _MAFITERATOR_H_
#define _MAFITERATOR_H_

#include "MafBlock.h"
#include "IterationListener.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief Interface to loop over maf alignment blocks.
 */
class MafIteratorInterface
{
public:
  virtual ~MafIteratorInterface() {}

public:
  /**
   * @brief Get the next available alignment block.
   *
   * @return A maf alignment block, or a null pointer if no more block is available.
   */
  virtual std::unique_ptr<MafBlock> nextBlock() = 0;

  virtual bool isVerbose() const = 0;

  virtual void setVerbose(bool yn) = 0;

  virtual void addIterationListener(std::unique_ptr<IterationListenerInterface> listener) = 0;
};


/**
 * @brief Interface to loop over removed blocks of a maf alignment.
 */
class MafTrashIteratorInterface
{
public:
  virtual ~MafTrashIteratorInterface() {}

public:
  /**
   * @brief Get the next available removed alignment block.
   *
   * @return A maf alignment block, or a null pointer if no more block is available.
   */
  virtual std::unique_ptr<MafBlock> nextRemovedBlock() = 0;
};
} // end of namespace bpp.

#endif // _MAFITERATOR_H_
