// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _ITERATIONLISTENER_H_
#define _ITERATIONLISTENER_H_

#include "MafBlock.h"

namespace bpp
{
/**
 * @brief Listener which enables to catch events when parsing a Maf file.
 */
class IterationListenerInterface
{
public:
  virtual ~IterationListenerInterface() {}

public:
  virtual void iterationStarts() = 0;
  virtual void iterationMoves(const MafBlock& currentBlock) = 0;
  virtual void iterationStops() = 0;
};
} // end of namespace bpp.

#endif // _ITERATIONLISTENER_H_
