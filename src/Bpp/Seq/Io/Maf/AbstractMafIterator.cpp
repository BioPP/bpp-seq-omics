// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "AbstractMafIterator.h"

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>

using namespace std;

void AbstractMafIterator::fireIterationStartSignal_()
{
  for (auto& it : iterationListeners_)
  {
    it->iterationStarts();
  }
}

void AbstractMafIterator::fireIterationMoveSignal_(const MafBlock& currentBlock)
{
  for (auto& it : iterationListeners_)
  {
    it->iterationMoves(currentBlock);
  }
}

void AbstractMafIterator::fireIterationStopSignal_()
{
  for (auto& it : iterationListeners_)
  {
    it->iterationStops();
  }
}
