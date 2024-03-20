// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "ChromosomeRenamingMafIterator.h"

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>

using namespace std;

unique_ptr<MafBlock> ChromosomeRenamingMafIterator::analyseCurrentBlock_()
{
  currentBlock_ = iterator_->nextBlock();
  if (currentBlock_)
  {
    for (size_t i = 0; i < currentBlock_->getNumberOfSequences(); ++i)
    {
      string chr = currentBlock_->sequence(i).getChromosome();
      auto tln = chrTranslation_.find(chr);
      if (tln != chrTranslation_.end())
      {
        // We force conversion to avoid unecessary recopy
        const_cast<MafSequence&>(currentBlock_->sequence(i)).setChromosome(tln->second);
        if (logstream_)
        {
          (*logstream_ << "CHROMOSOME RENAMING: renamed " << chr << " to " << tln->second << ".").endLine();
        }
      }
    }
  }

  return move(currentBlock_);
}
