//
// File: ChromosomeRenamingMafIterator.cpp
// Authors: Julien Dutheil
// Created: Thu Feb 03 2022
//

/*
   Copyright or © or Copr. Bio++ Development Team, (2022)

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
      if (tln != chrTranslation_.end()) {
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

