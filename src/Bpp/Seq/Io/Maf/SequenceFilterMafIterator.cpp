//
// File: SequenceFilterMafIterator.cpp
// Authors: Julien Dutheil
// Created: Tue Sep 07 2010
//

/*
   Copyright or © or Copr. Bio++ Development Team, (2010)

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

#include "SequenceFilterMafIterator.h"

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>

using namespace std;

MafBlock* SequenceFilterMafIterator::analyseCurrentBlock_()
{
  currentBlock_ = iterator_->nextBlock();
  while (currentBlock_)
  {
    map<string, unsigned int> counts;
    for (size_t i = currentBlock_->getNumberOfSequences(); i > 0; --i)
    {
      string species = currentBlock_->getMafSequence(i - 1).getSpecies();
      if (!VectorTools::contains(species_, species))
      {
        if (logstream_)
        {
          (*logstream_ << "SEQUENCE FILTER: remove sequence '" << species << "' from current block " << currentBlock_->getDescription() << ".").endLine();
        }
        if (!keep_)
        {
          currentBlock_->getAlignment().removeSequence(i - 1);
        }
      }
      else
      {
        counts[species]++;
      }
    }
    bool test = currentBlock_->getNumberOfSequences() == 0;
    // Avoid a memory leak:
    if (test)
    {
      if (logstream_)
      {
        (*logstream_ << "SEQUENCE FILTER: block " << currentBlock_->getDescription() << " is now empty. Try to get the next one.").endLine();
      }
      delete currentBlock_;
    }
    else
    {
      test = strict_ && (counts.size() != species_.size());
      if (test)
      {
        if (logstream_)
        {
          (*logstream_ << "SEQUENCE FILTER: block " << currentBlock_->getDescription() << " does not contain all species and will be ignored. Try to get the next one.").endLine();
        }
        delete currentBlock_;
      }
      else
      {
        if (rmDuplicates_)
        {
          test = false;
          map<string, unsigned int>::iterator it;
          for (it = counts.begin(); it != counts.end() && !(test = it->second > 1); it++)
          {}
          if (test)
          {
            if (logstream_)
            {
              (*logstream_ << "SEQUENCE FILTER: block " << currentBlock_->getDescription() << " has two sequences for species '" << it->first << "' and will be ignored. Try to get the next one.").endLine();
            }
            delete currentBlock_;
          }
          else
          {
            return currentBlock_;
          }
        }
        else
        {
          return currentBlock_;
        }
      }
    }

    // Look for the next block:
    currentBlock_ = iterator_->nextBlock();
  }

  return currentBlock_;
}
