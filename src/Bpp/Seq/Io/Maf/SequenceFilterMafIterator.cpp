// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "SequenceFilterMafIterator.h"

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>

using namespace std;

unique_ptr<MafBlock> SequenceFilterMafIterator::analyseCurrentBlock_()
{
  currentBlock_ = iterator_->nextBlock();
  while (currentBlock_)
  {
    map<string, unsigned int> counts;
    for (size_t i = currentBlock_->getNumberOfSequences(); i > 0; --i)
    {
      string species = currentBlock_->sequence(i - 1).getSpecies();
      if (!VectorTools::contains(species_, species))
      {
        if (!keep_)
        {
          if (logstream_)
          {
            (*logstream_ << "SEQUENCE FILTER: remove sequence '" << species << "' from current block " << currentBlock_->getDescription() << ".").endLine();
          }
          currentBlock_->removeSequence(i - 1);
        }
      }
      else
      {
        counts[species]++;
      }
    }
    bool test = currentBlock_->getNumberOfSequences() == 0;
    if (test)
    {
      if (logstream_)
      {
        (*logstream_ << "SEQUENCE FILTER: block " << currentBlock_->getDescription() << " is now empty. Try to get the next one.").endLine();
      }
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
            if (duplicatePolicy_ == DUPLICATE_POLICY_DISCARD)
            {
              if (logstream_)
              {
                (*logstream_ << "SEQUENCE FILTER: block " << currentBlock_->getDescription() << " has two sequences for species '" << it->first << "' and will be ignored. Try to get the next one.").endLine();
              }
              else
              {
                return std::move(currentBlock_);
              }
            }
            else if(duplicatePolicy_ == DUPLICATE_POLICY_SAMPLE)
            {
              map<string, unsigned int> counts2;
              for (size_t i = currentBlock_->getNumberOfSequences(); i > 0; --i)
              {
                string species = currentBlock_->sequence(i - 1).getSpecies();
		counts2[species]++;
		if (counts2[species] > 1)
		{
                  if (logstream_)
                  {
                    (*logstream_ << "SEQUENCE FILTER: remove duplicated sequence '" << species << "' from current block " << currentBlock_->getDescription() << ".").endLine();
                  }
                  currentBlock_->removeSequence(i - 1);
		}
              }
              return std::move(currentBlock_);
            }
            else if(duplicatePolicy_ == DUPLICATE_POLICY_REMOVE)
            {
              for (it = counts.begin(); it != counts.end(); it++)
              {
	        if (it->second > 1)
		{
                  currentBlock_->deleteSequencesForSpecies(it->first);
		}
	      }	       
              return std::move(currentBlock_);
            }
            else
            {
              throw Exception("SequenceFilterMafIterator::analyseCurrentBlock_: invalid DUPLICATE_POLICY.");
            }
          }
        }
        else
        {
          return std::move(currentBlock_);
        }
      }
    }

    // Look for the next block:
    currentBlock_ = iterator_->nextBlock();
  }

  return std::move(currentBlock_);
}
