// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "SequenceStatisticsMafIterator.h"

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>

using namespace std;

SequenceStatisticsMafIterator::SequenceStatisticsMafIterator(
    std::shared_ptr<MafIteratorInterface> iterator,
    const std::vector<shared_ptr<MafStatisticsInterface>> statistics) :
  AbstractFilterMafIterator(iterator),
  statistics_(statistics),
  results_(),
  names_()
{
  string name;
  for (size_t i = 0; i < statistics_.size(); ++i)
  {
    name = statistics_[i]->getShortName();
    vector<string> tags = statistics_[i]->getSupportedTags();
    if (tags.size() > 1)
    {
      for (size_t j = 0; j < tags.size(); ++j)
      {
        names_.push_back(name + "." + tags[j]);
      }
    }
    else
    {
      names_.push_back(name);
    }
  }
  results_.resize(names_.size());
}

unique_ptr<MafBlock> SequenceStatisticsMafIterator::analyseCurrentBlock_()
{
  vector<string> tags;
  currentBlock_ = iterator_->nextBlock();
  if (currentBlock_)
  {
    size_t k = 0;
    for (size_t i = 0; i < statistics_.size(); ++i)
    {
      statistics_[i]->compute(*currentBlock_);
      const MafStatisticsResult& result = statistics_[i]->getResult();
      tags = statistics_[i]->getSupportedTags();
      for (size_t j = 0; j < tags.size(); ++j)
      {
        if (result.hasValue(tags[j]))
        {
          results_[k].reset(result.getValue(tags[j]).clone());
        }
        else
        {
          results_[k] = 0;
        }
        k++;
      }
    }
  }
  return move(currentBlock_);
}
