//
// File: SequenceStatisticsMafIterator.cpp
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

#include "SequenceStatisticsMafIterator.h"

using namespace bpp;

//From the STL:
#include <string>
#include <numeric>

using namespace std;

SequenceStatisticsMafIterator::SequenceStatisticsMafIterator(MafIterator* iterator, const std::vector<MafStatistics*> statistics) :
  AbstractFilterMafIterator(iterator),
  statistics_(statistics),
  results_(),
  names_()
{
  string name;
  for (size_t i = 0; i < statistics_.size(); ++i) {
    name = statistics_[i]->getShortName();
    vector<string> tags = statistics_[i]->getSupportedTags();
    if (tags.size() > 1) {
      for (size_t j = 0; j < tags.size(); ++j) {
        names_.push_back(name + "." + tags[j]);
      }
    } else {
      names_.push_back(name);
    }
  }
  results_.resize(names_.size());
}

MafBlock* SequenceStatisticsMafIterator::analyseCurrentBlock_() throw (Exception)
{
  vector<string> tags;
  currentBlock_ = iterator_->nextBlock();
  if (currentBlock_) {
    size_t k = 0;
    for (size_t i = 0; i < statistics_.size(); ++i) {
      statistics_[i]->compute(*currentBlock_);
      const MafStatisticsResult& result = statistics_[i]->getResult();
      tags = statistics_[i]->getSupportedTags();
      for (size_t j = 0; j < tags.size(); ++j) {
        if (result.hasValue(tags[j])) {
          results_[k] = result.getValue(tags[j]);
        } else {
          results_[k] = NumConstants::NaN;
        }
        k++;
      }
    }
  }
  return currentBlock_;
}

