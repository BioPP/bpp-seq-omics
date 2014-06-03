//
// File: SequenceStatisticsMafIterator.h
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

#ifndef _SEQUENCESTATISTICSMAFITERATOR_H_
#define _SEQUENCESTATISTICSMAFITERATOR_H_

#include "MafIterator.h"
#include "MafStatistics.h"

//From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp {

/**
 * @brief Compute a series of sequence statistics for each block.
 *
 * Computed statistics are stored into a vector of double, which can be retrieved as well as statistics names.
 * Listeners can be set up to automatically analyse or write the output after iterations are over.
 *
 * The current implementation focuses on speed and memory efificiency, as it only stores in memory the current results of the statistics.
 * The only drawback of this, is that disk access might be high when writing the results,
 * although appropriate buffering should most likely circumvent the issue.
 * The code is easily extensible, however, to enable storage of all results into a matrix,
 * with writing only once at the end of iterations.
 */
class SequenceStatisticsMafIterator:
  public AbstractFilterMafIterator
{
  private:
    std::vector<MafStatistics*> statistics_;
    std::vector<const BppNumberI*> results_;
    std::vector<std::string> names_;

  public:
    /**
     * @param iterator The input iterator.
     * @param statistics A vector of pointers toward MafStatistics, to be computed simultaneously for each maf block.
     */
    SequenceStatisticsMafIterator(MafIterator* iterator, const std::vector<MafStatistics*> statistics);

  private:
    SequenceStatisticsMafIterator(const SequenceStatisticsMafIterator& iterator) :
      AbstractFilterMafIterator(0),
      statistics_(iterator.statistics_),
      results_(iterator.results_),
      names_(iterator.names_)
    {}
    
    SequenceStatisticsMafIterator& operator=(const SequenceStatisticsMafIterator& iterator)
    {
      statistics_ = iterator.statistics_;
      results_ = iterator.results_;
      names_ = iterator.names_;
      return *this;
    }

  public:
    const std::vector<const BppNumberI*>& getResults() const { return results_; }
    const std::vector<std::string>& getResultsColumnNames() const { return names_; }

  private:
    MafBlock* analyseCurrentBlock_() throw (Exception);

};

} // end of namespace bpp.

#endif //_SEQUENCESTATISTICSMAFITERATOR_H_
