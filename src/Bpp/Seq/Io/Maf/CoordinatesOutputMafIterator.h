//
// File: CoordinatesOutputMafIterator.h
// Authors: Julien Dutheil
// Created: Mon Jun 02 2014
//

/*
Copyright or © or Copr. Bio++ Development Team, (2014)

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

#ifndef _COORDINATESOUTPUTMAFITERATOR_H_
#define _COORDINATESOUTPUTMAFITERATOR_H_

#include "MafIterator.h"

//From the STL:
#include <iostream>
#include <string>

namespace bpp {

/**
 * @brief Output each sequence coordinates for each block.
 *
 * The set of species for which coordinates are output is provided as argument.
 * The current implementation outputs results as a table to a file. Later implementation
 * may involve a dedicated data structure for other application usage (file indexing and so one).
 */
class CoordinatesOutputMafIterator:
  public AbstractFilterMafIterator
{
  private:
    std::ostream* output_;
    std::vector<std::string> species_;
    bool includeSrcSize_;

  public:
    /**
     * @brief Creates a new CoordinatesOutputMafIterator object.
     *
     * @param iterator The input iterator.
     * @param out A pointer toward the output stream. The stream will not be own by this instance, and will not be copied neither destroyed.
     * @param species A vector of species names for which coordinates should be output. In case of missing species for one block, NA will be produced.
     * @param includeSrcSize Tell if source size should also be written (useful to convert coordinates on the negative strand).
     */
    CoordinatesOutputMafIterator(MafIterator* iterator, std::ostream* out, const std::vector<std::string>& species, bool includeSrcSize = false):
      AbstractFilterMafIterator(iterator), output_(out), species_(species), includeSrcSize_(includeSrcSize)
    {
      if (output_)
        writeHeader_(*output_);
    }

  private:
    CoordinatesOutputMafIterator(const CoordinatesOutputMafIterator& iterator) :
      AbstractFilterMafIterator(0),
      output_(iterator.output_),
      species_(iterator.species_),
      includeSrcSize_(iterator.includeSrcSize_)
    {}
    
    CoordinatesOutputMafIterator& operator=(const CoordinatesOutputMafIterator& iterator)
    {
      output_ = iterator.output_;
      species_ = iterator.species_;
      includeSrcSize_ = iterator.includeSrcSize_;
      return *this;
    }

  private:
    void writeHeader_(std::ostream& out) const;
    MafBlock* analyseCurrentBlock_();

};

} // end of namespace bpp.

#endif //_COORDINATESOUTPUTMAFITERATOR_H_
