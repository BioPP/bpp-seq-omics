//
// File: CoordinateTranslatorMafIterator.h
// Authors: Julien Dutheil
// Created: Thu Jan 28 2016
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

#ifndef _COORDINATETRANSLATORMAFITERATOR_H_
#define _COORDINATETRANSLATORMAFITERATOR_H_

#include "MafIterator.h"

//From the STL:
#include <iostream>
#include <string>
#include <deque>

//From bpp-core:
#include <Bpp/Numeric/DataTable.h>

namespace bpp {

/**
 * @brief Translate features coordinates from one species to another, based on the alignment.
 *
 * This filter is similar in principle to the UCSC "liftOver" utility and software alike.
 * For now, only write a text file with all coordinates from reference and corresponding target sequence.
 */
class CoordinateTranslatorMafIterator:
  public AbstractFilterMafIterator
{
  private:
    std::string referenceSpecies_;
    std::string targetSpecies_;
    std::map<std::string, SequenceFeatureSet*> inputFeaturesPerChr_;
    std::ostream& output_;
    bool outputClosestCoordinate_;

  public:
    /**
     * @brief Build a new CoordinateTranslator iterator.
     *
     * @param iterator The input iterator
     * @param referenceSpecies The reference species for feature coordinates
     * @param targetSpecies The target species for which features coordinates should be translated
     * @param features The set of features to lift over
     * @param output Output stream for translated coordinates
     * @param outputClosestCoordinate In case the target sequence has a gap at the corresponding position,
     *        tells if the previous non-gap position should be returned, or NA.
     */
    CoordinateTranslatorMafIterator(
        MafIterator* iterator,
        const std::string& referenceSpecies,
        const std::string& targetSpecies,
        const SequenceFeatureSet& features,
        std::ostream& output,
        bool outputClosestCoordinate = true) :
      AbstractFilterMafIterator(iterator),
      referenceSpecies_(referenceSpecies),
      targetSpecies_(targetSpecies),
      inputFeaturesPerChr_(),
      output_(output),
      outputClosestCoordinate_(outputClosestCoordinate)
    {
      //Sort features per chromosome for a faster access:
      std::set<std::string> seqIds = features.getSequences();
      for (std::set<std::string>::iterator it = seqIds.begin();
          it != seqIds.end();
          ++it) {
        {
          inputFeaturesPerChr_[*it] = features.getSubsetForSequence(*it);
        }
      }
      output_ << "chr.ref\tstrand.ref\tbegin.ref\tend.ref\tchr.target\tstrand.target\tbegin.target\tend.target" << std::endl;
    }

    virtual ~CoordinateTranslatorMafIterator() {
      //Clean sorted features.
      for (std::map<std::string, SequenceFeatureSet*>::iterator it = inputFeaturesPerChr_.begin();
          it != inputFeaturesPerChr_.end();
          it++) {
        delete it->second;
      }
    }

  private:
    MafBlock* analyseCurrentBlock_();

};

} // end of namespace bpp.

#endif //_COORDINATETRANSLATORMAFITERATOR_H_
