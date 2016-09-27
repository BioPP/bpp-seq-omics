//
// File: CoordinateTranslatorMafIterator.cpp
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

#include "CoordinateTranslatorMafIterator.h"

//From bpp-seq:
#include <Bpp/Seq/SequenceWalker.h>

using namespace bpp;

//From the STL:
#include <string>
#include <numeric>

using namespace std;

MafBlock* CoordinateTranslatorMafIterator::analyseCurrentBlock_() throw (Exception)
{
  unique_ptr<MafBlock> block(iterator_->nextBlock());
  if (!block.get()) return 0; //No more block.

  //Check if the block contains the reference and target species:
  if (!block->hasSequenceForSpecies(referenceSpecies_))
    return block.release();
  if (!block->hasSequenceForSpecies(targetSpecies_))
    return block.release();

  //Get the feature ranges for this block:
  const MafSequence& refSeq = block->getSequenceForSpecies(referenceSpecies_);
  const MafSequence& targetSeq = block->getSequenceForSpecies(targetSpecies_);

  //first check if there is one (for now we assume that features refer to the chromosome or contig name, with implicit species):
  std::map<std::string, SequenceFeatureSet*>::iterator mr = inputFeaturesPerChr_.find(refSeq.getChromosome());
  if (mr == inputFeaturesPerChr_.end())
    return block.release();
    
  //second get only features within this block:
  unique_ptr<SequenceFeatureSet> selectedFeatures(mr->second->getSubsetForRange(SeqRange(refSeq.getRange(true)), true)); 

  //test if there are some features to translate here:
  if (selectedFeatures->isEmpty())
    return block.release();

  //Get coordinate range sets:
  RangeSet<size_t> ranges;
  selectedFeatures->fillRangeCollection(ranges);

  //If the reference sequence is on the negative strand, then we have to correct the coordinates:
  if (refSeq.getStrand() == '-') {
    RangeSet<size_t> cRanges;
    for (set<Range<size_t>*>::iterator it = ranges.getSet().begin();
        it != ranges.getSet().end();
        ++it)
    {
      cRanges.addRange(SeqRange(refSeq.getSrcSize() - (**it).end(), refSeq.getSrcSize() - (**it).begin(), dynamic_cast<SeqRange*>(*it)->getStrand()));
    }
    ranges = cRanges;
  }

  //We will need to convert to alignment positions, using a sequence walker:
  SequenceWalker referenceWalker(refSeq);
  SequenceWalker targetWalker(targetSeq);

  //Now creates all blocks for all ranges:
  if (verbose_) {
    ApplicationTools::message->endLine();
    ApplicationTools::displayTask("Extracting annotations", true);
  }
  if (logstream_) {
    (*logstream_ << "COORDINATE CONVERTOR: lifting over " << ranges.getSet().size() << " features from block " << block->getDescription() << ".").endLine();
  }

  size_t i = 0;
  for (set<Range<size_t>*>::iterator it = ranges.getSet().begin();
      it !=  ranges.getSet().end();
      ++it)
  {
    if (verbose_) {
      ApplicationTools::displayGauge(i++, ranges.getSet().size() - 1, '=');
    }
    size_t a = referenceWalker.getAlignmentPosition((**it).begin() - refSeq.start());
    size_t b = referenceWalker.getAlignmentPosition((**it).end() - refSeq.start() - 1);
    size_t a2 = targetWalker.getSequencePosition(a) + targetSeq.start();
    size_t b2 = targetWalker.getSequencePosition(b) + targetSeq.start() + 1;
    if (targetSeq.getStrand() == '-') {
      a2 = targetSeq.getSrcSize() - a2;
      b2 = targetSeq.getSrcSize() - b2;
    }
    output_ << refSeq.getChromosome() << "\t" << refSeq.getStrand() << "\t" << (**it).begin() << "\t" << (**it).end() << "\t";
    output_ << targetSeq.getChromosome() << "\t" << targetSeq.getStrand() << "\t" << a2 << "\t" << b2 << endl;
  }
        
  if (verbose_)
    ApplicationTools::displayTaskDone();

  //Block is simply forwarded:
  return block.release();
}

