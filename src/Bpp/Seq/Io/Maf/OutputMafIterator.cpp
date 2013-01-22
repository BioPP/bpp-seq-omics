//
// File: OutputMafIterator.cpp
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

#include "OutputMafIterator.h"

//From bpp-seq:
#include <Bpp/Seq/SequenceWithAnnotationTools.h>
#include <Bpp/Seq/SequenceWithQuality.h>

using namespace bpp;

//From the STL:
#include <string>
#include <numeric>

using namespace std;

void OutputMafIterator::writeHeader(std::ostream& out) const
{
  out << "##maf version=1 program=Bio++" << endl << "#" << endl;
  //There are more options in the header that we may want to support...
}

void OutputMafIterator::writeBlock(std::ostream& out, const MafBlock& block) const
{
  out << "a";
  if (! isinf(block.getScore()))
    out << " score=" << block.getScore();
  if (block.getPass() > 0)
    out << " pass=" << block.getPass();
  out << endl;
  
  //Now we write sequences. First need to count characters for aligning blocks:
  size_t mxcSrc = 0, mxcStart = 0, mxcSize = 0, mxcSrcSize = 0;
  for (size_t i = 0; i < block.getNumberOfSequences(); i++) {
    const MafSequence* seq = &block.getSequence(i);
    size_t start = 0; //Maybe we should output sthg else here?
    if (seq->hasCoordinates())
      start = seq->start();
    mxcSrc     = max(mxcSrc    , seq->getName().size());
    mxcStart   = max(mxcStart  , TextTools::toString(start).size());
    mxcSize    = max(mxcSize   , TextTools::toString(seq->getGenomicSize()).size());
    mxcSrcSize = max(mxcSrcSize, TextTools::toString(seq->getSrcSize()).size());
  }
  //Now print each sequence:
  for (size_t i = 0; i < block.getNumberOfSequences(); i++) {
    const MafSequence* seq = &block.getSequence(i);
    out << "s ";
    out << TextTools::resizeRight(seq->getName(), mxcSrc, ' ') << " ";
    size_t start = 0; //Maybe we should output sthg else here?
    if (seq->hasCoordinates())
      start = seq->start();
    out << TextTools::resizeLeft(TextTools::toString(start), mxcStart, ' ') << " ";
    out << TextTools::resizeLeft(TextTools::toString(seq->getGenomicSize()), mxcSize, ' ') << " ";
    out << seq->getStrand() << " ";
    out << TextTools::resizeLeft(TextTools::toString(seq->getSrcSize()), mxcSrcSize, ' ') << " ";
    //Shall we write the sequence as masked?
    string seqstr = seq->toString();
    if (mask_ && seq->hasAnnotation(SequenceMask::MASK)) {
      const SequenceMask* mask = &dynamic_cast<const SequenceMask&>(seq->getAnnotation(SequenceMask::MASK));
      for (size_t j = 0; j < seqstr.size(); ++j) {
        char c = ((*mask)[j] ? TextTools::toLower(seqstr[j]) : seqstr[j]);
        out << c;
      }
    } else {
      out << seqstr;
    }
    out << endl;
    //Write quality scores if any:
    if (mask_ && seq->hasAnnotation(SequenceQuality::QUALITY_SCORE)) {
      const SequenceQuality* qual = &dynamic_cast<const SequenceQuality&>(seq->getAnnotation(SequenceQuality::QUALITY_SCORE));
      out << "q ";
      out << TextTools::resizeRight(seq->getName(), mxcSrc + mxcStart + mxcSize + mxcSrcSize + 5, ' ') << " ";
      string qualStr;
      for (size_t j = 0; j < seq->size(); ++j) {
        int s = (*qual)[j];
        if (s == -1) {
          qualStr += "-";
        } else if (s == -2) {
          qualStr += "?";
        } else if (s >=0 && s < 10) {
          qualStr += TextTools::toString(s);
        } else if (s == 10) {
          qualStr += "F";
        } else {
          throw Exception("MafAlignmentParser::writeBlock. Unsuported score value: " + TextTools::toString(s));
        }
      }
      out << qualStr << endl;
    }
  }
  out << endl;
}

