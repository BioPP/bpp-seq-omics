//
// File: MafAlignmentParser.h
// Authors: Julien Dutheil
// Created: Tue Apr 27 2010
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

#include "SequenceStreamToMafIterator.h"
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>

using namespace std;
using namespace bpp;

unique_ptr<MafBlock> SequenceStreamToMafIterator::analyseCurrentBlock_()
{
  auto block = make_unique<MafBlock>();

  auto mafSeq = make_unique<MafSequence>();
  if (stream_->eof())
    return nullptr;

  seqStream_->nextSequence(*stream_, *mafSeq);
  // Check if sequence name contains meta information:
  string meta = mafSeq->getName();
  StringTokenizer st(meta, ":");
  if (st.numberOfRemainingTokens() == 5)
  {
    string species = st.nextToken();
    string chr     = st.nextToken();
    unsigned int start   = TextTools::to<unsigned int>(st.nextToken());
    if (!zeroBasedCoords_)
      start--;
    string strand  = st.nextToken();
    unsigned int length  = TextTools::to<unsigned int>(st.nextToken());
    mafSeq->setName(species + "." + chr);
    mafSeq->setStrand(strand[0]);
    mafSeq->setStart(start);
    if (mafSeq->size() != length)
      throw Exception("SequenceStreamToMafIterator::analyseCurrentBlock_. Sequence size does not match its header specification: expected " + TextTools::toString(length) + " and found " + TextTools::toString(mafSeq->size()));
  }
  block->addSequence(mafSeq->getName(), mafSeq);

  return block;
}
