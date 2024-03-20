// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "SequenceStreamToMafIterator.h"
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>

using namespace std;
using namespace bpp;

unique_ptr<MafBlock> SequenceStreamToMafIterator::analyseCurrentBlock_()
{
  auto block = make_unique<MafBlock>();

  shared_ptr<const Alphabet> alpha = AlphabetTools::DNA_ALPHABET;
  auto seq = make_unique<Sequence>(alpha);
  if (stream_->eof())
    return nullptr;

  seqStream_->nextSequence(*stream_, *seq);
  // Check if sequence name contains meta information:
  string meta = seq->getName();
  StringTokenizer st(meta, ":");
  auto mafSeq = make_unique<MafSequence>(*seq);
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
  block->addSequence(mafSeq);

  return block;
}
