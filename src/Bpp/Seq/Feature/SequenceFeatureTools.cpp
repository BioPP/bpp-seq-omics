// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "SequenceFeatureTools.h"

// From bpp-seq:
#include <Bpp/Seq/SequenceTools.h>
#include <Bpp/Seq/Alphabet/CodonAlphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Alphabet/AlphabetExceptions.h>

// From STL
#include <vector>

using namespace bpp;
using namespace std;

/******************************************************************************/

void SequenceFeatureTools::extract(
    const SequenceInterface& seq,
    const SeqRange& range, 
    SequenceInterface& output)
{
  if (range.end() > seq.size())
    throw IndexOutOfBoundsException ("SequenceTools::extract: Invalid upper bound", range.end(), 0, seq.size());
  SequenceTools::subseq(seq, range.begin(), range.end() - 1, output);
  if (range.isNegativeStrand())
  {
    SequenceTools::invertComplement(output);
  }
}

/******************************************************************************/

unsigned int SequenceFeatureTools::getOrfs(
    const SequenceInterface& seq,
    SequenceFeatureSet& featSet,
    const GeneticCode& gCode)
{
  if (!AlphabetTools::isNucleicAlphabet(&seq.alphabet()))
  {
    throw AlphabetException("SequenceFeatureTools::getOrfs: Sequence alphabet must be nucleic!", seq.getAlphabet());
  }
  unsigned int orfCpt = 0;
  auto codonAlpha = gCode.getCodonAlphabet();
  std::vector< std::vector<size_t>> starts(3), stops(3);
  size_t phase = 0;
  for (size_t p = 0; p < seq.size() - 2; p++)
  {
    phase = p % 3;
    if (gCode.isStart(codonAlpha->getCodon(seq.getValue(p), seq.getValue(p + 1), seq.getValue(p + 2))))
    {
      starts[phase].push_back(p);
      // std::cerr << "Start: " << p << " (" << phase << ")" << std::endl;
    }
    else if (gCode.isStop(codonAlpha->getCodon(seq.getValue(p), seq.getValue(p + 1), seq.getValue(p + 2))))
    {
      stops[phase].push_back(p);
      // std::cerr << "Stop:  " << p << " (" << phase << ")" << std::endl;
    }
  }
  for (size_t i = 0; i < 3; ++i)
  {
    std::vector< size_t >::iterator start(starts[i].begin()), stop(stops[i].begin());
    while (stop != stops[i].end() && start != starts[i].end())
    {
      if (*stop < *start)
      {
        stop++;
      }
      else
      {
        orfCpt++;
        // std::cerr << "ORF:  " << *start << " - " << *stop + 2 << " (" << i << ")" << std::endl;
        bpp::BasicSequenceFeature feat("", seq.getName(), "Bio++", "CDS", *start, *stop + 2, '+');
        featSet.addFeature(feat);
        start++;
      }
    }
  }
  return orfCpt;
}

/******************************************************************************/
