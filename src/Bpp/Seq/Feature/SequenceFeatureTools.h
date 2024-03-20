// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _SEQUENCEFEATURETOOLS_H_
#define _SEQUENCEFEATURETOOLS_H_

#include "SequenceFeature.h"

// From bpp-seq:
#include <Bpp/Seq/Sequence.h>
#include <Bpp/Seq/SequenceTools.h>
#include <Bpp/Seq/GeneticCode/GeneticCode.h>

#include <memory>

namespace bpp
{
class SequenceFeatureTools
{
public:
  /**
   * @brief Extract a sub-sequence given a SeqRange.
   *
   * The sub-sequence is reverse-complemented if SeqRange is in negative
   * strand.
   *
   * @param seq The Sequence to trunc.
   * @param range The SeqRange to extract.
   * @return A new Sequence object with the given subsequence oriented
   * according to the SeqRange.
   * @author Sylvain Gaillard
   */
  template<class SequenceTypeOut>
  static std::unique_ptr<SequenceTypeOut> extract(
      const SequenceInterface& seq,
      const SeqRange& range)
  {
    if (range.end() > seq.size())
      throw IndexOutOfBoundsException ("SequenceTools::extract: Invalid upper bound", range.end(), 0, seq.size());
    auto sout = SequenceTools::subseq<SequenceTypeOut>(seq, range.begin(), range.end() - 1);
    if (range.isNegativeStrand())
    {
      SequenceTools::invertComplement(*sout);
    }
    return sout;
  }

  /**
   * @brief Extract a sub-sequence given a SeqRange.
   *
   * The sub-sequence is revese-complemented if SeqRange is in negative
   * strand.
   *
   * @param seq The Sequence to trunc.
   * @param range The SeqRange to extract.
   * @param output The Sequence object to be filled with the given subsequence
   * oriented according to the SeqRange.
   * @author Sylvain Gaillard
   */
  static void extract(
      const SequenceInterface& seq,
      const SeqRange& range,
      SequenceInterface& output);

  /**
   * @brief Get ORF features for a Sequence.
   *
   * @param seq The Sequence where to find ORF. Must be a nucleic sequence.
   * @param featSet A SequenceFeatureSet to fill with the annotations.
   * @param gCode The genetic code to use.
   * @return The number of ORF found.
   * @author Sylvain Gaillard
   */
  static unsigned int getOrfs(
      const SequenceInterface& seq,
      SequenceFeatureSet& featSet,
      const GeneticCode& gCode);
};
} // end of namespace bpp

#endif // _SEQUENCEFEATURETOOLS_H_
