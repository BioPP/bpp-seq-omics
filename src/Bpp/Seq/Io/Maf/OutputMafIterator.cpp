// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "OutputMafIterator.h"

// From bpp-seq:
#include <Bpp/Seq/SequenceWithAnnotationTools.h>
#include <Bpp/Seq/SequenceWithQuality.h>

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>

using namespace std;

void OutputMafIterator::writeHeader(std::ostream& out) const
{
  out << "##maf version=1 program=Bio++" << endl << "#" << endl;
  // There are more options in the header that we may want to support...
}

void OutputMafIterator::writeBlock(std::ostream& out, const MafBlock& block) const
{
  out << "a";
  if (!std::isinf(block.getScore()))
    out << " score=" << block.getScore();
  if (block.getPass() > 0)
    out << " pass=" << block.getPass();
  out << endl;

  // Now we write sequences. First need to count characters for aligning blocks:
  size_t mxcSrc = 0, mxcStart = 0, mxcSize = 0, mxcSrcSize = 0;
  for (size_t i = 0; i < block.getNumberOfSequences(); ++i)
  {
    const MafSequence& seq = block.sequence(i);
    size_t start = 0; // Maybe we should output sthg else here?
    if (seq.hasCoordinates())
      start = seq.start();
    mxcSrc     = max(mxcSrc, seq.getName().size());
    mxcStart   = max(mxcStart, TextTools::toString(start).size());
    mxcSize    = max(mxcSize, TextTools::toString(seq.getGenomicSize()).size());
    mxcSrcSize = max(mxcSrcSize, TextTools::toString(seq.getSrcSize()).size());
  }
  // Now print each sequence:
  for (size_t i = 0; i < block.getNumberOfSequences(); ++i)
  {
    const MafSequence& seq = block.sequence(i);
    out << "s ";
    out << TextTools::resizeRight(seq.getName(), mxcSrc, ' ') << " ";
    size_t start = 0; // Maybe we should output sthg else here?
    if (seq.hasCoordinates())
      start = seq.start();
    out << TextTools::resizeLeft(TextTools::toString(start), mxcStart, ' ') << " ";
    out << TextTools::resizeLeft(TextTools::toString(seq.getGenomicSize()), mxcSize, ' ') << " ";
    out << seq.getStrand() << " ";
    out << TextTools::resizeLeft(TextTools::toString(seq.getSrcSize()), mxcSrcSize, ' ') << " ";
    // Shall we write the sequence as masked?
    string seqstr = seq.toString();
    if (mask_ && seq.hasAnnotation(SequenceMask::MASK))
    {
      const SequenceMask& mask = dynamic_cast<const SequenceMask&>(seq.annotation(SequenceMask::MASK));
      for (size_t j = 0; j < seqstr.size(); ++j)
      {
        char c = ((mask)[j] ? static_cast<char>(tolower(static_cast<int>(seqstr[j]))) : seqstr[j]);
        out << c;
      }
    }
    else
    {
      out << seqstr;
    }
    out << endl;
    // Write quality scores if any:
    if (mask_ && seq.hasAnnotation(SequenceQuality::QUALITY_SCORE))
    {
      const SequenceQuality& qual = dynamic_cast<const SequenceQuality&>(seq.annotation(SequenceQuality::QUALITY_SCORE));
      out << "q ";
      out << TextTools::resizeRight(seq.getName(), mxcSrc + mxcStart + mxcSize + mxcSrcSize + 5, ' ') << " ";
      string qualStr;
      for (size_t j = 0; j < seq.size(); ++j)
      {
        int s = (qual)[j];
        if (s == -1)
        {
          qualStr += "-";
        }
        else if (s == -2)
        {
          qualStr += "?";
        }
        else if (s >= 0 && s < 10)
        {
          qualStr += TextTools::toString(s);
        }
        else if (s == 10)
        {
          qualStr += "F";
        }
        else
        {
          throw Exception("MafAlignmentParser::writeBlock. Unsupported score value: " + TextTools::toString(s));
        }
      }
      out << qualStr << endl;
    }
  }
  out << endl;
}
