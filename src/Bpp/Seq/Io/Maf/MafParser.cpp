// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "MafParser.h"
#include <Bpp/Seq/SequenceWithQuality.h>
#include <Bpp/Seq/SequenceWithAnnotationTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>

#include <algorithm>

using namespace std;
using namespace bpp;

std::unique_ptr<MafBlock> MafParser::analyseCurrentBlock_()
{
  unique_ptr<MafBlock> block = nullptr;

  string line;
  bool test = true;
  unique_ptr<MafSequence> currentSequence;

  while (test)
  {
    if (stream_->eof())
    {
      break;
    }
    getline(*stream_, line, '\n');
    if (TextTools::isEmpty(line))
    {
      if (firstBlock_)
        continue;
      if (currentSequence)
      {
        // Add previous sequence:
        block->addSequence(currentSequence);
      }

      // end of paragraph
      test = false;
      firstBlock_ = true; //This allows to have multiple blanck lines between blocks, not only at the start of the file.
    }
    else if (line[0] == 'a')
    {
      if (currentSequence)
      {
        // Add previous sequence:
        block->addSequence(currentSequence);
      }

      // New block.
      block = make_unique<MafBlock>();
      firstBlock_ = false;

      map<string, string> args;
      if (line.size() > 2)
      {
        KeyvalTools::multipleKeyvals(line.substr(2), args, " ");

        if (args.find("score") != args.end())
          if (args["score"] != "NA")
            block->setScore(TextTools::toDouble(args["score"]));

        if (args.find("pass") != args.end())
          block->setPass(TextTools::to<unsigned int>(args["pass"]));
      }
    }
    else if (line[0] == 's')
    {
      StringTokenizer st(line);
      st.nextToken(); // The 's' tag
      if (!st.hasMoreToken())
        throw IOException("Sequence description should include a source field.");
      string src = st.nextToken();
      if (!st.hasMoreToken())
        throw IOException("Sequence description should include a start field.");
      unsigned int start = TextTools::to<unsigned int>(st.nextToken());
      if (!st.hasMoreToken())
        throw IOException("Sequence description should include a size field.");
      unsigned int size = TextTools::to<unsigned int>(st.nextToken());
      if (!st.hasMoreToken())
        throw IOException("Sequence description should include a strand field.");
      string tmp = st.nextToken();
      if (tmp.size() != 1)
        throw Exception("MafAlignmentParser::nextBlock. Strand specification is incorrect, should be only one character long, found " + TextTools::toString(tmp.size()) + ".");
      char strand = tmp[0];

      if (!st.hasMoreToken())
        throw IOException("Sequence description should include a source size field.");
      unsigned int srcSize = TextTools::to<unsigned int>(st.nextToken());
      if (currentSequence)
      {
        // Add previous sequence:
        block->addSequence(currentSequence);
      }
      if (!st.hasMoreToken())
        throw IOException("Sequence description without a sequence.");
      string seq = st.nextToken();
      if (dotOption_ == DOT_ASGAP)
      {
        std::replace(seq.begin(), seq.end(), '.', '-');
      }
      if (dotOption_ == DOT_ASUNRES)
      {
        std::replace(seq.begin(), seq.end(), '.', 'N');
      }
      currentSequence.reset(new MafSequence(src, seq, start, strand, srcSize));
      if (currentSequence->getGenomicSize() != size)
      {
        if (checkSequenceSize_)
          throw Exception("MafAlignmentParser::nextBlock. Sequence found (" + src + ") does not match specified size: " + TextTools::toString(currentSequence->getGenomicSize()) + ", should be " + TextTools::toString(size) + ".");
        else
        {
          if (verbose_)
          {
            ApplicationTools::displayWarning("MafAlignmentParser::nextBlock. Sequence found (" + src + ") does not match specified size: " + TextTools::toString(currentSequence->getGenomicSize()) + ", should be " + TextTools::toString(size) + ".");
          }
        }
      }
      // Add mask:
      if (mask_)
      {
        vector<bool> mask(currentSequence->size());
        for (size_t i = 0; i < mask.size(); ++i)
        {
          mask[i] = cmAlphabet_.isMasked(seq[i]);
        }
        currentSequence->addAnnotation(make_shared<SequenceMask>(mask));
      }
    }
    else if (line[0] == 'q')
    {
      if (!currentSequence)
        throw Exception("MafAlignmentParser::nextBlock(). Quality scores found, but there is currently no sequence!");
      StringTokenizer st(line);
      st.nextToken(); // The 'q' tag
      string name = st.nextToken();
      if (name != currentSequence->getName())
        throw Exception("MafAlignmentParser::nextBlock(). Quality scores found, but with a different name from the previous sequence: " + name + ", should be " + currentSequence->getName() + ".");
      string qstr = st.nextToken();
      // Now parse the score string:
      auto seqQual = make_shared<SequenceQuality>(qstr.size());
      for (size_t i = 0; i < qstr.size(); ++i)
      {
        char c = qstr[i];
        if (c == '-')
        {
          seqQual->setScore(i, -1);
        }
        else if (c == '0' || c == '1' || c == '2' || c == '3' || c == '4' || c == '5' || c == '6' || c == '7' || c == '8' || c == '9')
        {
          seqQual->setScore(i, c - '0');
        }
        else if (c == 'F' || c == 'f')  // Finished
        {
          seqQual->setScore(i, 10);
        }
        else if (c == '?' || c == '.')
        {
          seqQual->setScore(i, -2);
        }
        else
        {
          throw Exception("MafAlignmentParser::nextBlock(). Invalid quality score: " + TextTools::toString(c) + ". Should be 0-9, F or '-'.");
        }
      }
      currentSequence->addAnnotation(seqQual);
    }
  }
  //// Final check and passing by results

  // In case last line in not empty:
  if (currentSequence)
  {
    // Add previous sequence:
    block->addSequence(currentSequence);
  }

  // Returning block:
  return block;
}
