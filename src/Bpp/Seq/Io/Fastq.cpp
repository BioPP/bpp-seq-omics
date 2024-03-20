// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "Fastq.h"
#include <Bpp/Seq/SequenceWithQuality.h>

#include <typeinfo>

using namespace bpp;

bool Fastq::nextSequence(std::istream& input, SequenceWithQuality& seq) const
{
  if (input && !input.eof())
  {
    // SequenceWithQuality& sq;
    std::string buffer;
    while (TextTools::isEmpty(buffer) && !input.eof())
    {
      getline(input, buffer);
    }
    if (input.eof())  // We hit the end of the file
    {
      return false;
    }
    // first line: seq name
    if (buffer[0] == '@')
    {
      seq.setName(std::string(buffer.begin() + 1, buffer.end()));
    }
    // second line: seq
    getline(input, buffer);
    seq.setContent(buffer);
    // third line: seq name (again)
    getline(input, buffer);
    if (repeatName())
    {
      std::string secName = std::string(buffer.begin() + 1, buffer.end());
      if (secName != seq.getName())
      {
        throw Exception("Names are not equivalent for sequence(@ line) and quality (+ line)");
      }
    }
    // fourth line: quality
    getline(input, buffer);
    try
    {
      SequenceWithQuality& sq = dynamic_cast<SequenceWithQuality&>(seq);
      for (size_t i = 0; i < buffer.size(); i++)
      {
        sq.setQuality(i, static_cast<int>(buffer[i]));
      }
    }
    catch (...)
    {}
    return true;
  }
  return false;
}

/******************************************************************************/

void Fastq::writeSequence(std::ostream& output, const SequenceWithQuality& seq) const
{
  std::string qual(seq.size(), static_cast<char>(33));
  try
  {
    const SequenceWithQuality& sq = dynamic_cast<const SequenceWithQuality&>(seq);
    for (size_t i = 0; i < sq.size(); i++)
    {
      char q = static_cast<char>(sq.getQuality(i));
      if (q < 33 || q > 126)
      {
        throw BadIntegerException("Quality must lie between 33 and 126", q);
      }
      qual[i] = q;
    }
  }
  catch (const std::bad_cast& e)
  {
    throw Exception("seq must be a SequenceWithQuality object");
  }
  output << "@" << seq.getName() << std::endl;
  output << seq.toString() << std::endl;
  output << "+";
  if (repeatName())
  {
    output << seq.getName();
  }
  output << std::endl;
  output << qual << std::endl;
}
