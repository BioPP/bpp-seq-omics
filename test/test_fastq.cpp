// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Seq/Io/Fastq.h>
#include <Bpp/Seq/SequenceWithQuality.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>

#include <iostream>
#include <fstream>
#include <memory>

using namespace bpp;
using namespace std;

int main ()
{
  try
  {
    std::string filename = "example.fastq";
    std::ifstream input(filename.c_str(), std::ios::in);
    if (!input)
    {
      std::cerr << "Could not open " << filename << std::endl;
      return 1;
    }
    Fastq fq;
    shared_ptr<const Alphabet> alpha = AlphabetTools::DNA_ALPHABET;
    SequenceWithQuality seq("", "", alpha);
    while (fq.nextSequence(input, seq))
    {
      std::cout << seq.getName() << " " << seq.size();
      std::cout << " " << VectorTools::min(seq.getQualities()) - 33;
      std::cout << " " << VectorTools::max(seq.getQualities()) - 33;
      std::cout << std::endl;
      fq.repeatName(true);
      fq.writeSequence(std::cout, seq);
      fq.repeatName(false);
    }
    return 0;
  }
  catch (std::exception& ex)
  {
    std::cerr << ex.what() << std::endl;
    return 1;
  }
}
