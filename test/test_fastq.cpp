// 
// File:    test_fastq.cpp
// Author:  Sylvain Gaillard
// Created: 22/11/2011 11:29:16
// 

/*
Copyright or © or Copr. Bio++ Development Team, (November 22, 2011)

This software is a computer program whose purpose is to provide classes
for numerical calculus. This file is part of the Bio++ project.

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

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Seq/Io/Fastq.h>
#include <Bpp/Seq/SequenceWithQuality.h>
#include <Bpp/Seq/Alphabet/DNA.h>

#include <iostream>
#include <fstream>

using namespace bpp;

int main () {
  try {
    std::string filename = "example.fastq";
    std::ifstream input(filename.c_str(), std::ios::in);
    if (!input) {
      std::cerr << "Could not open " << filename << std::endl;
      return 1;
    }
    Fastq fq;
    const Alphabet* alpha = new DNA();
    SequenceWithQuality seq("", "", alpha);
    while (fq.nextSequence(input, seq)) {
      std::cout << seq.getName() << " " << seq.size();
      std::cout << " " << VectorTools::min(seq.getQualities()) - 33;
      std::cout << " " << VectorTools::max(seq.getQualities()) - 33;
      std::cout << std::endl;
      fq.repeatName(true);
      fq.writeSequence(std::cout, seq);
      fq.repeatName(false);
    }
    return 0;
  } catch (std::exception& ex) {
    std::cerr << ex.what() << std::endl;
    return 1;
  }
}
