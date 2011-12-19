//
// File: test_sequences.cpp
// Created by: Julien Dutheil
// Created on: Mon Dec 130 17:10 2010
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

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

#include <Bpp/Seq/Alphabet.all>
#include <Bpp/Seq/SequenceTools.h>
#include <iostream>

using namespace bpp;
using namespace std;

int main() {
  //This is a very simple test that instanciate all alpahabet classes.
  BasicSequence seq1("test DNA", "ATTTCG---TCGTT-AAAGCACATGCATCGATC", &AlphabetTools::DNA_ALPHABET);
  BasicSequence motif1("motif", "ATTT", &AlphabetTools::DNA_ALPHABET);
  BasicSequence motif2("motif", "TCG", &AlphabetTools::DNA_ALPHABET);
  BasicSequence motif3("motif", "GATC", &AlphabetTools::DNA_ALPHABET);
  BasicSequence motif4("motif", "CGTC", &AlphabetTools::DNA_ALPHABET);
  unsigned int pos;

  pos = SequenceTools::findFirstOf(seq1, motif1);
  if (pos != 0) return 1;
  cout << motif1.toString() << ": " << pos << endl;
  
  pos = SequenceTools::findFirstOf(seq1, motif2);
  if (pos != 3) return 1;
  cout << motif2.toString() << ": " << pos << endl;

  pos = SequenceTools::findFirstOf(seq1, motif3);
  if (pos != 29) return 1;
  cout << motif3.toString() << ": " << pos << endl;

  pos = SequenceTools::findFirstOf(seq1, motif4);
  if (pos != 33) return 1;
  cout << motif4.toString() << ": " << pos << endl;

  return (0);
}
