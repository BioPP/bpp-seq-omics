// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Seq/Feature/Gff/GffFeatureReader.h>

#include <iostream>
#include <fstream>

using namespace bpp;
using namespace std;

int main() {
  try {
  ifstream input("example.gff", ios::in);
  GffFeatureReader reader(input);
  while (reader.hasMoreFeature()) {
    BasicSequenceFeature feature = reader.nextFeature();
    cout << "Found feature " << feature.getId() << " of type " << feature.getType() << ", starting at " << feature.getStart() << " and ending at " << feature.getEnd() << endl;
  }
  return 0;
  } catch (exception& ex) {
    cerr << ex.what() << endl;
    return 1;
  }
}
