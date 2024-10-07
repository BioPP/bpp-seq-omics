// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "EstSfsOutputMafIterator.h"

// From bpp-seq:
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/SequenceWalker.h>

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>

using namespace std;

void EstSfsOutputMafIterator::writeBlock_(std::ostream& out, const MafBlock& block) const
{
  auto ingroupContainer = block.getAlignment(ingroup_);
  auto outgroup1Container = block.getAlignment(outgroup1_);
  auto outgroup2Container = block.getAlignment(outgroup2_);
  auto outgroup3Container = block.getAlignment(outgroup3_);

  for (size_t i = 0; i < ingroupContainer->getNumberOfSites(); ++i)
  {
    if (SymbolListTools::isComplete(ingroupContainer->site(i)))
    {
      map<int, size_t> counts;
      SiteTools::getCounts(ingroupContainer->site(i), counts);
      // Alphabet states are in alaphabetical order
      out << counts[0] << "," << counts[1] << "," << counts[2] << "," << counts[3];
    }
    else
    {
      out << "0,0,0,0";
    }

    if (SymbolListTools::isComplete(outgroup1Container->site(i)))
    {
      map<int, size_t> counts;
      SiteTools::getCounts(outgroup1Container->site(i), counts);
      out << " " << counts[0] << "," << counts[1] << "," << counts[2] << "," << counts[3];
    }
    else
    {
      out << " 0,0,0,0";
    }

    if (outgroup2_.size() > 0)
    {
      // Second outgroup is present
      if (SymbolListTools::isComplete(outgroup2Container->site(i)))
      {
        map<int, size_t> counts;
        SiteTools::getCounts(outgroup2Container->site(i), counts);
        out << " " << counts[0] << "," << counts[1] << "," << counts[2] << "," << counts[3];
      }
      else
      {
        out << " 0,0,0,0";
      }
    }

    if (outgroup3_.size() > 0)
    {
      // Third outgroup is present
      if (SymbolListTools::isComplete(outgroup3Container->site(i)))
      {
        map<int, size_t> counts;
        SiteTools::getCounts(outgroup3Container->site(i), counts);
        out << " " << counts[0] << "," << counts[1] << "," << counts[2] << "," << counts[3];
      }
      else
      {
        out << " 0,0,0,0";
      }
    }

    out << endl;
  }
}
