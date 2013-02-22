//
// File: MafStatistics.cpp
// Authors: Julien Dutheil
// Created: Mon Jun 25 2012
//

/*
Copyright or © or Copr. Bio++ Development Team, (2010)

This software is a computer program whose purpose is to provide classes
for sequences analysis.

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

#include "MafStatistics.h"
#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/SiteTools.h>

//From bpp-core:
#include <Bpp/Numeric/NumConstants.h>

//From the STL:
#include <cmath>
#include <map>

using namespace bpp;
using namespace std;

void PairwiseDivergenceMafStatistics::compute(const MafBlock& block)
{
  vector<const MafSequence*> seqs1 = block.getSequencesForSpecies(species1_);
  vector<const MafSequence*> seqs2 = block.getSequencesForSpecies(species2_);
  if (seqs1.size() > 1 || seqs2.size() > 1)
    throw Exception("PairwiseDivergenceMafStatistics::compute. Duplicated sequence for species " + species1_ + "or " + species2_ + ".");
  if (seqs1.size() == 0 || seqs2.size() == 0)
    result_.setValue(NumConstants::NaN());
  else
    result_.setValue(100. - SequenceTools::getPercentIdentity(*seqs1[0], *seqs2[0], true));
}

vector<string> CharacterCountsMafStatistics::getSupportedTags() const
{
  vector<string> tags;
  for (int i = 0; i < static_cast<int>(alphabet_->getSize()); ++i) {
    tags.push_back(alphabet_->intToChar(i));
  }
  tags.push_back("Gap");
  tags.push_back("Unresolved");

  return tags;
}

void CharacterCountsMafStatistics::compute(const MafBlock& block)
{
  std::map<int, int> counts;
  SequenceContainerTools::getCounts(block.getAlignment(), counts); 
  for (int i = 0; i < static_cast<int>(alphabet_->getSize()); ++i) {
    result_.setValue(alphabet_->intToChar(i), counts[i]);
  }
  result_.setValue("Gap", counts[alphabet_->getGapCharacterCode()]);
  double countUnres = 0;
  for (map<int, int>::iterator it = counts.begin(); it != counts.end(); ++it) {
    if (alphabet_->isUnresolved(it->first))
      countUnres += it->second;
  }
  result_.setValue("Unresolved", countUnres);
}

SiteContainer* AbstractSpeciesSelectionMafStatistics::getSiteContainer(const MafBlock& block)
{
  VectorSiteContainer* alignment = new VectorSiteContainer(block.getAlignment().getAlphabet());
  for (size_t i = 0; i < species_.size(); ++i) {
    if (block.hasSequenceForSpecies(species_[i])) {
      vector<const MafSequence*> selection = block.getSequencesForSpecies(species_[i]);
      for (size_t j = 0; j < selection.size(); ++j) {
        alignment->addSequence(*selection[j]);
      }
    }
  }
  return alignment;
}

vector<string> SiteFrequencySpectrumMafStatistics::getSupportedTags() const
{
  vector<string> tags;
  for (size_t i = 0; i < categorizer_.getNumberOfCategories(); ++i) {
    tags.push_back("Bin" + TextTools::toString(i + 1)); 
  }
  tags.push_back("Unresolved");
  tags.push_back("Saturated");
  tags.push_back("Ignored");
  return tags;
}

void SiteFrequencySpectrumMafStatistics::compute(const MafBlock& block)
{
  unsigned int nbUnresolved = 0;
  unsigned int nbSaturated = 0;
  unsigned int nbIgnored = 0;
  counts_.assign(categorizer_.getNumberOfCategories(), 0);
  int state;
  auto_ptr<SiteContainer> alignment(getSiteContainer(block));
  if (alignment->getNumberOfSequences() > 0) {
    for (size_t i = 0; i < block.getNumberOfSites(); ++i) {
      //Note: we do not rely on SiteTool::getCounts as it would be unefficient to count everything.
      const Site& site = alignment->getSite(i);
      map<int, unsigned int> counts;
      for (unsigned int j = 0; j < site.size(); ++j) {
        state = site[j];
        if (alphabet_->isGap(state) || alphabet_->isUnresolved(state)) {
          nbUnresolved++;
          break;
        } else {
          counts[state]++;
          if (counts.size() > 2) {
            nbSaturated++;
            break;
          }
        }
      }
      if (counts.size() > 0) {
        //Determine frequency class:
        double count;
        if (counts.size() == 1) {
          count = 0;
        } else {
          map<int, unsigned int>::iterator it = counts.begin();
          unsigned int count1 = it->second;
          it++;
          unsigned int count2 = it->second;
          count = min(count1, count2);
        }
        try {
          counts_[categorizer_.getCategory(count) - 1]++;
        } catch (OutOfRangeException& oof) {
          nbIgnored++;
        }
      }
    }
  }
  result_.setValue("Unresolved", nbUnresolved);
  result_.setValue("Saturated", nbSaturated);
  result_.setValue("Ignored", nbIgnored);
  for (size_t i = 0; i < counts_.size(); ++i) {
    result_.setValue("Bin" + TextTools::toString(i + 1), counts_[i]);
  }
}

vector<string> SiteMafStatistics::getSupportedTags() const
{
  vector<string> tags;
  tags.push_back("NbWithoutGap");
  tags.push_back("NbComplete");
  tags.push_back("NbParsimonyInformative");
  return tags;
}

void SiteMafStatistics::compute(const MafBlock& block)
{
  auto_ptr<SiteContainer> alignment(getSiteContainer(block));
  unsigned int nbNg = 0;
  unsigned int nbCo = 0;
  unsigned int nbPi = 0;
  if (alignment->getNumberOfSequences() > 0) {
    for (unsigned int i = 0; i < alignment->getNumberOfSites(); ++i) {
      if (!SiteTools::hasGap(alignment->getSite(i)))
        nbNg++;
      if (SiteTools::isComplete(alignment->getSite(i)))
        nbCo++;
      if (SiteTools::isParsimonyInformativeSite(alignment->getSite(i)))
        nbPi++;
    }
  }
  result_.setValue("NbWithoutGap", nbNg);
  result_.setValue("NbComplete", nbCo);
  result_.setValue("NbParsimonyInformative", nbPi);
}

