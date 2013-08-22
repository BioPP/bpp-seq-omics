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

SiteContainer* AbstractSpeciesSelectionMafStatistics::getSiteContainer_(const MafBlock& block)
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

AbstractSpeciesMultipleSelectionMafStatistics::AbstractSpeciesMultipleSelectionMafStatistics(const std::vector< std::vector<std::string> >& species):
  species_(species)
{
  size_t n = VectorTools::vectorUnion(species).size();
  size_t m = 0;
  for (size_t i =  0; i < species.size(); ++i)
    m += species_[i].size();
  if (m != n)
    throw Exception("AbstractSpeciesMultipleSelectionMafStatistics (constructor). Species selections must be fully distinct.");
}

vector<SiteContainer*> AbstractSpeciesMultipleSelectionMafStatistics::getSiteContainers_(const MafBlock& block)
{
  vector<SiteContainer*> alignments;
  for (size_t k = 0; k < species_.size(); ++k) {
    VectorSiteContainer* alignment = new VectorSiteContainer(block.getAlignment().getAlphabet());
    for (size_t i = 0; i < species_[k].size(); ++i) {
      if (block.hasSequenceForSpecies(species_[k][i])) {
        vector<const MafSequence*> selection = block.getSequencesForSpecies(species_[k][i]);
        for (size_t j = 0; j < selection.size(); ++j) {
          alignment->addSequence(*selection[j]);
        }
      }
    }
    alignments.push_back(alignment);
  }
  return alignments;
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
  bool hasOutgroup = (outgroup_ != "");
  bool isAnalyzable;
  auto_ptr<SiteContainer> alignment;
  const Sequence* outgroupSeq = 0;
  if (hasOutgroup) {
    isAnalyzable = (block.hasSequenceForSpecies(outgroup_) && block.getNumberOfSequences() > 1);
    if (isAnalyzable) {
      //We need to extract the outgroup sequence:
      outgroupSeq = &block.getSequenceForSpecies(outgroup_); //Here we assume there is only one! Otherwise we take the first one...
      alignment.reset(getSiteContainer_(block));
    }
  } else {
    isAnalyzable = (block.getNumberOfSequences() > 0);
    if (isAnalyzable)
      alignment.reset(getSiteContainer_(block));
  }
  if (isAnalyzable) {
    for (size_t i = 0; i < block.getNumberOfSites(); ++i) {
      //Note: we do not rely on SiteTool::getCounts as it would be unefficient to count everything.
      const Site& site = alignment->getSite(i);
      map<int, unsigned int> counts;
      bool isUnresolved = false;
      bool isSaturated = false;
      for (size_t j = 0; !isUnresolved && !isSaturated && j < site.size(); ++j) {
        state = site[j];
        if (alphabet_->isGap(state) || alphabet_->isUnresolved(state)) {
          isUnresolved = true;
        } else {
          counts[state]++;
          if (counts.size() > 2) {
            isSaturated = true;
          }
        }
      }
      if (isUnresolved) {
        nbUnresolved++;
      } else if (isSaturated) {
        nbSaturated++;
      } else if (hasOutgroup && (
          alignment->getAlphabet()->isGap((*outgroupSeq)[i]) ||
          alignment->getAlphabet()->isUnresolved((*outgroupSeq)[i]))) {
        nbUnresolved++;
      } else {
        //Determine frequency class:
        double count;
        if (counts.size() == 1) {
          if (hasOutgroup) {
            if (counts.begin()->first == (*outgroupSeq)[i])
              count = 0; //This is the ancestral state.
            else
              count = counts.begin()->second; //This is a derived state.
          } else {
            count = 0; //In this case we do not know, so we put 0.
          }
        } else {
          map<int, unsigned int>::iterator it = counts.begin();
          unsigned int count1 = it->second;
          it++;
          unsigned int count2 = it->second;
          if (hasOutgroup) {
            if (counts.begin()->first == (*outgroupSeq)[i])
              count = counts.rbegin()->second; //This is the ancestral state, therefore we other one is the derived state.
            else if (counts.rbegin()->first == (*outgroupSeq)[i])
              count = counts.begin()->second; //The second state is the ancestral one, therefore the first one is the derived state.
            else {
              //None of the two states are ancestral! The position is therefore discarded.
              isSaturated = true;
            }
          } else {
            count = min(count1, count2); //In this case we do not know, so we take the minimum of the two values.
          }
        }
        if (isSaturated)
          nbSaturated++;
        else {
          try {
            counts_[categorizer_.getCategory(count) - 1]++;
          } catch (OutOfRangeException& oof) {
            nbIgnored++;
          }
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

vector<string> FourSpeciesPatternCountsMafStatistics::getSupportedTags() const
{
  vector<string> tags;
  tags.push_back("f1100");
  tags.push_back("f0110");
  tags.push_back("f1010");
  tags.push_back("Ignored");
  return tags;
}

void FourSpeciesPatternCountsMafStatistics::compute(const MafBlock& block)
{
  counts_.assign(6, 0);
  auto_ptr<SiteContainer> alignment(getSiteContainer_(block));
  if (alignment->getNumberOfSequences() == 4) {
    unsigned int nbIgnored = 0;
    for (size_t i = 0; i < block.getNumberOfSites(); ++i) {
      const Site& site = alignment->getSite(i);
      if (SiteTools::isComplete(site)) {
        if (site[0] == site[1] &&
            site[2] != site[1] &&
            site[3] == site[2])
          counts_[0]++;
        else if (site[1] == site[2] &&
            site[1] != site[0] &&
            site[3] == site[0])
          counts_[1]++;
        else if (site[0] == site[2] &&
            site[1] != site[0] &&
            site[3] == site[1])
          counts_[2]++;
      } else {
        nbIgnored++;
      }
    }
    result_.setValue("f1100", counts_[0]);
    result_.setValue("f0110", counts_[1]);
    result_.setValue("f1010", counts_[2]);
    result_.setValue("Ignored", nbIgnored);
  } else {
    result_.setValue("f1100", 0);
    result_.setValue("f0110", 0);
    result_.setValue("f1010", 0);
    result_.setValue("Ignored", static_cast<double>(block.getNumberOfSites()));
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
  auto_ptr<SiteContainer> alignment(getSiteContainer_(block));
  unsigned int nbNg = 0;
  unsigned int nbCo = 0;
  unsigned int nbPi = 0;
  if (alignment->getNumberOfSequences() > 0) {
    for (size_t i = 0; i < alignment->getNumberOfSites(); ++i) {
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

vector<string> PolymorphismMafStatistics::getSupportedTags() const
{
  vector<string> tags;
  tags.push_back("F");
  tags.push_back("P");
  tags.push_back("FP");
  tags.push_back("PF");
  tags.push_back("FF");
  tags.push_back("X");
  tags.push_back("FX");
  tags.push_back("PX");
  tags.push_back("XF");
  tags.push_back("XP");
  return tags;
}

vector<int> PolymorphismMafStatistics::getPatterns_(const SiteContainer& sites)
{
  vector<int> patterns(sites.getNumberOfSites());
  for (size_t i = 0; i < sites.getNumberOfSites(); ++i) {
    const Site& site = sites.getSite(i);
    int s = -1; //Unresolved
    if (SiteTools::isComplete(site)) {
      if (SiteTools::isConstant(site)) {
        s = site[0]; //The fixed state
      } else {
        s = -10; //Polymorphic.
      }
    }
    patterns[i] = s;
  }
  return patterns;
}

void PolymorphismMafStatistics::compute(const MafBlock& block)
{
  vector<SiteContainer*> alignments(getSiteContainers_(block));
  unsigned int nbF = 0;
  unsigned int nbP = 0;
  unsigned int nbFF = 0;
  unsigned int nbFP = 0;
  unsigned int nbPF = 0;
  unsigned int nbX = 0;
  unsigned int nbFX = 0;
  unsigned int nbPX = 0;
  unsigned int nbXF = 0;
  unsigned int nbXP = 0;
  //Get all patterns:
  vector<int> patterns1(block.getNumberOfSites(), -1);
  vector<int> patterns2(block.getNumberOfSites(), -1); 
  if (alignments[0]->getNumberOfSequences() > 0) {
    patterns1 = getPatterns_(*alignments[0]);
  }
  if (alignments[1]->getNumberOfSequences() > 0) {
    patterns2 = getPatterns_(*alignments[1]);
  }
  //Compare patterns:
  for (size_t i = 0; i < block.getNumberOfSites(); ++i) {
    int p1 = patterns1[i];  
    int p2 = patterns2[i];
    switch (p1) {
      case -1 :
        switch (p2) {
          case -1 :
            nbX++;
            break;
          case -10 :
            nbXP++;
            break;
          default :
            nbXF++;
        }
        break;            
      
      case -10 :
        switch (p2) {
          case -1 :
            nbPX++;
            break;
          case -10 :
            nbP++;
            break;
          default :
            nbPF++;
        }
        break;            

      default :
        switch (p2) {
          case -1 :
            nbFX++;
            break;
          case -10 :
            nbFP++;
            break;
          default :
            if (p1 == p2)
              nbF++;
            else
              nbFF++;
        }
    }
  }

  //Set results:
  result_.setValue("F", nbF);
  result_.setValue("P", nbP);
  result_.setValue("FF", nbFF);
  result_.setValue("FP", nbFP);
  result_.setValue("PF", nbPF);
  result_.setValue("X", nbX);
  result_.setValue("FX", nbFX);
  result_.setValue("PX", nbPX);
  result_.setValue("XF", nbXF);
  result_.setValue("XP", nbXP);
}

vector<string> SequenceDiversityMafStatistics::getSupportedTags() const
{
  vector<string> tags;
  tags.push_back("NbSeggregating");
  tags.push_back("WattersonTheta");
  return tags;
}

void SequenceDiversityMafStatistics::compute(const MafBlock& block)
{
  auto_ptr<SiteContainer> alignment(getSiteContainer_(block));
  unsigned int nbSeg = 0;
  unsigned int nbTot = 0;
  if (alignment->getNumberOfSequences() > 0) {
    for (size_t i = 0; i < alignment->getNumberOfSites(); ++i) {
      const Site& site = alignment->getSite(i);
      if (SiteTools::isComplete(site)) {
        nbTot++;
        if (!SiteTools::isConstant(site))
          nbSeg++;
      }
    }
  }
  double wt = 0;
  if (nbSeg > 0) {
    size_t n = alignment->getNumberOfSequences();
    double hf = 0;
    for (double i = 1; i < n; ++i)
      hf += 1. / i;
    wt = static_cast<double>(nbSeg) / (static_cast<double>(nbTot) * hf);
  }
  result_.setValue("NbSeggregating", nbSeg);
  result_.setValue("WattersonTheta", wt);
}

