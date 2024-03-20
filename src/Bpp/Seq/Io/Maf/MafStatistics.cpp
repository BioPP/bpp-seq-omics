// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "MafStatistics.h"
#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SiteTools.h>

// From bpp-core:
#include <Bpp/Numeric/NumConstants.h>

// From the STL:
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

unique_ptr<SiteContainerInterface> AbstractSpeciesSelectionMafStatistics::getSiteContainer_(const MafBlock& block)
{
  auto alignment = make_unique<VectorSiteContainer>(block.getAlphabet());
  if (noSpeciesMeansAllSpecies_ && species_.size() == 0)
  {
    for (size_t i = 0; i < block.getNumberOfSequences(); ++i)
    {
      auto tmpSeq = make_unique<Sequence>(block.sequence(i));
      alignment->addSequence(tmpSeq->getName(), tmpSeq);
    }
  }
  // Otherwise, we select species:
  for (size_t i = 0; i < species_.size(); ++i)
  {
    if (block.hasSequenceForSpecies(species_[i]))
    {
      vector<const MafSequence*> selection = block.getSequencesForSpecies(species_[i]);
      for (size_t j = 0; j < selection.size(); ++j)
      {
  auto tmpSeq = make_unique<Sequence>(*selection[j]);
        alignment->addSequence(tmpSeq->getName(), tmpSeq);
      }
    }
  }
  return alignment;
}

AbstractSpeciesMultipleSelectionMafStatistics::AbstractSpeciesMultipleSelectionMafStatistics(const std::vector< std::vector<std::string> >& species) :
  species_(species)
{
  size_t n = VectorTools::vectorUnion(species).size();
  size_t m = 0;
  for (size_t i =  0; i < species.size(); ++i)
  {
    m += species_[i].size();
  }
  if (m != n)
    throw Exception("AbstractSpeciesMultipleSelectionMafStatistics (constructor). Species selections must be fully distinct.");
}

vector<unique_ptr<SiteContainerInterface>> AbstractSpeciesMultipleSelectionMafStatistics::getSiteContainers_(const MafBlock& block)
{
  vector<unique_ptr<SiteContainerInterface>> alignments;
  for (size_t k = 0; k < species_.size(); ++k)
  {
    auto alignment = make_unique<VectorSiteContainer>(block.getAlphabet());
    for (size_t i = 0; i < species_[k].size(); ++i)
    {
      if (block.hasSequenceForSpecies(species_[k][i]))
      {
        vector<const MafSequence*> selection = block.getSequencesForSpecies(species_[k][i]);
        for (size_t j = 0; j < selection.size(); ++j)
        {
    auto tmpSeq = make_unique<Sequence>(*selection[j]);
          alignment->addSequence(tmpSeq->getName(), tmpSeq);
        }
      }
    }
    alignments.push_back(move(alignment));
  }
  return alignments;
}

vector<string> CharacterCountsMafStatistics::getSupportedTags() const
{
  vector<string> tags;
  for (int i = 0; i < static_cast<int>(alphabet_->getSize()); ++i)
  {
    tags.push_back(alphabet_->intToChar(i));
  }
  tags.push_back("Gap");
  tags.push_back("Unresolved");

  return tags;
}

void CharacterCountsMafStatistics::compute(const MafBlock& block)
{
  std::map<int, unsigned int> counts;
  auto sites = getSiteContainer_(block);
  SequenceContainerTools::getCounts(*sites, counts);
  for (int i = 0; i < static_cast<int>(alphabet_->getSize()); ++i)
  {
    result_.setValue(alphabet_->intToChar(i), counts[i]);
  }
  result_.setValue("Gap", counts[alphabet_->getGapCharacterCode()]);
  double countUnres = 0;
  for (auto& it : counts)
  {
    if (alphabet_->isUnresolved(it.first))
      countUnres += it.second;
  }
  result_.setValue("Unresolved", countUnres);
}

vector<string> SiteFrequencySpectrumMafStatistics::getSupportedTags() const
{
  vector<string> tags;
  for (size_t i = 0; i < categorizer_.getNumberOfCategories(); ++i)
  {
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
  unique_ptr<SiteContainerInterface> alignment;
  const SequenceInterface* outgroupSeq = nullptr;
  if (hasOutgroup)
  {
    isAnalyzable = (block.hasSequenceForSpecies(outgroup_) && block.getNumberOfSequences() > 1);
    if (isAnalyzable)
    {
      // We need to extract the outgroup sequence:
      outgroupSeq = &block.sequenceForSpecies(outgroup_); // Here we assume there is only one! Otherwise we take the first one...
      alignment = getSiteContainer_(block);
    }
  }
  else
  {
    isAnalyzable = (block.getNumberOfSequences() > 0);
    if (isAnalyzable)
      alignment = getSiteContainer_(block);
  }
  if (isAnalyzable)
  {
    for (size_t i = 0; i < alignment->getNumberOfSites(); ++i)
    {
      // Note: we do not rely on SiteTool::getCounts as it would be unefficient to count everything.
      const Site& site = alignment->site(i);
      map<int, unsigned int> counts;
      bool isUnresolved = false;
      bool isSaturated = false;
      for (size_t j = 0; !isUnresolved && !isSaturated && j < site.size(); ++j)
      {
        state = site[j];
        if (alphabet_->isGap(state) || alphabet_->isUnresolved(state))
        {
          isUnresolved = true;
        }
        else
        {
          counts[state]++;
          if (counts.size() > 2)
          {
            isSaturated = true;
          }
        }
      }
      if (isUnresolved)
      {
        nbUnresolved++;
      }
      else if (isSaturated)
      {
        nbSaturated++;
      }
      else if (hasOutgroup && (
                 alignment->getAlphabet()->isGap((*outgroupSeq)[i]) ||
                 alignment->getAlphabet()->isUnresolved((*outgroupSeq)[i])))
      {
        nbUnresolved++;
      }
      else
      {
        // Determine frequency class:
        double count;
        if (counts.size() == 1)
        {
          if (hasOutgroup)
          {
            if (counts.begin()->first == (*outgroupSeq)[i])
              count = 0; // This is the ancestral state.
            else
              count = counts.begin()->second; // This is a derived state.
          }
          else
          {
            count = 0; // In this case we do not know, so we put 0.
          }
        }
        else
        {
          map<int, unsigned int>::iterator it = counts.begin();
          unsigned int count1 = it->second;
          it++;
          unsigned int count2 = it->second;
          if (hasOutgroup)
          {
            if (counts.begin()->first == (*outgroupSeq)[i])
              count = counts.rbegin()->second; // This is the ancestral state, therefore we other one is the derived state.
            else if (counts.rbegin()->first == (*outgroupSeq)[i])
              count = counts.begin()->second; // The second state is the ancestral one, therefore the first one is the derived state.
            else
            {
              // None of the two states are ancestral! The position is therefore discarded.
              isSaturated = true;
            }
          }
          else
          {
            count = min(count1, count2); // In this case we do not know, so we take the minimum of the two values.
          }
        }
        if (isSaturated)
          nbSaturated++;
        else
        {
          try
          {
            counts_[categorizer_.getCategory(count) - 1]++;
          }
          catch (OutOfRangeException& oof)
          {
            nbIgnored++;
          }
        }
      }
    }
  }
  result_.setValue("Unresolved", nbUnresolved);
  result_.setValue("Saturated", nbSaturated);
  result_.setValue("Ignored", nbIgnored);
  for (size_t i = 0; i < counts_.size(); ++i)
  {
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
  auto alignment = getSiteContainer_(block);
  if (alignment->getNumberOfSequences() == 4)
  {
    unsigned int nbIgnored = 0;
    for (size_t i = 0; i < block.getNumberOfSites(); ++i)
    {
      const Site& site = alignment->site(i);
      if (SiteTools::isComplete(site))
      {
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
      }
      else
      {
        nbIgnored++;
      }
    }
    result_.setValue("f1100", counts_[0]);
    result_.setValue("f0110", counts_[1]);
    result_.setValue("f1010", counts_[2]);
    result_.setValue("Ignored", nbIgnored);
  }
  else
  {
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
  tags.push_back("NbConstant");
  tags.push_back("NbBiallelic");
  tags.push_back("NbTriallelic");
  tags.push_back("NbQuadriallelic");
  tags.push_back("NbParsimonyInformative");
  return tags;
}

void SiteMafStatistics::compute(const MafBlock& block)
{
  auto alignment = getSiteContainer_(block);
  unsigned int nbNg = 0;
  unsigned int nbCo = 0;
  unsigned int nbPi = 0;
  unsigned int nbP1 = 0;
  unsigned int nbP2 = 0;
  unsigned int nbP3 = 0;
  unsigned int nbP4 = 0;
  if (alignment->getNumberOfSequences() > 0)
  {
    for (size_t i = 0; i < alignment->getNumberOfSites(); ++i)
    {
      if (!SiteTools::hasGap(alignment->site(i)))
        nbNg++;
      if (SiteTools::isComplete(alignment->site(i)))
      {
        nbCo++;
        map<int, size_t> counts;
        SiteTools::getCounts(alignment->site(i), counts);
        switch (counts.size())
        {
        case 1: nbP1++; break;
        case 2: nbP2++; break;
        case 3: nbP3++; break;
        case 4: nbP4++; break;
        default: throw Exception("The impossible happened. Probably a distortion in the Minkowski space.");
        }
      }
      if (SiteTools::isParsimonyInformativeSite(alignment->site(i)))
        nbPi++;
    }
  }
  result_.setValue("NbWithoutGap", nbNg);
  result_.setValue("NbComplete", nbCo);
  result_.setValue("NbConstant", nbP1);
  result_.setValue("NbBiallelic", nbP2);
  result_.setValue("NbTriallelic", nbP3);
  result_.setValue("NbQuadriallelic", nbP4);
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

vector<int> PolymorphismMafStatistics::getPatterns_(const SiteContainerInterface& sites)
{
  vector<int> patterns(sites.getNumberOfSites());
  for (size_t i = 0; i < sites.getNumberOfSites(); ++i)
  {
    const Site& site = sites.site(i);
    int s = -1; // Unresolved
    if (SiteTools::isComplete(site))
    {
      if (SiteTools::isConstant(site))
      {
        s = site[0]; // The fixed state
      }
      else
      {
        s = -10; // Polymorphic.
      }
    }
    patterns[i] = s;
  }
  return patterns;
}

void PolymorphismMafStatistics::compute(const MafBlock& block)
{
  vector<unique_ptr<SiteContainerInterface>> alignments = getSiteContainers_(block);
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
  // Get all patterns:
  vector<int> patterns1(block.getNumberOfSites(), -1);
  vector<int> patterns2(block.getNumberOfSites(), -1);
  if (alignments[0]->getNumberOfSequences() > 0)
  {
    patterns1 = getPatterns_(*alignments[0]);
  }
  if (alignments[1]->getNumberOfSequences() > 0)
  {
    patterns2 = getPatterns_(*alignments[1]);
  }
  // Compare patterns:
  for (size_t i = 0; i < block.getNumberOfSites(); ++i)
  {
    int p1 = patterns1[i];
    int p2 = patterns2[i];
    switch (p1)
    {
    case -1:
      switch (p2)
      {
      case -1:
        nbX++;
        break;
      case -10:
        nbXP++;
        break;
      default:
        nbXF++;
      }
      break;

    case -10:
      switch (p2)
      {
      case -1:
        nbPX++;
        break;
      case -10:
        nbP++;
        break;
      default:
        nbPF++;
      }
      break;

    default:
      switch (p2)
      {
      case -1:
        nbFX++;
        break;
      case -10:
        nbFP++;
        break;
      default:
        if (p1 == p2)
          nbF++;
        else
          nbFF++;
      }
    }
  }

  // Set results:
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
  tags.push_back("TajimaPi");
  tags.push_back("TajimaD");
  return tags;
}

void SequenceDiversityMafStatistics::compute(const MafBlock& block)
{
  unique_ptr<SiteContainerInterface> alignment = getSiteContainer_(block);
  // Get only complete sites:
  unique_ptr<VectorSiteContainer> alignment2 = SiteContainerTools::getCompleteSites<Site, Sequence>(*alignment);

  double S = 0;
  size_t nbTot = 0;
  size_t n = alignment2->getNumberOfSequences();
  if (n > 0)
  {
    for (size_t i = 0; i < alignment2->getNumberOfSites(); ++i)
    {
      const Site& site = alignment2->site(i);
      if (SiteTools::isComplete(site))
      {
        nbTot++;
        if (!SiteTools::isConstant(site))
          S++;
      }
    }
  }
  else
  {
    nbTot = alignment2->getNumberOfSites();
  }

  double a1 = 0;
  double a2 = 0;
  double dn = static_cast<double>(n);
  for (double i = 1; i < dn; ++i)
  {
    a1 += 1. / i;
    a2 += 1. / (i * i);
  }
  double wt = S / (static_cast<double>(nbTot) * a1);
  double b1 = (dn + 1) / (3 * (dn - 1));
  double b2 = 2 * (dn * dn + dn + 3) / (9 * dn * (dn - 1));
  double c1 = b1 - 1. / a1;
  double c2 = b2 - (dn + 2) / (a1 * dn) + a2 / (a1 * a1);
  double e1 = c1 / a1;
  double e2 = c2 / (a1 * a1 + a2);

  // Compute pairwise heterozigosity:
  double pi = 0;
  for (size_t i = 0; i < n - 1; ++i)
  {
    for (size_t j = i + 1; j < n; ++j)
    {
      pi += SiteContainerTools::computeSimilarity(
        alignment2->sequence(i),
        alignment2->sequence(j),
        true,
        SiteContainerTools::SIMILARITY_NOGAP,
        true);
    }
  }
  pi /= static_cast<double>((n - 1) * n / 2);

  // Compute Tajima's D:
  double tajd = static_cast<double>(nbTot) * (pi - wt) / sqrt(e1 * S + e2 * S * (S - 1));

  result_.setValue("NbSeggregating", S);
  result_.setValue("WattersonTheta", wt);
  result_.setValue("TajimaPi", pi);
  result_.setValue("TajimaD", tajd);
}
