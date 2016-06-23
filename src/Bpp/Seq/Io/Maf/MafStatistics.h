//
// File: MafStatistics.h
// Authors: Julien Dutheil
// Created: Mon Jun 25 2012
//

/*
Copyright or © or Copr. Bio++ Development Team, (2012)

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

#ifndef _MAFSTATISTICS_H_
#define _MAFSTATISTICS_H_

#include "MafBlock.h"

//From bpp-core:
#include <Bpp/Utils/MapTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Number.h>

//From the STL:
#include <map>
#include <string>

namespace bpp {

/**
 * @brief General interface for storing statistical results.
 *
 * @author Julien Dutheil
 * @see MafStatistics
 */
class MafStatisticsResult
{
  protected:
    mutable std::map<std::string, BppNumberI*> values_;

  public:
    MafStatisticsResult(): values_() {}
    virtual ~MafStatisticsResult() {}

    MafStatisticsResult(const MafStatisticsResult& msr): values_()
    {
      for (std::map<std::string, BppNumberI*>::const_iterator it = msr.values_.begin();
          it != msr.values_.end();
          ++it) {
        values_[it->first] = it->second->clone();
      }
    }

  public:
    virtual const BppNumberI& getValue(const std::string& tag) const throw (Exception) {
      std::map<std::string, BppNumberI*>::iterator it = values_.find(tag);
      if (it != values_.end())
        return *it->second;
      else
        throw Exception("MafStatisticsResult::getValue(). No value found for tag: " + tag + ".");
    }

    /**
     * @brief Associate a value to a certain tag. Any existing tag will be overwritten
     *
     * @param tag The name of the value to associate.
     * @param value The value to associate to the tag.
     */
    virtual void setValue(const std::string& tag, double value) throw (Exception) {
      if (values_[tag]) {
        delete values_[tag];
      }
      values_[tag] = new BppDouble(value);
    }

    /**
     * @brief Associate a value to a certain tag. Any existing tag will be overwritten
     *
     * @param tag The name of the value to associate.
     * @param value The value to associate to the tag.
     */
    virtual void setValue(const std::string& tag, int value) throw (Exception) {
      if (values_[tag]) {
        delete values_[tag];
      }
      values_[tag] = new BppInteger(value);
    }

    /**
     * @brief Associate a value to a certain tag. Any existing tag will be overwritten
     *
     * @param tag The name of the value to associate.
     * @param value The value to associate to the tag.
     */
    virtual void setValue(const std::string& tag, unsigned int value) throw (Exception) {
      if (values_[tag]) {
        delete values_[tag];
      }
      values_[tag] = new BppUnsignedInteger(value);
    }

    /**
     * @return A boolean saying whether a value is available for the given tag.
     * @param tag The name of the value to associate.
     */
    virtual bool hasValue(const std::string& tag) const {
      return (values_.find(tag) != values_.end()); 
    }

    /**
     * @return A vector with all available tags.
     */
    std::vector<std::string> getAvailableTags() const { return MapTools::getKeys(values_); }
};

/**
 * @brief A simple maf statistics result, with only one value.
 */
class SimpleMafStatisticsResult:
  public virtual MafStatisticsResult
{
  private:
    std::string name_;

  public:
    SimpleMafStatisticsResult(const std::string& name): MafStatisticsResult(), name_(name) {
      setValue(name, 0);
    }
    virtual ~SimpleMafStatisticsResult() {}

  public:
    virtual const BppNumberI& getValue(const std::string& tag) const throw (Exception) { return MafStatisticsResult::getValue(tag); }
    
    virtual const BppNumberI& getValue() const { return *values_[name_]; }

    virtual void setValue(const std::string& tag, double value) throw (Exception) {
      if (tag == name_)
        setValue(value);
      else
        throw Exception("SimpleMafStatisticsResult::setValue(). Unvalid tag name: " + tag + ".");
    }
 
    virtual void setValue(const std::string& tag, int value) throw (Exception) {
      if (tag == name_)
        setValue(value);
      else
        throw Exception("SimpleMafStatisticsResult::setValue(). Unvalid tag name: " + tag + ".");
    }
 
    virtual void setValue(const std::string& tag, unsigned int value) throw (Exception) {
      if (tag == name_)
        setValue(value);
      else
        throw Exception("SimpleMafStatisticsResult::setValue(). Unvalid tag name: " + tag + ".");
    }
   
    virtual void setValue(double value) {
      if (values_[name_]) delete values_[name_];
      values_[name_] = new BppDouble(value);
    }

    virtual void setValue(int value) {
      if (values_[name_]) delete values_[name_];
      values_[name_] = new BppInteger(value);
    }

    virtual void setValue(unsigned int value) {
      if (values_[name_]) delete values_[name_];
      values_[name_] = new BppUnsignedInteger(value);
    }

};

/**
 * @brief General interface for computing statistics based on a Maf block.
 *
 * @author Julien Dutheil
 * @see MafBlock
 */
class MafStatistics
{
  public:
    MafStatistics() {}
    virtual ~MafStatistics() {}

  public:
    virtual std::string getShortName() const = 0;
    virtual std::string getFullName() const = 0;
    virtual const MafStatisticsResult& getResult() const = 0;
    virtual void compute(const MafBlock& block) = 0;

    /**
     * @return A vector with all available tags.
     */
    virtual std::vector<std::string> getSupportedTags() const = 0;

};

/**
 * @brief Partial implementation of MafStatistics, for convenience.
 */
class AbstractMafStatistics:
  public virtual MafStatistics
{
  protected:
    MafStatisticsResult result_;

  public:
    AbstractMafStatistics(): result_() {}
    virtual ~AbstractMafStatistics() {}

  public:
    const MafStatisticsResult& getResult() const { return result_; }
};

/**
 * @brief Partial implementation of MafStatistics, for convenience.
 */
class AbstractMafStatisticsSimple:
  public MafStatistics
{
  protected:
    SimpleMafStatisticsResult result_;

  public:
    AbstractMafStatisticsSimple(const std::string& name): result_(name) {}
    virtual ~AbstractMafStatisticsSimple() {}

  public:
    const SimpleMafStatisticsResult& getResult() const { return result_; }
    std::vector<std::string> getSupportedTags() const { return result_.getAvailableTags(); }
};

/**
 * @brief Computes the pairwise divergence for a pair of sequences in a maf block.
 */
class PairwiseDivergenceMafStatistics:
  public AbstractMafStatisticsSimple
{
  private:
    std::string species1_;
    std::string species2_;

  public:
    PairwiseDivergenceMafStatistics(const std::string& species1, const std::string& species2):
      AbstractMafStatisticsSimple("Divergence"), species1_(species1), species2_(species2) {}

    ~PairwiseDivergenceMafStatistics() {}

  public:
    std::string getShortName() const { return "Div." + species1_ + "-" + species2_; }
    std::string getFullName() const { return "Pairwise divergence between " + species1_ + " and " + species2_ + "."; }
    void compute(const MafBlock& block);

};

/**
 * @brief Computes the number of sequences in a maf block.
 */
class BlockSizeMafStatistics:
  public AbstractMafStatisticsSimple
{
  public:
    BlockSizeMafStatistics(): AbstractMafStatisticsSimple("BlockSize") {}
    ~BlockSizeMafStatistics() {}

  public:
    std::string getShortName() const { return "BlockSize"; }
    std::string getFullName() const { return "Number of sequences."; }
    void compute(const MafBlock& block) {
      result_.setValue(static_cast<double>(block.getNumberOfSequences()));
    }
};

/**
 * @brief Computes the number of columns in a maf block.
 */
class BlockLengthMafStatistics:
  public AbstractMafStatisticsSimple
{
  public:
    BlockLengthMafStatistics(): AbstractMafStatisticsSimple("BlockLength") {}
    ~BlockLengthMafStatistics() {}

  public:
    std::string getShortName() const { return "BlockLength"; }
    std::string getFullName() const { return "Number of sites."; }
    void compute(const MafBlock& block) {
      result_.setValue(static_cast<double>(block.getNumberOfSites()));
    }
};

/**
 * @brief Retrieve the sequence length (number of nucleotides) for a given species in a maf block.
 *
 * If no sequence is found for the current block, 0 is returned.
 * If several sequences are found for a given species, an exception is thrown.
 */
class SequenceLengthMafStatistics:
  public AbstractMafStatisticsSimple
{
  private:
    std::string species_;

  public:
    SequenceLengthMafStatistics(const std::string& species): AbstractMafStatisticsSimple("BlockSize"), species_(species) {}
    ~SequenceLengthMafStatistics() {}

  public:
    std::string getShortName() const { return "SequenceLengthFor" + species_; }
    std::string getFullName() const { return "Sequence length for species " + species_; }
    void compute(const MafBlock& block) {
      std::vector<const MafSequence*> seqs = block.getSequencesForSpecies(species_);
      if (seqs.size() == 0)
        result_.setValue(0.);
      else if (seqs.size() == 1)
        result_.setValue(static_cast<double>(SequenceTools::getNumberOfSites(*seqs[0])));
      else
        throw Exception("SequenceLengthMafStatistics::compute. More than one sequence found for species " + species_ + " in current block.");
    }
};


/**
 * @brief Retrieves the alignment score of a maf block.
 */
class AlignmentScoreMafStatistics:
  public AbstractMafStatisticsSimple
{
  public:
    AlignmentScoreMafStatistics(): AbstractMafStatisticsSimple("AlnScore") {}
    ~AlignmentScoreMafStatistics() {}

  public:
    std::string getShortName() const { return "AlnScore"; }
    std::string getFullName() const { return "Alignment score."; }
    void compute(const MafBlock& block) {
      result_.setValue(block.getScore());
    }
};



/**
 * @brief Partial implementation of MafStatistics for method working on a subset of species, in a site-wise manner.
 *
 * This class stores a selection of species and create for each block the corresponding SiteContainer instance.
 */
class AbstractSpeciesSelectionMafStatistics:
  public virtual MafStatistics
{
  private:
    std::vector<std::string> species_;
    bool noSpeciesMeansAllSpecies_;

  protected:
    std::string suffix_;

  public:
    AbstractSpeciesSelectionMafStatistics(
        const std::vector<std::string>& species,
        bool noSpeciesMeansAllSpecies = false,
        const std::string& suffix = ""):
      species_(species),
      noSpeciesMeansAllSpecies_(noSpeciesMeansAllSpecies),
      suffix_(suffix)
    {}

  protected:
    SiteContainer* getSiteContainer_(const MafBlock& block);

};


/**
 * @brief Partial implementation of MafStatistics for method working on multiple distinct subsets of species, in a site-wise manner.
 *
 * This class stores two non-overlapping selections of species and create for each block the corresponding SiteContainer instances.
 */
class AbstractSpeciesMultipleSelectionMafStatistics:
  public virtual MafStatistics
{
  private:
    std::vector< std::vector<std::string> > species_;

  public:
    AbstractSpeciesMultipleSelectionMafStatistics(const std::vector< std::vector<std::string> >& species);

  protected:
    std::vector<SiteContainer*> getSiteContainers_(const MafBlock& block);

};



/**
 * @brief Compute the base frequencies of a maf block.
 *
 * For each block, provides the following numbers (with their corresponding tags):
 * - A: total counts of A
 * - C: total counts of C
 * - G: total counts of G
 * - T [or U]: total counts of T/U
 * - Gap: total counts of gaps
 * - Unresolved: total counts of unresolved characters
 * The sum of all characters should equal BlockSize x BlockLength 
 */
class CharacterCountsMafStatistics:
  public AbstractMafStatistics,
  public AbstractSpeciesSelectionMafStatistics
{
  private:
    const Alphabet* alphabet_;

  public:
    CharacterCountsMafStatistics(const Alphabet* alphabet, const std::vector<std::string>& species, const std::string suffix):
      AbstractMafStatistics(),
      AbstractSpeciesSelectionMafStatistics(species, true, suffix),
      alphabet_(alphabet) {}

    CharacterCountsMafStatistics(const CharacterCountsMafStatistics& stats):
      AbstractMafStatistics(stats),
      AbstractSpeciesSelectionMafStatistics(stats),
      alphabet_(stats.alphabet_) {}
    
    CharacterCountsMafStatistics& operator=(const CharacterCountsMafStatistics& stats) {
      AbstractMafStatistics::operator=(stats);
      AbstractSpeciesSelectionMafStatistics::operator=(stats);
      alphabet_ = stats.alphabet_;
      return *this;
    }

    virtual ~CharacterCountsMafStatistics() {}

  public:
    std::string getShortName() const { return "Counts" + suffix_; }
    std::string getFullName() const { return "Character counts (" + suffix_ + ")."; }
    void compute(const MafBlock& block);
    std::vector<std::string> getSupportedTags() const;
};



/**
 * @brief Compute the Site Frequency Spectrum of a maf block.
 *
 * If no outgroup is provided, the ancestral states are considered as unknown
 * and the unfolded spectrum is computed, so that 10000 and 11110 sites are treated equally.
 */
class SiteFrequencySpectrumMafStatistics:
  public AbstractMafStatistics,
  public AbstractSpeciesSelectionMafStatistics
{
  private:
    class Categorizer {
      private:
        std::vector<double> bounds_;

      public:
        Categorizer(const std::vector<double>& bounds):
          bounds_(bounds) {
          std::sort(bounds_.begin(), bounds_.end());
        }

      public:
        size_t getNumberOfCategories() const { return (bounds_.size() - 1); }

        //Category numbers start at 1!
        size_t getCategory(double value) const throw (OutOfRangeException) {
          if (value < bounds_[0])
            throw OutOfRangeException("SiteFrequencySpectrumMafStatistics::Categorizer::getCategory.", value, *bounds_.begin(), *bounds_.rbegin());
          for (size_t i = 1; i < bounds_.size(); ++i) {
            if (value < bounds_[i])
              return i;
          }
          throw OutOfRangeException("SiteFrequencySpectrumMafStatistics::Categorizer::getCategory.", value, *bounds_.begin(), *bounds_.rbegin());
        }
    };

  private:
    const Alphabet* alphabet_;
    Categorizer categorizer_;
    std::vector<unsigned int> counts_;
    std::string outgroup_;

  public:
    SiteFrequencySpectrumMafStatistics(const Alphabet* alphabet, const std::vector<double>& bounds, const std::vector<std::string>& ingroup, const std::string outgroup = ""):
      AbstractMafStatistics(),
      AbstractSpeciesSelectionMafStatistics(ingroup),
      alphabet_(alphabet),
      categorizer_(bounds),
      counts_(bounds.size() - 1),
      outgroup_(outgroup)
    {}

    SiteFrequencySpectrumMafStatistics(const SiteFrequencySpectrumMafStatistics& stats):
      AbstractMafStatistics(),
      AbstractSpeciesSelectionMafStatistics(stats),
      alphabet_(stats.alphabet_),
      categorizer_(stats.categorizer_),
      counts_(stats.counts_),
      outgroup_(stats.outgroup_)
    {}

    SiteFrequencySpectrumMafStatistics& operator=(const SiteFrequencySpectrumMafStatistics& stats) {
      AbstractMafStatistics::operator=(stats);
      AbstractSpeciesSelectionMafStatistics::operator=(stats);
      alphabet_    = stats.alphabet_;
      categorizer_ = stats.categorizer_;
      counts_      = stats.counts_;
      outgroup_    = stats.outgroup_;
      return *this;
    }

    virtual ~SiteFrequencySpectrumMafStatistics() {}

  public:
    std::string getShortName() const { return "SiteFrequencySpectrum"; }
    std::string getFullName() const { return "Site frequency spectrum."; }
    void compute(const MafBlock& block);
    std::vector<std::string> getSupportedTags() const;
};


/**
 * @brief Compute the frequency of site patterns for a quadruplet of species.
 *
 * Only parsimony informative sites are categorized.
 * Species: A B C D
 * P1       1 1 0 0
 * P2       0 1 1 0
 * P3       1 0 1 0
 * Sites with more than two states are ignored, as well as sites containing gaps or unresolved characters.
 */
class FourSpeciesPatternCountsMafStatistics:
  public AbstractMafStatistics,
  public AbstractSpeciesSelectionMafStatistics
{
  private:
    const Alphabet* alphabet_;
    std::vector<unsigned int> counts_;

  public:
    FourSpeciesPatternCountsMafStatistics(
        const Alphabet* alphabet,
        const std::vector<std::string>& species):
      AbstractMafStatistics(),
      AbstractSpeciesSelectionMafStatistics(species),
      alphabet_(alphabet),
      counts_(6)
    {
      if (species.size() != 4)
        throw Exception("FourSpeciesPatternCountsMafStatistics, constructor: 4 species should be provided.");
      if (VectorTools::unique(species).size() != 4)
        throw Exception("FourSpeciesPatternCountsMafStatistics, constructor: duplicated species name!");
    }

    FourSpeciesPatternCountsMafStatistics(const FourSpeciesPatternCountsMafStatistics& stats):
      AbstractMafStatistics(),
      AbstractSpeciesSelectionMafStatistics(stats),
      alphabet_(stats.alphabet_),
      counts_(stats.counts_)
    {}

    FourSpeciesPatternCountsMafStatistics& operator=(const FourSpeciesPatternCountsMafStatistics& stats) {
      AbstractMafStatistics::operator=(stats);
      AbstractSpeciesSelectionMafStatistics::operator=(stats);
      alphabet_    = stats.alphabet_;
      counts_      = stats.counts_;
      return *this;
    }

    virtual ~FourSpeciesPatternCountsMafStatistics() {}

  public:
    std::string getShortName() const { return "FourSpeciesPatternCounts"; }
    std::string getFullName() const { return "FourSpecies pattern counts."; }
    void compute(const MafBlock& block);
    std::vector<std::string> getSupportedTags() const;
};



/**
 * @brief Compute a few site statistics in a maf block.
 *
 * Computed statistics include:
 * - Number of sites without gaps
 * - Number of complete sites (no gap, no unresolved)
 * - Number of complete sites with only one state (constant sites)
 * - Number of complete biallelic sites
 * - Number of complete triallelic sites
 * - Number of complete Quadriallelic sites
 * - Number of parsimony informative sites
 */
class SiteMafStatistics:
  public AbstractMafStatistics,
  public AbstractSpeciesSelectionMafStatistics
{
  public:
    SiteMafStatistics(const std::vector<std::string>& species):
      AbstractMafStatistics(),
      AbstractSpeciesSelectionMafStatistics(species)
    {}

    virtual ~SiteMafStatistics() {}

  public:
    std::string getShortName() const { return "SiteStatistics"; }
    std::string getFullName() const { return "Site statistics."; }
    void compute(const MafBlock& block);
    std::vector<std::string> getSupportedTags() const;
};


/**
 * @brief Counts number of polymorphic / fixed sites in two populations.
 *
 * The two populations are defined as two distinct sets of species.
 * The following counts are computed and returned:
 * - P: number of sites polymorphic in both populations
 * - F: number of sites fixed in both populations
 * - FF: number of sites fixed in both populations, but with distinct states
 * - PF / FP: number of sites polymorphic in one species and fixed in the other.
 * - X: unresolved (because of gap or generic character)
 * - FX / PX / XF / XP: unresolved in one population
 */
class PolymorphismMafStatistics:
  public AbstractMafStatistics,
  public AbstractSpeciesMultipleSelectionMafStatistics
{
  public:
    PolymorphismMafStatistics(const std::vector< std::vector<std::string> >& species):
      AbstractMafStatistics(),
      AbstractSpeciesMultipleSelectionMafStatistics(species)
    {
      if (species.size() != 2)
        throw Exception("PolymorphismStatistics: exactly two species selection should be provided.");
    }

    virtual ~PolymorphismMafStatistics() {}

  public:
    std::string getShortName() const { return "PolymorphismStatistics"; }
    std::string getFullName() const { return "Polymorphism statistics."; }
    void compute(const MafBlock& block);
    std::vector<std::string> getSupportedTags() const;

  private:
    static std::vector<int> getPatterns_(const SiteContainer& sites); 
};



/**
 * @brief Provide estimates of sequence diversity.
 *
 * - Number of segregating sites
 * - Watterson's theta (per site)
 * - Tajima's pi (average pairwise differences)
 * - Tajima's D
 *
 * Only fully resolved sites are analyzed (no gap, no generic character).
 */
class SequenceDiversityMafStatistics:
  public AbstractMafStatistics,
  public AbstractSpeciesSelectionMafStatistics
{
  public:
    SequenceDiversityMafStatistics(const std::vector<std::string>& ingroup):
      AbstractMafStatistics(),
      AbstractSpeciesSelectionMafStatistics(ingroup)
    {}

    virtual ~SequenceDiversityMafStatistics() {}

  public:
    std::string getShortName() const { return "SequenceDiversityStatistics"; }
    std::string getFullName() const { return "Sequence diversity statistics."; }
    void compute(const MafBlock& block);
    std::vector<std::string> getSupportedTags() const;

  private:
    static std::vector<int> getPatterns_(const SiteContainer& sites); 
};


} // end of namespace bpp

#endif //_MAFSTATISTICS_H_

