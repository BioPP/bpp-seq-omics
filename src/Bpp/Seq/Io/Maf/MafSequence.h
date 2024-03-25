// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _MAFSEQUENCE_H_
#define _MAFSEQUENCE_H_

#include "../../Feature/SequenceFeature.h"

#include <Bpp/Seq/SequenceWithAnnotation.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/SequenceTools.h>

namespace bpp
{
/**
 * @brief A sequence class which is used to store data from MAF files.
 *
 * It extends the SequenceWithAnnotation class to store MAF-specific features,
 * like the chromosome position. The sequence is its own listener,
 * and recomputes its "genomic" site by using the SequenceTools::getNumberOfSites
 * function when a content modification is performed.
 * Tags like begin and stop, however, have to be set by hand.
 *
 * A MAF sequence is necessarily a DNA sequence.
 */
class MafSequence :
  public SequenceWithAnnotation
{
private:
  bool hasCoordinates_;
  size_t begin_;
  std::string species_;
  std::string chromosome_;
  char strand_;
  size_t size_;
  size_t srcSize_;

public:
  MafSequence(std::shared_ptr<const Alphabet> alphabet = AlphabetTools::DNA_ALPHABET) :
    AbstractTemplateSymbolList<int>(alphabet),
    SequenceWithAnnotation(alphabet),
    hasCoordinates_(false),
    begin_(0),
    species_(""),
    chromosome_(""),
    strand_(0),
    size_(0),
    srcSize_(0)
  {}

  MafSequence(
      const std::string& name,
      const std::string& sequence,
      bool parseName = true,
      std::shared_ptr<const Alphabet> alphabet = AlphabetTools::DNA_ALPHABET) :
    AbstractTemplateSymbolList<int>(alphabet),
    SequenceWithAnnotation(name, sequence, alphabet),
    hasCoordinates_(false),
    begin_(0),
    species_(""),
    chromosome_(""),
    strand_(0),
    size_(0),
    srcSize_(0)
  {
    size_ = SequenceTools::getNumberOfSites(*this);
    if (parseName)
      splitNameIntoSpeciesAndChromosome(name, species_, chromosome_);
  }

  MafSequence(
      const std::string& name,
      const std::string& sequence,
      size_t begin,
      char strand,
      size_t srcSize,
      bool parseName = true,
      std::shared_ptr<const Alphabet> alphabet = AlphabetTools::DNA_ALPHABET) :
    AbstractTemplateSymbolList<int>(alphabet),
    SequenceWithAnnotation(name, sequence, alphabet),
    hasCoordinates_(true),
    begin_(begin),
    species_(""),
    chromosome_(""),
    strand_(strand),
    size_(0),
    srcSize_(srcSize)
  {
    size_ = SequenceTools::getNumberOfSites(*this);
    if (parseName)
      splitNameIntoSpeciesAndChromosome(name, species_, chromosome_);
  }

  MafSequence(const MafSequence& mafSeq) :
    AbstractTemplateSymbolList<int>(mafSeq),
    SequenceWithAnnotation(mafSeq),
    hasCoordinates_(mafSeq.hasCoordinates_),
    begin_(mafSeq.begin_),
    species_(mafSeq.species_),
    chromosome_(mafSeq.chromosome_),
    strand_(mafSeq.strand_),
    size_(mafSeq.size_),
    srcSize_(mafSeq.srcSize_)
  {}

  MafSequence& operator=(const MafSequence& mafSeq)
  {
    SequenceWithAnnotation::operator=(mafSeq);
    hasCoordinates_ = mafSeq.hasCoordinates_;
    begin_          = mafSeq.begin_;
    species_        = mafSeq.species_;
    chromosome_     = mafSeq.chromosome_;
    strand_         = mafSeq.strand_;
    size_           = mafSeq.size_;
    srcSize_        = mafSeq.srcSize_;
    return *this;
  }

  MafSequence(const SequenceInterface& seq, bool parseName = true) :
    AbstractTemplateSymbolList<int>(seq),
    SequenceWithAnnotation(seq),
    hasCoordinates_(false),
    begin_(0),
    species_(""),
    chromosome_(""),
    strand_(0),
    size_(0),
    srcSize_(0)
  {
    size_ = SequenceTools::getNumberOfSites(*this);
    if (parseName)
      splitNameIntoSpeciesAndChromosome(seq.getName(), species_, chromosome_);
  }

  MafSequence* clone() const override
  {
    return new MafSequence(*this);
  }

  MafSequence* cloneMeta() const
  {
    return new MafSequence(getName(), std::string(), begin_, strand_, srcSize_, true);
  }

  virtual ~MafSequence() {}

public:
  bool hasCoordinates() const { return hasCoordinates_; }

  void removeCoordinates() { hasCoordinates_ = false; begin_ = 0; }

  size_t start() const
  {
    if (hasCoordinates_) return begin_;
    else throw Exception("MafSequence::start(). Sequence " + getName() + " does not have coordinates.");
  }

  size_t stop() const
  {
    if (hasCoordinates_) return begin_ + size_;
    else throw Exception("MafSequence::stop(). Sequence " + getName() + " does not have coordinates.");
  }

  /**
   * @return A range with coordinates from this sequence.
   * @param origin Tell if coordinates according to original sequence should be used.
   * If 'yes' and the sequence is on the negative strand, the returned range will be computed as [SrcSize-Stop, SrcSize-Start[
   */
  Range<size_t> getRange(bool origin = true) const
  {
    if (hasCoordinates_)
    {
      if (origin && strand_ == '-')
      {
        return Range<size_t>(srcSize_ - stop(), srcSize_ - start());
      }
      else
      {
        return Range<size_t>(start(), stop());
      }
    }
    else throw Exception("MafSequence::getRange(). Sequence " + getName() + " does not have coordinates.");
  }

  void setName(const std::string& name) override
  {
    try
    {
      splitNameIntoSpeciesAndChromosome(name, species_, chromosome_);
    }
    catch (Exception& e)
    {
      species_ = "";
      chromosome_ = "";
    }
    SequenceWithAnnotation::setName(name);
  }

  static void splitNameIntoSpeciesAndChromosome(const std::string& name, std::string& species, std::string& chr)
  {
    size_t pos = name.find(".");
    if (pos != std::string::npos)
    {
      chr     = name.substr(pos + 1);
      species = name.substr(0, pos);
    }
    else
    {
      throw Exception("MafSequence::splitNameIntoSpeciesAndChromosome(). Invalid sequence name: " + name);
    }
  }

  const std::string& getSpecies() const { return species_; }

  const std::string& getChromosome() const { return chromosome_; }

  char getStrand() const { return strand_; }

  size_t getGenomicSize() const { return size_; }

  size_t getSrcSize() const { return srcSize_; }

  void setStart(size_t begin) { begin_ = begin; hasCoordinates_ = true; }

  void setChromosome(const std::string& chr)
  {
    chromosome_ = chr;
    SequenceWithAnnotation::setName(species_ + "." + chromosome_);
  }

  void setSpecies(const std::string& species)
  {
    species_ = species;
    SequenceWithAnnotation::setName(species_ + "." + chromosome_);
  }

  void setStrand(char s) { strand_ = s; }

  void setSrcSize(size_t srcSize) { srcSize_ = srcSize; }

  std::string getDescription() const { return getName() + strand_ + ":" + (hasCoordinates_ ? TextTools::toString(start()) + "-" + TextTools::toString(stop()) : "?-?"); }

  /**
   * @brief Extract a sub-sequence.
   *
   * @return A subsequence.
   * @param startAt Beginning of sub-sequence.
   * @param length  the length of the sub-sequence.
   */
  std::unique_ptr<MafSequence> subSequence(size_t startAt, size_t length) const;

private:
  void beforeSequenceChanged(const IntSymbolListEditionEvent& event) override {}
  void afterSequenceChanged(const IntSymbolListEditionEvent& event) override
  {
    size_ = SequenceTools::getNumberOfSites(*this);
  }
  void beforeSequenceInserted(const IntSymbolListInsertionEvent& event) override {}
  void afterSequenceInserted(const IntSymbolListInsertionEvent& event) override
  {
    size_ = SequenceTools::getNumberOfSites(*this);
  }
  void beforeSequenceDeleted(const IntSymbolListDeletionEvent& event) override {}
  void afterSequenceDeleted(const IntSymbolListDeletionEvent& event) override
  {
    size_ = SequenceTools::getNumberOfSites(*this);
  }
  void beforeSequenceSubstituted(const IntSymbolListSubstitutionEvent& event) override {}
  void afterSequenceSubstituted(const IntSymbolListSubstitutionEvent& event) override {}
};
} // end of namespace bpp.

#endif // _MAFSEQUENCE_H_
