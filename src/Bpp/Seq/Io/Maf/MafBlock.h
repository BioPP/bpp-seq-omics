// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _MAFBLOCK_H_
#define _MAFBLOCK_H_

#include "MafSequence.h"
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>

#include <Bpp/Clonable.h>

namespace bpp
{
/**
 * @brief A synteny block data structure, the basic unit of a MAF alignment file.
 *
 * This class basically contains a AlignedSequenceContainer made of MafSequence objects.
 */
class MafBlock :
  protected TemplateAlignedSequenceContainer<MafSequence, Site>
{
private:
  double score_;
  unsigned int pass_;
  std::map<std::string, std::unique_ptr<Clonable>> properties_;
  unsigned int idCounter_;

public:
  MafBlock() :
    TemplateAlignedSequenceContainer(AlphabetTools::DNA_ALPHABET),
    score_(log(0)),
    pass_(0),
    properties_(),
    idCounter_(0)
  {}

  MafBlock(const MafBlock& block) :
    TemplateAlignedSequenceContainer(block),
    score_(block.score_),
    pass_(block.pass_),
    properties_(),
    idCounter_(block.idCounter_)
  {
    for (const auto& it : block.properties_)
    {
      properties_[it.first].reset(it.second->clone());
    }
  }

  MafBlock& operator=(const MafBlock& block)
  {
    TemplateAlignedSequenceContainer::operator=(block),
    score_     = block.score_;
    pass_      = block.pass_;
    deleteProperties_();
    for (const auto& it : block.properties_)
    {
      properties_[it.first].reset(it.second->clone());
    }
    idCounter_ = block.idCounter_;
    return *this;
  }

  MafBlock* clone() const override { return new MafBlock(*this); }

  virtual ~MafBlock() { deleteProperties_(); }

public:
  void setScore(double score) { score_ = score; }
  void setPass(unsigned int pass) { pass_ = pass; }

  double getScore() const { return score_; }
  unsigned int getPass() const { return pass_; }

  /**
   * @brief Export the content of the block to an AlignedSequenceContainer object.
   *
   * All non-sequence information will not be passed to the new container.
   *
   * @return A pointer toward a new AlignedSequenceContainer object.
   */
  std::unique_ptr<AlignedSequenceContainer> getAlignment() const
  {
    auto aln = std::make_unique<AlignedSequenceContainer>(AlphabetTools::DNA_ALPHABET);
    for (size_t i = 0; i < getNumberOfSequences(); ++i)
    {
      auto copySeq = std::make_unique<Sequence>(sequence(i));
      aln->addSequence(sequenceKey(i), copySeq);
    }
    return aln;
  }

  /**
   * @brief Export a selection of the content of the block to an AlignedSequenceContainer object.
   * Only sequences for the given set of species will be copied in the new container.
   * All non-sequence information will not be passed to the new container.
   *
   * @return A pointer toward a new AlignedSequenceContainer object.
   */
  std::unique_ptr<AlignedSequenceContainer> getAlignment(const std::vector<std::string>& species) const
  {
    auto aln = std::make_unique<AlignedSequenceContainer>(AlphabetTools::DNA_ALPHABET);
    for (size_t i = 0; i < getNumberOfSequences(); ++i)
    {
      if (VectorTools::contains(species, sequence(i).getSpecies()))
      {
        auto copySeq = std::make_unique<Sequence>(sequence(i));
        aln->addSequence(sequenceKey(i), copySeq);
      }
    }
    return aln;
  }

  void addSequence(std::unique_ptr<MafSequence>& sequence) override
  {
    std::string key = "maf_seq_" + TextTools::toString(idCounter_++);
    TemplateAlignedSequenceContainer::addSequence(key, sequence);
  }

  using TemplateAlignedSequenceContainer::getAlphabet;

  using TemplateAlignedSequenceContainer::alphabet;

  using TemplateAlignedSequenceContainer::getSequenceNames;

  using TemplateAlignedSequenceContainer::getNumberOfSequences;

  using TemplateAlignedSequenceContainer::getNumberOfSites;

  using TemplateAlignedSequenceContainer::site;

  using TemplateAlignedSequenceContainer::deleteSite;

  using TemplateAlignedSequenceContainer::deleteSites;


  using TemplateAlignedSequenceContainer::hasSequence;

  using TemplateAlignedSequenceContainer::sequence;

  using TemplateAlignedSequenceContainer::removeSequence;

  using TemplateAlignedSequenceContainer::clear;

  bool hasSequenceForSpecies(const std::string& species) const
  {
    for (size_t i = 0; i < getNumberOfSequences(); ++i)
    {
      const MafSequence& seq = sequence(i);
      if (seq.getSpecies() == species)
        return true;
    }
    return false;
  }

  // Return the first sequence with the species name.
  const MafSequence& sequenceForSpecies(const std::string& species) const
  {
    for (size_t i = 0; i < getNumberOfSequences(); ++i)
    {
      const MafSequence& seq = sequence(i);
      if (seq.getSpecies() == species)
        return seq;
    }
    throw SequenceNotFoundException("MafBlock::sequenceForSpecies. No sequence with the given species name in this block.", species);
  }

  // Return all sequences with the species name.
  std::vector<const MafSequence*> getSequencesForSpecies(const std::string& species) const
  {
    std::vector<const MafSequence*> selection;
    for (size_t i = 0; i < getNumberOfSequences(); ++i)
    {
      const MafSequence* seq = &sequence(i);
      if (seq->getSpecies() == species)
        selection.push_back(seq);
    }
    return selection;
  }

  // Return the first sequence with the species name.
  std::unique_ptr<MafSequence> removeSequenceForSpecies(const std::string& species)
  {
    for (size_t i = 0; i < getNumberOfSequences(); ++i)
    {
      const MafSequence& seq = sequence(i);
      if (seq.getSpecies() == species)
        return removeSequence(i);
    }
    throw SequenceNotFoundException("MafBlock::removeSequenceForSpecies. No sequence with the given species name in this block.", species);
  }

  /**
   * @return The species names for all sequences in the container.
   */
  std::vector<std::string> getSpeciesList() const
  {
    std::vector<std::string> lst;
    for (size_t i = 0; i < getNumberOfSequences(); ++i)
    {
      lst.push_back(sequence(i).getSpecies());
    }
    return lst;
  }

  void removeCoordinatesFromSequence(size_t i)
  {
    // This is a bit of a trick, but avoid useless recopies.
    // It is safe here because the AlignedSequenceContainer is fully encapsulated.
    // It would not work if a VectorSiteContainer was used.
    const_cast<MafSequence&>(sequence(i)).removeCoordinates();
  }

  /**
   * @brief Add a new annotation to a sequence.
   *
   * @param i Sequence position in the container.
   * @param anno The annotation object to be added. 
   * @throw Exception If the annotation is not valid for this sequence.
   */
  void addAnnotationToSequence(size_t i, std::shared_ptr<SequenceAnnotation> anno)
  {
    sequence_(i).addAnnotation(anno);
  }
  
  /**
   * @brief Add a new annotation to a sequence for a given species.
   *
   * @param species Species name. In case several sequences are available, the first one will be used.
   * @param anno The annotation object to be added. 
   * @throw Exception If the annotation is not valid for this sequence.
   */
  void addAnnotationToSequenceForSpecies(const std::string& species, std::shared_ptr<SequenceAnnotation> anno)
  {
    sequenceForSpecies_(species).addAnnotation(anno);
  }

  std::string getDescription() const
  {
    std::string desc;
    desc += TextTools::toString(getNumberOfSequences()) + "x" + TextTools::toString(getNumberOfSites());
    return desc;
  }

  /**
   * @return True or False, if data are associated to the given property.
   * @param property The name of the property to look for.
   */
  bool hasProperty(const std::string& property) const
  {
    auto it = properties_.find(property);
    return it != properties_.end();
  }

  /**
   * @brief Get the data associated to a query property.
   *
   * @param property The property to look for.
   * @return The data associated to the given property.
   * @throw Exception if no data is associated to the given property.
   */
  const Clonable& getProperty(const std::string& property) const
  {
    auto it = properties_.find(property);
    if (it == properties_.end())
      throw Exception("MafBlock::getProperty. No data for property: " + property + " in block.");
    return *it->second;
  }

  /**
   * @brief Delete the data associated to a query property.
   *
   * @param property The property to look for.
   * @throw Exception if no data is associated to the given property.
   */
  void deleteProperty(const std::string& property)
  {
    auto it = properties_.find(property);
    if (it == properties_.end())
      throw Exception("MafBlock::deleteProperty. No data for property: " + property + " in block.");
    properties_.erase(it);
  }

  /**
   * @brief Set the data associated to a query property.
   *
   * An existing data associated to this property will be deleted and replaced by the new one.
   * @param property The property to look for.
   * @param data The data to associate to this property.
   * @throw Exception if the pointer toward the input data is NULL.
   */
  void setProperty(const std::string& property, std::unique_ptr<Clonable> data)
  {
    if (!data)
      throw Exception("MafBlock::setProperty. Pointer to data is nullptr.");
    if (hasProperty(property))
      deleteProperty(property);
    properties_[property] = std::move(data);
  }

private:
  using TemplateAlignedSequenceContainer::addSequence;

  // Return the first sequence with the species name.
  MafSequence& sequenceForSpecies_(const std::string& species)
  {
    for (size_t i = 0; i < getNumberOfSequences(); ++i)
    {
      MafSequence& seq = sequence_(i);
      if (seq.getSpecies() == species)
        return seq;
    }
    throw SequenceNotFoundException("MafBlock::sequenceForSpecies. No sequence with the given species name in this block.", species);
  }

  void deleteProperties_()
  {
    properties_.clear();
  }
};
} // end of namespace bpp.

#endif // _MAFBLOCK_H_
