//
// File: OutputAlignmentMafIterator.h
// Authors: Julien Dutheil
// Created: Tue Sep 07 2010
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

#ifndef _OUTPUTALIGNMENTMAFITERATOR_H_
#define _OUTPUTALIGNMENTMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From bpp-seq:
#include <Bpp/Seq/Io/OSequence.h>

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief This iterator forward the iterator given as input after having printed its content to an alignment file.
 * The syntax for ENSEMBL meta data is used.
 */
class OutputAlignmentMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::shared_ptr<std::ostream> output_;
  std::string file_;
  bool mask_;
  bool outputCoordinates_;
  bool addLDHatHeader_;
  std::unique_ptr<OAlignment> writer_;
  unsigned int currentBlockIndex_;
  std::string refSpecies_;

public:
  /**
   * @brief Creates a new OutputAlignmentMafIterator object.
   *
   * All block will be printed as separate alignment, yet one after the other on the stream.
   * Be aware that not all format will recognize the resulting file as a multiple alignment file (Mase and Clustal will for instance, not Fasta).
   * @param iterator The input iterator
   * @param out A pointer toward the output stream. The stream will not be own by this instance, and will not be copied neither destroyed.
   * @param writer A pointer toward an alignment writer object which specifies the format to use when writing sequences.
   * The underlying object will be own by this instance, and destroyed when this object is deleted.
   * @param mask Tell if sequences should be printed masked (if applicable).
   * @param outputCoordinates Tell if coordinates should be written in sequence headers, if any.
   * @param addLDHatHeader Tell if first line of file should contain number and lenght of sequences (for instance for use with LDhat/convert).
   * @param reference [optional] specify a reference species which can be used to configure file names
   * (for instance using coordinates information).
   */
  OutputAlignmentMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      std::shared_ptr<std::ostream> out,
      std::unique_ptr<OAlignment> writer,
      bool mask = true,
      bool outputCoordinates = true,
      bool addLDHatHeader = false,
      const std::string& reference = "") :
    AbstractFilterMafIterator(iterator),
    output_(out),
    file_(),
    mask_(mask),
    outputCoordinates_(outputCoordinates),
    addLDHatHeader_(addLDHatHeader),
    writer_(move(writer)),
    currentBlockIndex_(0),
    refSpecies_(reference)
  {
    if (!writer)
      throw Exception("OutputAlignmentMafIterator (constructor 1): sequence writer should not be a NULL pointer!");
  }

  /**
   * @brief Creates a new OutputAlignmentMafIterator object.
   *
   * All block will be printed as separate alignment, yet one after the other on the stream.
   * Be aware that not all format will recognize the resulting file as a multiple alignment file (Mase and Clustal will for instance, not Fasta).
   * @param iterator The input iterator
   * @param file A string describing the path to the output files. Each block will be written to a distinct file.
   * If "file" is a fixed string, it will only contain the last block. Using the %i code in the file name allows to generate one file per block, %i denoting the block index.
   * @param writer A pointer toward an alignment writer object which specifies the format to use when writing sequences.
   * The underlying object will be own by this instance, and destroyed when this object is deleted.
   * @param mask Tell if sequences should be printed masked (if applicable).
   * @param outputCoordinates Tell if coordinates should be written in sequence headers, if any.
   * @param addLDHatHeader Tell if first line of file should contain number and lenght of sequences (for instance for use with LDhat/convert).
   * @param reference [optional] specify a reference species which can be used to configure file names
   * (for instance using coordinates information).
   */
  OutputAlignmentMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      const std::string& file,
      std::unique_ptr<OAlignment> writer,
      bool mask = true,
      bool outputCoordinates = true,
      bool addLDHatHeader = false,
      const std::string& reference = "") :
    AbstractFilterMafIterator(iterator),
    output_(0),
    file_(file),
    mask_(mask),
    outputCoordinates_(outputCoordinates),
    addLDHatHeader_(addLDHatHeader),
    writer_(move(writer)),
    currentBlockIndex_(0),
    refSpecies_(reference)
  {
    if (!writer)
      throw Exception("OutputAlignmentMafIterator (constructor 2): sequence writer should not be a NULL pointer!");
  }

  virtual ~OutputAlignmentMafIterator() {}

private:
  OutputAlignmentMafIterator(const OutputAlignmentMafIterator& iterator) :
    AbstractFilterMafIterator(0),
    output_(iterator.output_),
    file_(iterator.file_),
    mask_(iterator.mask_),
    outputCoordinates_(iterator.outputCoordinates_),
    addLDHatHeader_(iterator.addLDHatHeader_),
    writer_(),
    currentBlockIndex_(iterator.currentBlockIndex_),
    refSpecies_(iterator.refSpecies_)
  {}

  OutputAlignmentMafIterator& operator=(const OutputAlignmentMafIterator& iterator)
  {
    output_ = iterator.output_;
    file_   = iterator.file_;
    mask_   = iterator.mask_;
    outputCoordinates_ = iterator.outputCoordinates_;
    addLDHatHeader_ = iterator.addLDHatHeader_;
    writer_.release();
    currentBlockIndex_ = iterator.currentBlockIndex_;
    refSpecies_ = iterator.refSpecies_;
    return *this;
  }

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_();

  void writeBlock(std::ostream& out, const MafBlock& block) const;
};
} // end of namespace bpp.

#endif//_OUTPUTALIGNMENTMAFITERATOR_H_
