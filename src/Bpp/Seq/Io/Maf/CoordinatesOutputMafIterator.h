// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _COORDINATESOUTPUTMAFITERATOR_H_
#define _COORDINATESOUTPUTMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>

namespace bpp
{
/**
 * @brief Output each sequence coordinates for each block.
 *
 * The set of species for which coordinates are output is provided as argument.
 * The current implementation outputs results as a table to a file. Later implementation
 * may involve a dedicated data structure for other application usage (file indexing and so one).
 */
class CoordinatesOutputMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::shared_ptr<std::ostream> output_;
  std::vector<std::string> species_;
  bool includeSrcSize_;

public:
  /**
   * @brief Creates a new CoordinatesOutputMafIterator object.
   *
   * @param iterator The input iterator.
   * @param out A pointer toward the output stream. The stream will not be own by this instance, and will not be copied neither destroyed.
   * @param species A vector of species names for which coordinates should be output. In case of missing species for one block, NA will be produced.
   * @param includeSrcSize Tell if source size should also be written (useful to convert coordinates on the negative strand).
   */
  CoordinatesOutputMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      std::shared_ptr<std::ostream> out,
      const std::vector<std::string>& species,
      bool includeSrcSize = false) :
    AbstractFilterMafIterator(iterator),
      output_(out),
	species_(species),
	includeSrcSize_(includeSrcSize)
  {
    if (output_)
      writeHeader_(*output_);
  }

private:
  CoordinatesOutputMafIterator(const CoordinatesOutputMafIterator& iterator) :
    AbstractFilterMafIterator(0),
    output_(iterator.output_),
    species_(iterator.species_),
    includeSrcSize_(iterator.includeSrcSize_)
  {}

  CoordinatesOutputMafIterator& operator=(const CoordinatesOutputMafIterator& iterator)
  {
    output_ = iterator.output_;
    species_ = iterator.species_;
    includeSrcSize_ = iterator.includeSrcSize_;
    return *this;
  }

private:
  void writeHeader_(std::ostream& out) const;
  std::unique_ptr<MafBlock> analyseCurrentBlock_();
};
} // end of namespace bpp.

#endif//_COORDINATESOUTPUTMAFITERATOR_H_
