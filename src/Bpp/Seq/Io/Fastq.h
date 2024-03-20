// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _BPP_SEQ_IO_FASTQ_H_
#define _BPP_SEQ_IO_FASTQ_H_

#include <string>
#include <Bpp/Seq/Io/ISequenceStream.h>
#include <Bpp/Seq/Io/OSequenceStream.h>
#include <Bpp/Seq/Sequence.h>
#include <Bpp/Seq/SequenceWithQuality.h>

namespace bpp
{
/**
 * @brief The fastq sequence file format.
 *
 * @author Sylvain Gaillard
 */
class Fastq :
  public virtual ISequenceWithQualityStream,
  public virtual OSequenceWithQualityStream
{
private:
  bool repeatName_;

public:
  /**
   * @brief Build a new Fastq object.
   *
   * @param repName Tell if the names in the file is repeated (tested
   * on input) or must be repeated (for output).
   */
  Fastq(bool repName = false) : repeatName_(repName) {}

  // Class destructor
  virtual ~Fastq() {}

public:
  /**
   * @name The IOSequence interface.
   *
   * @{
   */
  const std::string getFormatName() const override { return "FASTQ file"; }
  const std::string getFormatDescription() const override
  {
    return "Sequence with quality";
  }
  const std::string getDataType() const override { return "Sequence with quality"; }
  /** @} */
  bool repeatName() const { return repeatName_; }
  void repeatName(bool yn) { repeatName_ = yn; }

  /**
   * @name The ISequenceStream interface.
   *
   * @{
   */
  /**
   * @copydoc bpp::ISequenceStream::nextSequence()
   * @author Sylvain Gaillard
   *
   * @par Usage
   *
   * @code
   * // Creating a SequenceWithQuality object
   * DNA alpha;
   * SequenceWithQuality seq(&alpha);
   *
   * // Create a FastQ parser
   * Fastq fq;
   *
   * // Opening the file
   * std::ifstream in("reads.fastq");
   *
   * // Read the sequences
   * while (fq.nextSequence(in, seq)) {
   *   // ... do something with the sequence ...
   * }
   *
   * // Close the file
   * in.close();
   * @endcode
   */
  bool nextSequence(std::istream& input, SequenceWithQuality& seq) const override;
  /** @} */

  /**
   * @name The OSequenceStream interface.
   *
   * @{
   */
  /**
   * @copydoc OSequenceStream::writeSequence()
   * @author Sylvain Gaillard
   */
  void writeSequence(std::ostream& output, const SequenceWithQuality& seq) const override;
  /** @} */
};
}

#endif // _BPP_SEQ_IO_FASTQ_H_
