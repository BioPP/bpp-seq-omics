// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _ABSTRACTMAFITERATOR_H_
#define _ABSTRACTMAFITERATOR_H_

#include "MafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>
#include <memory>

namespace bpp
{

/**
 * @brief Partial implementation of the MafIterator interface.
 *
 * This implements the listener parts.
 */
class AbstractMafIterator :
  public virtual MafIteratorInterface
{
protected:
  std::vector<std::unique_ptr<IterationListenerInterface>> iterationListeners_;
  bool started_;
  bool verbose_;

public:
  AbstractMafIterator() :
    iterationListeners_(),
    started_(false),
    verbose_(true)
  {}
  
  virtual ~AbstractMafIterator() {}

protected:
  AbstractMafIterator(const AbstractMafIterator& it):
    iterationListeners_(),
    started_(false),
    verbose_(it.verbose_)
  {}

  AbstractMafIterator& operator=(const AbstractMafIterator& it)
  {
    iterationListeners_.clear();
    started_ = false;
    verbose_ = it.verbose_;
    return *this;
  }

public:
  void addIterationListener(std::unique_ptr<IterationListenerInterface> listener)
  {
    iterationListeners_.push_back(move(listener));
  }

  std::unique_ptr<MafBlock> nextBlock()
  {
    if (!started_)
    {
      fireIterationStartSignal_();
      started_ = true;
    }
    auto block = analyseCurrentBlock_();
    if (block)
      fireIterationMoveSignal_(*block);
    else
      fireIterationStopSignal_();
    return block;
  }

  bool isVerbose() const { return verbose_; }
  void setVerbose(bool yn) { verbose_ = yn; }

protected:
  virtual std::unique_ptr<MafBlock> analyseCurrentBlock_() = 0;
  virtual void fireIterationStartSignal_();
  virtual void fireIterationMoveSignal_(const MafBlock& currentBlock);
  virtual void fireIterationStopSignal_();
};



/**
 * @brief Helper class for developping filter for maf blocks.
 */
class AbstractFilterMafIterator :
  public AbstractMafIterator
{
protected:
  std::shared_ptr<MafIteratorInterface> iterator_;
  std::unique_ptr<MafBlock> currentBlock_;
  std::shared_ptr<OutputStream> logstream_;

public:
  AbstractFilterMafIterator(std::shared_ptr<MafIteratorInterface> iterator) :
    AbstractMafIterator(),
    iterator_(std::move(iterator)), currentBlock_(nullptr),
    logstream_(ApplicationTools::message) {}

private:
  AbstractFilterMafIterator(const AbstractFilterMafIterator& it) :
    AbstractMafIterator(it),
    iterator_(it.iterator_), currentBlock_(nullptr),
    logstream_(it.logstream_) {}

  AbstractFilterMafIterator& operator=(const AbstractFilterMafIterator& it)
  {
    AbstractMafIterator::operator=(it);
    currentBlock_ = nullptr;
    iterator_  = it.iterator_;
    logstream_ = it.logstream_;
    return *this;
  }

public:
  void setLogStream(std::shared_ptr<OutputStream> logstream) { logstream_ = logstream; }
};


class TrashIteratorAdapter :
  public AbstractMafIterator
{
private:
  std::shared_ptr<MafTrashIteratorInterface> iterator_;

public:
  TrashIteratorAdapter(std::shared_ptr<MafTrashIteratorInterface> iterator) :
    iterator_(iterator) {}

private:
  TrashIteratorAdapter(const TrashIteratorAdapter& iterator) :
    iterator_(iterator.iterator_) {}

  TrashIteratorAdapter& operator=(const TrashIteratorAdapter& iterator)
  {
    iterator_ = iterator.iterator_;
    return *this;
  }

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_()
  {
    return iterator_->nextRemovedBlock();
  }
};


/**
 * @brief This special iterator synchronizes two adaptors.
 *
 * It takes as input a main iterator and a secondary one. The nextBlock method of the secondary iterator will be
 * called immediately after the one of the primary one. The resulting block of the main iterator will be forwarded,
 * while the one of the secondary iterator will be destroyed.
 */
class MafIteratorSynchronizer :
  public AbstractFilterMafIterator
{
private:
  std::shared_ptr<MafIteratorInterface> secondaryIterator_;

public:
  MafIteratorSynchronizer(
      std::shared_ptr<MafIteratorInterface> primaryIterator,
      std::shared_ptr<MafIteratorInterface> secondaryIterator) :
    AbstractFilterMafIterator(primaryIterator),
    secondaryIterator_(secondaryIterator)
  {}

private:
  MafIteratorSynchronizer(const MafIteratorSynchronizer& iterator) :
    AbstractFilterMafIterator(0),
    secondaryIterator_(iterator.secondaryIterator_)
  {}

  MafIteratorSynchronizer& operator=(const MafIteratorSynchronizer& iterator)
  {
    secondaryIterator_ = iterator.secondaryIterator_;
    return *this;
  }

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_()
  {
    currentBlock_ = iterator_->nextBlock();
    secondaryIterator_->nextBlock();
    return std::move(currentBlock_);
  }
};

} // end of namespace bpp.

#endif//_ABSTRACTMAFITERATOR_H_
