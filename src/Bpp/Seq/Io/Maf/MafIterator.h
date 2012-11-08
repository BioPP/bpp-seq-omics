//
// File: MafIterator.h
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

#ifndef _MAFITERATOR_H_
#define _MAFITERATOR_H_

#include "MafBlock.h"

//From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp {

//Forward declaration:
class IterationListener;

/**
 * @brief Interface to loop over maf alignment blocks.
 */
class MafIterator
{
  public:
    virtual ~MafIterator() {}

  public:
    /**
     * @brief Get the next available alignment block.
     *
     * @return A maf alignment block, or a null pointer if no more block is available.
     */
    virtual MafBlock* nextBlock() throw (Exception) = 0;

    virtual bool verbose() const = 0;
    
    virtual void verbose(bool yn) = 0;
    
    virtual void addIterationListener(IterationListener* listener) = 0;
    

};

/**
 * @brief Partial implementation of the MafIterator interface.
 *
 * This implements the listener parts.
 */
class AbstractMafIterator:
  public virtual MafIterator
{
  protected:
    std::vector<IterationListener*> iterationListeners_;
    bool started_;
    bool verbose_;

  public:
    AbstractMafIterator(): iterationListeners_(), started_(false), verbose_(true) {}
    
    virtual ~AbstractMafIterator() {}

  public:
    void addIterationListener(IterationListener* listener) {
      iterationListeners_.push_back(listener);
    }

    MafBlock* nextBlock() throw (Exception) {
      if (!started_) {
        fireIterationStartSignal_();
        started_ = true;
      }
      MafBlock* block = analyseCurrentBlock_();
      if (block)
        fireIterationMoveSignal_(*block);
      else
        fireIterationStopSignal_();
      return block;
    }
    
    bool verbose() const { return verbose_; }
    void verbose(bool yn) { verbose_ = yn; }

  protected:
    virtual MafBlock* analyseCurrentBlock_() = 0;
    virtual void fireIterationStartSignal_();
    virtual void fireIterationMoveSignal_(const MafBlock& currentBlock);
    virtual void fireIterationStopSignal_();

};

/**
 * @brief Interface to loop over removed blocks of a maf alignment.
 */
class MafTrashIterator
{
  public:
    virtual ~MafTrashIterator() {}

  public:
    /**
     * @brief Get the next available removed alignment block.
     *
     * @return A maf alignment block, or a null pointer if no more block is available.
     */
    virtual MafBlock* nextRemovedBlock() throw (Exception) = 0;
    
};


/**
 * @brief Helper class for developping filter for maf blocks.
 */
class AbstractFilterMafIterator:
  public AbstractMafIterator
{
  protected:
    MafIterator* iterator_;
    MafBlock* currentBlock_;
    OutputStream* logstream_;

  public:
    AbstractFilterMafIterator(MafIterator* iterator) :
      AbstractMafIterator(),
      iterator_(iterator), currentBlock_(0),
      logstream_(ApplicationTools::message) {}

  private:
    AbstractFilterMafIterator(const AbstractFilterMafIterator& it):
      AbstractMafIterator(it), 
      iterator_(it.iterator_), currentBlock_(0),
      logstream_(it.logstream_) {}

    AbstractFilterMafIterator& operator=(const AbstractFilterMafIterator& it) {
      AbstractMafIterator::operator=(it);
      currentBlock_ = 0;
      iterator_  = it.iterator_;
      logstream_ = it.logstream_;
      return *this;
    }

  public:
    void setLogStream(OutputStream* logstream) { logstream_ = logstream; }

};


class TrashIteratorAdapter:
  public AbstractMafIterator
{
  private:
    MafTrashIterator* iterator_;

  public:
    TrashIteratorAdapter(MafTrashIterator* iterator) :
      iterator_(iterator) {}

  private:
    TrashIteratorAdapter(const TrashIteratorAdapter& iterator) :
      iterator_(iterator.iterator_) {}
    
    TrashIteratorAdapter& operator=(const TrashIteratorAdapter& iterator) {
      iterator_ = iterator.iterator_;
      return *this;
    }

  private:
    MafBlock* analyseCurrentBlock_() throw (Exception) {
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
class MafIteratorSynchronizer:
  public AbstractFilterMafIterator
{
  private:
    MafIterator* secondaryIterator_;

  public:
    MafIteratorSynchronizer(MafIterator* primaryIterator, MafIterator* secondaryIterator) :
      AbstractFilterMafIterator(primaryIterator), secondaryIterator_(secondaryIterator)
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
    MafBlock* analyseCurrentBlock_() throw (Exception) {
      currentBlock_ = iterator_->nextBlock();
      MafBlock* secondBlock = secondaryIterator_->nextBlock();
      if (secondBlock)
        delete secondBlock;
      return currentBlock_;
    }

};


} // end of namespace bpp.

#endif //_MAFITERATOR_H_
