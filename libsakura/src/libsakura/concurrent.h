/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2016
 * National Astronomical Observatory of Japan
 * 2-21-1, Osawa, Mitaka, Tokyo, 181-8588, Japan.
 * 
 * This file is part of Sakura.
 * 
 * Sakura is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the 
 * Free Software Foundation, either version 3 of the License, or (at your 
 * option) any later version.
 * 
 * Sakura is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
 * License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License 
 * along with Sakura.  If not, see <http://www.gnu.org/licenses/>.
 * @SAKURA_LICENSE_HEADER_END@
 */
//
// C++ Interface: concurrent
//
// Description:
//
//
// Author: Kohji Nakamura <k.nakamura@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef LIBSAKURA_LIBSAKURA_CONCURRENT_H_
#define LIBSAKURA_LIBSAKURA_CONCURRENT_H_

#include <assert.h>
#include <pthread.h>

namespace concurrent {
/**
 * Producer/Consumer Exception
 */
class PCException {
public:
	virtual ~PCException() {
	}
	virtual void Raise() const /*throw(PCException)*/{
		throw *this;
	}
};

//typedef void (*produce_t)(void *context) throw(PCException);
//typedef void (*consume_t)(void *context) throw(PCException);

class Mutex {
public:
	Mutex() throw (int);
	virtual ~Mutex();
	void Lock() throw (int);
	/**
	 * Returns true if this thread could lock.
	 * Returns false if already locked.
	 */
	bool TryLock() throw (int);
	void Unlock() throw (int);
private:
	Mutex(Mutex const &other);
	Mutex &operator =(Mutex const &other);

	pthread_mutex_t mutex_;
};

class Semaphore {
public:
	explicit Semaphore(unsigned initial = 0U) throw (int);
	virtual ~Semaphore();
	void Up(unsigned amount = 1U) throw (int);
	void Down(unsigned amount = 1U) throw (int);
private:
	Semaphore(Semaphore const &other);
	Semaphore &operator =(Semaphore const &other);

	pthread_mutex_t mutex_;
	pthread_cond_t condition_;
	unsigned semaphore_;
};

class FIFOException {
public:
	virtual ~FIFOException() {
	}
	virtual void Raise() const {
		throw *this;
	}
};

class EmptyException {
public:
	virtual ~EmptyException() {
	}
	virtual void Raise() const {
		throw *this;
	}
};

class FullException {
public:
	virtual ~FullException() {
	}
	virtual void Raise() const {
		throw *this;
	}
};

template<class T, size_t N>
class FIFO {
public:
	FIFO() {
		assert(N > 0);
		Reset();
	}
	virtual ~FIFO() {
		mutex_.Lock(); // wait until current lock is released.
		mutex_.Unlock();
	}

	virtual void Lock() {
		mutex_.Lock();
	}

	/**
	 * @return true if this thread could lock. false if already locked.
	 */
	virtual bool TryLock() {
		return mutex_.TryLock();
	}

	virtual void Unlock() {
		mutex_.Unlock();
	}

	virtual void Clear() {
		Reset();
	}

	virtual void Put(T const &value) throw (FullException) {
		size_t new_tail = Wrap(tail_ + 1);
		if (head_ == new_tail) {
			throw FullException();
		}
		elements_[tail_] = value;
		tail_ = new_tail;
	}

	virtual T Get() throw (EmptyException) {
		if (head_ == tail_) {
			throw EmptyException();
		}
		T result = elements_[head_];
		head_ = Wrap(head_ + 1);
		return result;
	}

	/**
	 * Returns capacity size.
	 */
	virtual size_t Size() const {
		return N;
	}

	/**
	 * Returns number of elements in this FIFO.
	 */
	virtual size_t Length() const {
		size_t result = tail_ - head_;
		if (head_ > tail_) {
			result = N + 1 - (head_ - tail_);
		}
		return result;
	}
private:
	FIFO(FIFO const &other);
	FIFO &operator =(FIFO const &other);

	size_t Wrap(size_t n) {
		return n % (N + 1);
	}

	void Reset() {
		head_ = tail_ = 0;
	}

	T elements_[N + 1]; // +1 to make an implementation simple.
	Mutex mutex_;
	size_t head_;
	size_t tail_;

};

class Broker {
public:
	Broker(bool (*producer)(void *context) throw (PCException),
	void (*consumer)(void *context) throw (PCException));
	virtual ~Broker();
	static void EnableNested();
	static void DisableNested();
	static void SetNestedState(bool nested);
	static bool GetNestedState();

	virtual void Run(void *context, unsigned do_ahead = 1) throw (PCException);
	virtual void RunProducerAsMasterThread(void *context, unsigned do_ahead = 1)
			throw (PCException);
	virtual void RunConsumerAsMasterThread(void *context, unsigned do_ahead = 1)
			throw (PCException);
	virtual void RunSequential(void *context) throw (PCException);
protected:
	enum ThreadSpec {
		kProdAsMaster, kConsAsMaster, kUnspecified
	};
	bool (*producer_)(void *context) throw (PCException);
	void (*consumer_)(void *context) throw (PCException);
	virtual void _Run(void *context, unsigned do_ahead, ThreadSpec thread_spec)
			throw (PCException);
};

#if 1
template<class Context, class Product>
class Producer {
public:
	virtual Product Produce(Context *ctx) throw (PCException) = 0;
	virtual ~Producer() {
	}
};

template<class Context, class Product>
class Consumer {
public:
	virtual void Consume(Context *ctx, Product const*product)
			throw (PCException) = 0;
	virtual ~Consumer() {
	}
};

class ProdCons {
public:
	virtual ~ProdCons() {
	}

	virtual void RunProducerAsMasterThread(void *context) throw (PCException) = 0;
	virtual void RunConsumerAsMasterThread(void *context) throw (PCException) = 0;
	virtual void Produce(void *context) throw (PCException) = 0;
	virtual void Consume(void *context) throw (PCException) = 0;

	/**
	 * @ref Produce() should  call this method to
	 * let consumer to know it is Ready.
	 */
	virtual void Ready() = 0;

	/**
	 * Returns true if ready,
	 * otherwise, finished or there was an error, returns false.
	 *
	 * @ref Consume() should call this method
	 * to check if it is ready or not.
	 * If false is returned, @ref IsError() or/and
	 * @ref IsFinished() should be called to see the reason.
	 */
	virtual bool WaitForReady() = 0;

	/**
	 * Returns true if @ref Produce() was finished,
	 * otherwise returns false.
	 */
	virtual bool IsFinished() = 0;

	/**
	 * @ref Produce() / @ref Consume() should call this method
	 * to report error to @ref Consume() / @ref Produce() .
	 */
	virtual void ReportError(void *error_info) = 0;

	/**
	 * Returns true and set error_info if @ref ReportError() is called,
	 * otherwise returns false and set error_info to NULL.
	 */
	virtual bool IsError(void **error_info) = 0;
};
#endif

} // namespace concurrent

#endif	/* LIBSAKURA_LIBSAKURA_CONCURRENT_H_ */
