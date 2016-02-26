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
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <errno.h>
#include <assert.h>
#include <limits.h>

#include "libsakura/config.h"
#include "libsakura/localdef.h"
#include "libsakura/concurrent.h"

//#define LOG(x) do {}while(false)
#define LOG(x) fprintf(stderr, "Error: %d\n", (x))

namespace concurrent {
/* ======================= Mutex ======================= */
Mutex::Mutex() throw (int) {
	int result = pthread_mutex_init(&mutex_, NULL);
	if (result != 0) {
		LOG(result);
		throw result;
	}
}

Mutex::~Mutex() {
	int result = pthread_mutex_destroy(&mutex_);
	if (result != 0) {
		LOG(result);
	}
}

void Mutex::Lock() throw (int) {
	int result = pthread_mutex_lock(&mutex_);
	if (result != 0) {
		LOG(result);
		throw result;
	}
}

bool Mutex::TryLock() throw (int) {
	int result = pthread_mutex_trylock(&mutex_);
	if (result == 0) {
		return true;
	}
	if (result == EBUSY) {
		return false;
	}
	LOG(result);
	throw result;
}

void Mutex::Unlock() throw (int) {
	int result = pthread_mutex_unlock(&mutex_);
	if (result != 0) {
		LOG(result);
		throw result;
	}
}

/* ======================= Semaphore ======================= */
Semaphore::Semaphore(unsigned initial) throw (int) {
	//mutex_ = PTHREAD_MUTEX_INITIALIZER;
	//condition_ = PTHREAD_COND_INITIALIZER;
	semaphore_ = initial;

	int result = pthread_mutex_init(&mutex_, NULL);
	if (result != 0) {
		LOG(result);
		throw result;
	}

	result = pthread_cond_init(&condition_, NULL);
	if (result != 0) {
		pthread_mutex_destroy(&mutex_);
		LOG(result);
		throw result;
	}
}

Semaphore::~Semaphore() {
	int result = pthread_mutex_destroy(&mutex_);
	MARK_AS_USED(result);
	result = pthread_cond_destroy(&condition_);
}

void Semaphore::Up(unsigned amount) throw (int) {
	assert(0 < amount && amount <= UINT_MAX - semaphore_);
	int result = pthread_mutex_lock(&mutex_);
	if (result == 0) {
		semaphore_ += amount;
		result = pthread_cond_signal(&condition_);
		int result2 = pthread_mutex_unlock(&mutex_);
		if (result == 0 && result2 == 0) {
			return;
		}
		if (result != 0 && result2 != 0) {
			LOG(result);
			result = 0;
		}
		if (result == 0) {
			result = result2;
		}
	}
	LOG(result);
	throw result;
}

void Semaphore::Down(unsigned amount) throw (int) {
	assert(0 < amount);
	int result = pthread_mutex_lock(&mutex_);
	if (result == 0) {
		while (semaphore_ < amount) {
			result = pthread_cond_wait(&condition_, &mutex_);
			if (result != 0) {
				LOG(result);
				break;
			}
		}
		if (semaphore_ >= amount) {
			semaphore_ -= amount;
		}
		int result2 = pthread_mutex_unlock(&mutex_);
		if (result == 0 && result2 == 0) {
			return;
		}
		if (result != 0 && result2 != 0) {
			LOG(result);
			result = 0;
		}
		if (result == 0) {
			result = result2;
		}
	}
	LOG(result);
	throw result;
}

/* ======================= Broker ======================= */
Broker::Broker(bool (*producer)(void *context) throw (PCException),
void (*consumer)(void *context) throw (PCException)) {
	this->producer_ = producer;
	this->consumer_ = consumer;
}

Broker::~Broker() {
}

void Broker::EnableNested() {
#ifdef _OPENMP
	omp_set_nested(1);
#endif
}

void Broker::DisableNested() {
#ifdef _OPENMP
	omp_set_nested(0);
#endif
}

void Broker::SetNestedState(bool nested) {
#ifdef _OPENMP
	omp_set_nested(static_cast<int>(nested));
#endif
}

bool Broker::GetNestedState() {
#ifdef _OPENMP
	return omp_get_nested() ? true : false;
#else
	return false;
#endif
}

void Broker::RunProducerAsMasterThread(void *context, unsigned do_ahead)
		throw (PCException) {
	_Run(context, do_ahead, kProdAsMaster);
}

void Broker::RunConsumerAsMasterThread(void *context, unsigned do_ahead)
		throw (PCException) {
	_Run(context, do_ahead, kConsAsMaster);
}

void Broker::Run(void *context, unsigned do_ahead) throw (PCException) {
	_Run(context, do_ahead, kUnspecified);
}

void Broker::_Run(void *context, unsigned do_ahead, ThreadSpec thread_spec)
		throw (PCException) {
	assert(do_ahead > 0);
#if defined(_OPENMP)
	PCException const *prod_ex = NULL;
	PCException const *cons_ex = NULL;
	unsigned queued_jobs = 0;
	int consumer_terminated = 0;
	Semaphore semaphore_for_consumer;
	Semaphore semaphore_for_producer(do_ahead);

#pragma omp parallel num_threads(2) \
		shared(semaphore_for_consumer, semaphore_for_producer, \
			   consumer_terminated, queued_jobs)
	{
		//fprintf(stderr, "run: %p %d\n", context, omp_get_thread_num());
		bool run_prod = true;
		if (thread_spec == kUnspecified) {
#pragma omp single
			{
				run_prod = false;
			}
		} else {
			bool is_master = false;
#pragma omp master
			{
				is_master = true;
			}
			if (thread_spec == kProdAsMaster) {
				run_prod = is_master;
			} else { // ConsAsMaster
				run_prod = ! is_master;
			}
		}

		if (run_prod) { // producer
			for (;;) {
				semaphore_for_producer.Down();
				int consumers_dead = 0;
#pragma omp atomic
				consumers_dead += consumer_terminated;
				if (consumers_dead) {
					break;
				}
				try {
					bool produced = producer_(context);
					if (! produced) {
						break;
					}
				} catch (PCException const &e) {
					prod_ex = &e;
					break;
				}
#pragma omp atomic
				++queued_jobs;
				semaphore_for_consumer.Up();
			}
			// additional 'up' to give consumer a chance to terminate.
			semaphore_for_consumer.Up();
		} else { // consumer
			for (;;) {
				semaphore_for_consumer.Down();
				unsigned remaining_jobs = 0U;
#pragma omp atomic
				remaining_jobs += queued_jobs;
				if (remaining_jobs == 0U) {
					break;
				}
#pragma omp atomic
				--queued_jobs;
				try {
					consumer_(context);
				} catch (PCException const &e) {
					cons_ex = &e;
					break;
				}
				semaphore_for_producer.Up();
			}
#pragma omp atomic
			++consumer_terminated;
			// additional 'up' to give producer a chance to terminate.
			semaphore_for_producer.Up();
		}
	}
	if (prod_ex) {
		prod_ex->Raise();
	} else if (cons_ex) {
		cons_ex->Raise();
	}
#else
	RunSequential(context);
#endif
}

void Broker::RunSequential(void *context) throw (PCException) {
	for (;;) {
		bool produced = producer_(context);
		if (!produced) {
			break;
		}
		consumer_(context);
	}
}

} // namespace
