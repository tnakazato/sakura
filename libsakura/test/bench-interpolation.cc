/*
 * bench-template.cc
 *
 *  g++ -std=c++11 -O2 -o bench -Iinstalled/include ../test/bench-template.cc -rdynamic -L installed/lib -Wl,-rpath,installed/lib -lsakura -lpthread
 */

#include <libsakura/sakura.h>
#include <unistd.h>
#include <iostream>
#include <algorithm>
#include <array>
#include <memory>
#include <deque>
#include <chrono>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <sys/types.h>
#include <pthread.h>
#include <sched.h>

#define ELEMENTSOF(x) (sizeof(x) / sizeof((x)[0]))

namespace {

double GetCurrentTime() noexcept {
	using namespace std::chrono;
	auto now = system_clock::now().time_since_epoch();
	return duration_cast<duration<double, seconds::period>>(now).count();
}

template<typename T>
T *alloc_obj() {
	void *memptr = 0;
	int result = posix_memalign(&memptr, 32, sizeof(T));
	if (result != 0) {
		throw std::bad_alloc();
	}
	auto p = static_cast<T*>(memptr);
	return new (p) T();
}

template<typename T>
void free_obj(T* ptr) {
	if (ptr) {
		ptr->~T();
		free(ptr);
	}
}

template<size_t SIZE>
struct InterpolateXBench {
	struct WSType {
		static constexpr size_t num_interpolated = SIZE;
		static constexpr size_t num_base = 100;
		static constexpr size_t num_array = 1;
		static constexpr size_t N_in = num_base * num_array;
		static constexpr size_t N_out = num_interpolated * num_array;
		std::unique_ptr<std::array<double, num_base>,
				decltype(&free_obj<std::array<double, num_base> >)> base_position;
		std::unique_ptr<std::array<float, N_in>,
				decltype(&free_obj<std::array<float, N_in> >)> in_data;
		std::unique_ptr<std::array<bool, N_in>,
				decltype(&free_obj<std::array<bool, N_in> >)> in_mask;
		std::unique_ptr<std::array<double, num_interpolated>,
				decltype(&free_obj<std::array<double, num_interpolated> >)> interpolated_position;
		std::unique_ptr<std::array<float, N_out>,
				decltype(&free_obj<std::array<float, N_out> >)> out_data;
		std::unique_ptr<std::array<bool, N_out>,
				decltype(&free_obj<std::array<bool, N_out> >)> out_mask;
		WSType() :
				base_position(alloc_obj<std::array<double, num_base> >(),
						free_obj<std::array<double, num_base> >), in_data(
						alloc_obj<std::array<float, N_in> >(),
						free_obj<std::array<float, N_in> >), in_mask(
						alloc_obj<std::array<bool, N_in> >(),
						free_obj<std::array<bool, N_in> >), interpolated_position(
						alloc_obj<std::array<double, num_interpolated> >(),
						free_obj<std::array<double, num_interpolated> >), out_data(
						alloc_obj<std::array<float, N_out> >(),
						free_obj<std::array<float, N_out> >), out_mask(
						alloc_obj<std::array<bool, N_out> >(),
						free_obj<std::array<bool, N_out> >) {

		}
	};

	static WSType *createWS(unsigned id) {
		WSType *ws = new WSType();
		constexpr float amplitude = 10.0f;
		const double delta_x = 101.0/static_cast<double>(ws->num_interpolated-1);
		for (size_t i = 0; i < ws->num_base; ++i) {
			ws->base_position->data()[i] = static_cast<double>(i+1);
		}
		for (size_t i = 0; i < ws->num_interpolated; ++i) {
			ws->interpolated_position->data()[i] = static_cast<double>(i)*delta_x;
		}
		for (size_t j=0; j < ws->N_in; ++j) {
			size_t i = j % ws->num_base;
			float sign = j % 2 == 0 ? +1.0f : -1.0f;
			ws->in_data->data()[j] = static_cast<float>(i+1)+sign*amplitude;
		}
		ws->in_mask->fill(true);
		if (false) {// debug output
			if (id==0) {
				std::cout << "base_position = [" ;
				for (size_t i = 0; i < ws->base_position->size(); ++i)
					std::cout << ws->base_position->data()[i] << ", ";
				std::cout << "]"<< std::endl;
				std::cout << "interpolated_position = [";
				for (size_t i = 0; i < ws->interpolated_position->size(); ++i)
					std::cout << ws->interpolated_position->data()[i] << ", ";
				std::cout << "]"<< std::endl;
				std::cout << "in_data = [";
				for (size_t i = 0; i < ws->base_position->size(); ++i)
					std::cout << ws->in_data->data()[i] << ", ";
				std::cout << "]"<< std::endl;
			}
		}
		return ws;
	}

	static void destroyWS(unsigned id, WSType *ws) {
		delete ws;
	}

	struct ThreadLocalType { // typical member is thread local context object
		unsigned thread_id;
	};

	static ThreadLocalType *createThreadLocal(unsigned tread_id) {
		ThreadLocalType *tl = new ThreadLocalType;
		tl->thread_id = tread_id;
		return tl;
	}

	static void destroyThreadLocal(unsigned tread_id,
			ThreadLocalType *threadLocal) {
		delete threadLocal;
	}

	static void run(ThreadLocalType *tl, WSType *ws, unsigned id) {
		LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method = LIBSAKURA_SYMBOL(InterpolationMethod_kLinear);
		LIBSAKURA_SYMBOL(Status) status =
		LIBSAKURA_SYMBOL(InterpolateXAxisFloat)(interpolation_method, 0, ws->base_position->size(),
				ws->base_position->data(), ws->in_data->size()/ws->base_position->size(),
				ws->in_data->data(), ws->in_mask->data(), ws->interpolated_position->size(),
				ws->interpolated_position->data(),
				ws->out_data->data(), ws->out_mask->data());
		assert(status == LIBSAKURA_SYMBOL(Status_kOK));
	}
};

class BenchBase {
protected:
	unsigned numWorkingSet;
	unsigned numThreads;
	BenchBase(unsigned numWorkingSet, unsigned numThreads,
			unsigned iterPerThread) :
			numWorkingSet(numWorkingSet), numThreads(numThreads) {
		assert(numWorkingSet >= numThreads);
	}
public:
	virtual ~BenchBase() {
	}
	virtual void run() = 0;
};

class Protect {
	pthread_mutex_t *mutex;
public:
	Protect(pthread_mutex_t *mutex) :
			mutex(mutex) {
		lock();
	}
	~Protect() {
		unlock();
	}
	void lock() {
		int result = pthread_mutex_lock(mutex);
		assert(result == 0);
	}
	void unlock() {
		int result = pthread_mutex_unlock(mutex);
		if (result != 0) {
			char buf[100];
			auto ptr = strerror_r(result, buf, sizeof buf);
			std::cerr << "Error: " << buf << std::endl;
			throw "error";
		}
	}
}
;

template<typename T>
class Bench: public BenchBase {
	struct Context {
		Bench<T> *self;
		unsigned thread_id;
	};
	std::unique_ptr<typename T::WSType *[]> ws;
	std::unique_ptr<typename T::ThreadLocalType *[]> tl;
	std::unique_ptr<pthread_t[]> threads;
	std::unique_ptr<Context[]> contexts;
	pthread_mutex_t queue_mutex;
	std::deque<unsigned> queue;
	pthread_mutex_t count_mutex;
	unsigned long long numIters;
	unsigned long long count;
	volatile bool stop;

	void incr() {
		Protect protect(&count_mutex);
		++count;
	}
	unsigned long long counter() {
		Protect protect(&count_mutex);
		return count;
	}
	unsigned get() {
		Protect protect(&queue_mutex);
		while (queue.empty()) {
			protect.unlock();
			std::cerr << "queue stalled" << std::endl;
			usleep(10);
			protect.lock();
		}
		auto entry = queue.front();
		queue.pop_front();
		return entry;
	}
	void put(unsigned v) {
		Protect protect(&queue_mutex);
		queue.push_back(v);
	}
	static void *job(Context *arg) {
		char buf[100];
		auto self = arg->self;
		if (false) {
			sprintf(buf, "thr %u: started\n", arg->thread_id);
			std::cout << buf;
		}
		try {
			while (!self->stop) {
				auto working_set = self->get();
				if (false) {
					sprintf(buf, "thr %u: %u\n", arg->thread_id, working_set);
					std::cout << buf;
				}
				T::run(self->tl[arg->thread_id], self->ws[working_set],
						working_set);
				self->put(working_set);
				self->incr();
				if (false) {
					sprintf(buf, "thr end %u: %u\n", arg->thread_id,
							working_set);
					std::cout << buf;
				}
			}
		} catch (...) {
			std::cerr << "Exception was thrown\n";
			exit(1);
		}
		return nullptr;
	}
	void init() {
		for (unsigned i = 0; i < numWorkingSet; ++i) {
			ws[i] = T::createWS(i);
			queue.push_back(i);
		}

		for (unsigned i = 0; i < numThreads; ++i) {
			tl[i] = T::createThreadLocal(i);
			contexts[i].self = this;
			contexts[i].thread_id = i;
		}
	}
public:
	Bench(unsigned numWorkingSet, unsigned numThreads,
			unsigned long long numIters) :
			BenchBase(numWorkingSet, numThreads, numIters), ws(
					new typename T::WSType *[numWorkingSet]), tl(
					new typename T::ThreadLocalType*[numThreads]), threads(
					new pthread_t[numThreads]), contexts(
					new Context[numThreads]), queue(), numIters(numIters), count(
					0), stop(
			false) {
		int result = pthread_mutex_init(&queue_mutex, nullptr);
		assert(result == 0);
		result = pthread_mutex_init(&count_mutex, nullptr);
		assert(result == 0);
		init();
	}
	~Bench() {
		for (unsigned i = 0; i < numWorkingSet; ++i) {
			T::destroyWS(i, ws[i]);
		}

		for (unsigned i = 0; i < numThreads; ++i) {
			T::destroyThreadLocal(i, tl[i]);
		}
	}
	void run() {
		auto start = GetCurrentTime();
		auto cnt = counter();
		for (unsigned i = 0; i < numThreads; ++i) {
			pthread_create(&threads[i], nullptr,
					reinterpret_cast<void*(*)(void *)>(job), &contexts[i]);
		}

		while (cnt < numIters) {
			sleep(1);
			auto now = GetCurrentTime();
			cnt = counter();
			std::cout << (cnt / (now - start)) << " workingSets/sec\n";
		}
		stop = true;
		if (false) {// debug output
			std::cout << "out_data = [";
			for (size_t i = 0; i < ws[0]->interpolated_position->size(); ++i)
				std::cout << ws[0]->out_data->data()[i] << ", ";
			std::cout << "]"<< std::endl;
		}
		for (unsigned i = 0; i < numThreads; ++i) {
			void *result;
			auto r = pthread_join(threads[i], &result);
			(void) r;
			assert(r == 0);
		}
	}
};

}

int main(int argc, char const * const argv[]) {
	sakura_Status result = sakura_Initialize(nullptr, nullptr);
	{
		Bench<InterpolateXBench<1000000UL> > bench(8, 4, 20000ULL);
		bench.run();
	}
	sakura_CleanUp();
	return 0;
}
