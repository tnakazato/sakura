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
#include <numeric>
#include <chrono>
#include <cassert>
#include <cstdlib>
#include <cmath>
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
struct StatBench {
	struct WSType {
		static constexpr size_t N = SIZE;
		std::unique_ptr<std::array<float, N>,
				decltype(&free_obj<std::array<float, N> >)> data;
		std::unique_ptr<std::array<bool, N>,
				decltype(&free_obj<std::array<bool, N> >)> is_valid;
		WSType() :
				data(alloc_obj<std::array<float, N> >(),
						free_obj<std::array<float, N> >), is_valid(
						alloc_obj<std::array<bool, N> >(),
						free_obj<std::array<bool, N> >) {

		}
	};

	static WSType *createWS(unsigned id) {
		WSType *ws = new WSType();
		ws->data->fill(0.f);
		ws->is_valid->fill(true);
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
		LIBSAKURA_SYMBOL(StatisticsResultFloat) result;
		LIBSAKURA_SYMBOL(Status) status =
		LIBSAKURA_SYMBOL(ComputeStatisticsFloat)(ws->data->size(),
				ws->data->data(), ws->is_valid->data(), &result);
		assert(status == LIBSAKURA_SYMBOL(Status_kOK));
	}
};

template<size_t SIZE>
struct LSQFitBench {
	static constexpr uint16_t ORDER = 14U;
	static constexpr size_t NUM_COEFF = ORDER + 1;
	struct WSType {
		static constexpr size_t N = SIZE;
		std::unique_ptr<std::array<float, N>,
				decltype(&free_obj<std::array<float, N> >)> data;
		std::unique_ptr<std::array<bool, N>,
				decltype(&free_obj<std::array<bool, N> >)> is_valid;
		std::unique_ptr<std::array<double, NUM_COEFF>,
				decltype(&free_obj<std::array<double, NUM_COEFF> >)> coeff;
		std::unique_ptr<std::array<bool, N>,
				decltype(&free_obj<std::array<bool, N> >)> is_valid_out;
		WSType() :
				data(alloc_obj<std::array<float, N> >(),
						free_obj<std::array<float, N> >), is_valid(
						alloc_obj<std::array<bool, N> >(),
						free_obj<std::array<bool, N> >), coeff(
						alloc_obj<std::array<double, NUM_COEFF> >(),
						free_obj<std::array<double, NUM_COEFF> >), is_valid_out(
						alloc_obj<std::array<bool, N> >(),
						free_obj<std::array<bool, N> >) {
		}
	};
	static WSType *createWS(unsigned id) {
		WSType *ws = new WSType();
		constexpr float period1 = WSType::N;
		constexpr float amplitude1 = 5.0f;

		constexpr float period2 = 0.695f * WSType::N;
		constexpr float amplitude2 = 5.0f;
		constexpr float pi = M_PI;
		for (size_t x = 0; x < WSType::N; ++x) {
			ws->data->data()[x] = amplitude1
					* std::sin(2.0f * pi * static_cast<float>(x) / period1)
					+ amplitude2
							* std::cos(
									2.0f * pi * static_cast<float>(x)
											/ period2);
		}
		constexpr size_t edge = static_cast<size_t>(0.1 * WSType::N);
		assert(0 < edge && edge < WSType::N);
		float y0 = std::accumulate(ws->data->data(), &ws->data->data()[edge],
				0.f) / static_cast<float>(edge);
		assert(edge <= WSType::N && 0 < WSType::N);
		float y1 = std::accumulate(&ws->data->data()[WSType::N - edge],
				&ws->data->data()[WSType::N], 0.f) / static_cast<float>(edge);
		float k = (y1 - y0) / static_cast<float>(WSType::N - 1 - 0);
		for (size_t x = 0; x < WSType::N; ++x) {
			ws->data->data()[x] -= y0 + k * static_cast<float>(x);
		}
		ws->is_valid->fill(true);
		return ws;
	}

	static void destroyWS(unsigned id, WSType *ws) {
		delete ws;
	}

	struct ThreadLocalType { // typical member is thread local context object
		unsigned thread_id;
		std::unique_ptr<sakura_LSQFitContextFloat,
				decltype(&LIBSAKURA_SYMBOL(DestroyLSQFitContextFloat))> context;
		ThreadLocalType() :
				thread_id(0), context(nullptr,
						LIBSAKURA_SYMBOL(DestroyLSQFitContextFloat)) {
		}
	};

	static ThreadLocalType *createThreadLocal(unsigned tread_id) {
		ThreadLocalType *tl = new ThreadLocalType;
		tl->thread_id = tread_id;
		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context = nullptr;
		auto status = LIBSAKURA_SYMBOL(CreateLSQFitContextPolynomialFloat)(
				LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), ORDER, WSType::N,
				&context);
		assert(status == LIBSAKURA_SYMBOL(Status_kOK));
		assert(context != nullptr);
		tl->context.reset(context);
		return tl;
	}

	static void destroyThreadLocal(unsigned tread_id,
			ThreadLocalType *threadLocal) {
		assert(threadLocal->context != nullptr);
	}

	static void run(ThreadLocalType *tl, WSType *ws, unsigned id) {
		LIBSAKURA_SYMBOL(LSQFitStatus) result;
		float rms = 0.f;
		LIBSAKURA_SYMBOL(Status) status =
		LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(tl->context.get(), ORDER,
				WSType::N, ws->data->data(), ws->is_valid->data(), 0.0001f, 1U,
				NUM_COEFF, ws->coeff->data(), nullptr, nullptr,
				ws->is_valid_out->data(), &rms, &result);
		assert(status == LIBSAKURA_SYMBOL(Status_kOK));
		assert(result == LIBSAKURA_SYMBOL(LSQFitStatus_kOK));
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

	if (false) {
		Bench<StatBench<1000000UL> > bench(8, 4, 20000ULL);
		bench.run();
	} else if (true) {
		Bench<LSQFitBench<1000000UL> > bench(8, 1, 20000ULL);
		bench.run();
	}
	sakura_CleanUp();
	return 0;
}
