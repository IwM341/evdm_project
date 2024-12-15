#ifndef PROGRESS_BAR_HPP
#define PROGRESS_BAR_HPP

#include <functional>
#include <atomic>
#include <mutex>
#ifdef _OPENMP
#include <omp.h>
#endif
namespace evdm {

	struct omp_locker {
#ifdef _OPENMP
		omp_lock_t m_lock;
		omp_locker() {
			omp_init_lock(&m_lock);
		}
		void lock() {
			omp_set_lock(&m_lock);
		}
		void unlock() {
			omp_unset_lock(&m_lock);
		}
#else
		void lock() {}
		void unlock() {}
#endif // _OPENMP
	};

	template <typename IntType = int>
	struct progress_omp_function {
		bool null_func;
		std::function<void(IntType, IntType)> m_next_func;
		progress_omp_function() :
			null_func(true), m_next_func([](auto, auto) {}) {}
		
		template <typename FuncType>
		progress_omp_function(FuncType m_next_func) :
			null_func(false), m_next_func(std::move(m_next_func)) {}
		inline void operator ()(IntType m_current, IntType m_full) {
			if (!null_func) {
				m_next_func(m_current, m_full);
			}
		}
	};

	template <typename IntType = int>
	struct progress_omp_bar {
		size_t full;
		size_t lapse;
		size_t m_current;
		omp_locker m_lock;
		progress_omp_function<IntType> m_update;
		

		progress_omp_bar(
			progress_omp_function<IntType> m_update = 
			progress_omp_function<IntType>(),
			size_t full = 100, size_t lapse = 1) :
			full(full), lapse(lapse), m_current(0),
			m_update(m_update) 
		{
			m_update(0, full);
		}

		inline void next() {
			if (!m_update.null_func) {
				m_lock.lock();
				size_t tmp = ++m_current;
				if (tmp % lapse == 0 || tmp == full) {
					m_update(tmp, full);
				}
				m_lock.unlock();
			}
		}
	};
};

#endif//PROGRESS_BAR_HPP