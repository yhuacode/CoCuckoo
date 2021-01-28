#ifndef _SPINLOCK_HPP
#define _SPINLOCK_HPP

#include <stdio.h>
#include <atomic>


 class spinlock 
 {
        std::atomic_flag lock_;
    public:
        spinlock() {
			lock_.clear();
        }

        inline void lock() {
            while (lock_.test_and_set(std::memory_order_acquire));
        }

        inline void unlock() {
            lock_.clear(std::memory_order_release);
        }

        inline bool try_lock() {
            return !lock_.test_and_set(std::memory_order_acquire);
        }

} __attribute__((aligned(64)));

#endif
