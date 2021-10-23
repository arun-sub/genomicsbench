#ifndef ALIGNED_ALLOCATOR_H
#define ALIGNED_ALLOCATOR_H

#include <stdexcept>
#include <xmmintrin.h>

/**
 * Allocator for aligned data.
 *
 * Modified from the Mallocator from Stephan T. Lavavej.
 * <http://blogs.msdn.com/b/vcblog/archive/2008/08/28/the-mallocator.aspx>
 * Modified from Gist version by Don
 * https://gist.githubusercontent.com/donny-dont/1471329/raw/8f063f5f4326d301b14fe6b781495be54ca48941/aligned_allocator.cpp
 */

// T is the element type
// Alighment is the alignment, in bytes
// Offset \in [0,Alighment-1] is the offset to the first aligned byte
template <typename T, std::size_t Alignment, std::size_t Offset = 0>
class Aligned_allocator
{
	public:

		// The following will be the same for virtually all allocators.
		typedef T * pointer;
		typedef const T * const_pointer;
		typedef T& reference;
		typedef const T& const_reference;
		typedef T value_type;
		typedef std::size_t size_type;
		typedef std::ptrdiff_t difference_type;

		T * address(T& r) const {
			return &r;
		}

		const T * address(const T& s) const {
			return &s;
		}

		std::size_t max_size() const {
			// The following has been carefully written to be independent of
			// the definition of size_t and to avoid signed/unsigned warnings.
			return (static_cast<std::size_t>(0) - static_cast<std::size_t>(1)) / sizeof(T);
		}


		// The following must be the same for all allocators.
		template <typename U> struct rebind {
			typedef Aligned_allocator<U, Alignment, Offset> other;
		};

		bool operator!=(const Aligned_allocator& other) const {
			return !(*this == other);
		}

		void construct(T * const p, const T& t) const {
			void * const pv = reinterpret_cast<void *>(p);
			new (pv) T(t);
		}

		void destroy(T * const p) const {
			p->~T();
		}

		// Returns true if and only if storage allocated from *this
		// can be deallocated from other, and vice versa.
		// Always returns true for stateless allocators.
		bool operator==(const Aligned_allocator& other) const {
			return true;
		}


		// Default constructor, copy constructor, rebinding constructor, and destructor.
		// Empty for stateless allocators.
		Aligned_allocator() { }

		Aligned_allocator(const Aligned_allocator&) { }

		template <typename U> Aligned_allocator(const Aligned_allocator<U, Alignment>&) { }

		~Aligned_allocator() { }


		// The following will be different for each allocator.
		T * allocate(const std::size_t n) const {
			// The return value of allocate(0) is unspecified.
			// Mallocator returns NULL in order to avoid depending
			// on malloc(0)'s implementation-defined behavior
			// (the implementation can define malloc(0) to return NULL,
			// in which case the bad_alloc check below would fire).
			// All allocators can return NULL in this case.
			if (n == 0) {
				return NULL;
			}

			// All allocators should contain an integer overflow check.
			// The Standardization Committee recommends that std::length_error
			// be thrown in the case of integer overflow.
			if (n > max_size()) {
				throw std::length_error("Aligned_allocator<T>::allocate() - Integer overflow.");
			}

			// Mallocator wraps malloc().
			const std::size_t o = Offset > 0? Alignment - Offset: 0;
			void * const pv = _mm_malloc(n * sizeof(T) + o, Alignment);
			// Allocators should throw std::bad_alloc in the case of memory allocation failure.
			if (pv == NULL) {
				throw std::bad_alloc();
			}
#ifdef TRACE_ALLOC
std::cerr << "Alloc " << n << '\n';
#endif
			return reinterpret_cast<T *>(reinterpret_cast<char *>(pv)+o);
		}

		void deallocate(T * const p, const std::size_t n) const {
            (void) n;
			const std::size_t o = Offset > 0? Alignment - Offset: 0;
#ifdef TRACE_ALLOC
std::cerr << "Free " << n << '\n';
#endif
			_mm_free(reinterpret_cast<char *>(p)-o);
		}


		// The following will be the same for all allocators that ignore hints.
		template <typename U>
		T * allocate(const std::size_t n, const U * /* const hint */) const {
			return allocate(n);
		}
		// Allocators are not required to be assignable, so
		// all allocators should have a private unimplemented
		// assignment operator. Note that this will trigger the
		// off-by-default (enabled under /Wall) warning C4626
		// "assignment operator could not be generated because a
		// base class assignment operator is inaccessible" within
		// the STL headers, but that warning is useless.
	private:
		Aligned_allocator& operator=(const Aligned_allocator&);
};

#endif

#if 0
int main()
{
	std::vector<float, Aligned_allocator<float, 16> > aligned_16_0;
	std::vector<float, Aligned_allocator<float, 16, 4> > aligned_16_4;
	std::vector<float, Aligned_allocator<float, 16, 8> > aligned_16_8;
	std::vector<float, Aligned_allocator<float, 16, 12> > aligned_16_12;
    for (int i = 0; i < 1024*1024; i++) {
        aligned_16_0.push_back(float(i));
        aligned_16_4.push_back(float(i));
        aligned_16_8.push_back(float(i));
        aligned_16_12.push_back(float(i));
    }
    std::cerr << reinterpret_cast<unsigned long int>(&aligned_16_0[0]) % 16 << '\n';
    std::cerr << reinterpret_cast<unsigned long int>(&aligned_16_4[1]) % 16 << '\n';
    std::cerr << reinterpret_cast<unsigned long int>(&aligned_16_8[2]) % 16 << '\n';
    std::cerr << reinterpret_cast<unsigned long int>(&aligned_16_12[3]) % 16 << '\n';
    return 0;
}
#endif


