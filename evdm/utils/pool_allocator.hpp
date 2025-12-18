#pragma once
#include <list>
#include <cstdlib>
#include <cstddef>
namespace evdm {
 

template <typename T>
class PoolAllocator {
private:
    static_assert(alignof(T) <= alignof(std::max_align_t), "PoolAllocator<T>: alignof(T) > alignof(std::max_align_t)");

    struct Block {
        void* ptr;
        size_t capacity;
        size_t used;
        
        Block(void* p, size_t cap) : ptr(p), capacity(cap), used(0) {}
    };

    struct FBlock { 
        void* ptr;
        size_t capacity;
        
        FBlock(void* p, size_t cap) : ptr(p), capacity(cap) {}
    };

    // Non-full blocks (used < capacity)
    std::list<Block> nonFullBlocks;
    
    // Full blocks (used == capacity)
    std::list<FBlock> fullBlocks;
    
    size_t minSize;
    
    void* allocateBlock(size_t size) {
        void* ptr = std::malloc (size);
        if (!ptr) {
            throw std::bad_alloc();
        }
        return ptr;
    }

public:
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;
    using reference = T&;
    using const_reference = const T&;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    template<typename U>
    struct rebind {
        using other = PoolAllocator<U>;
    };

    // Constructor with dynamic min_size parameter
    explicit PoolAllocator(size_t size_min = 1024) : minSize(size_min) {
        if (minSize < sizeof(T)) {
            minSize = sizeof(T);
        }
    }

    // Copy constructor
    PoolAllocator(const PoolAllocator& other) noexcept 
        : minSize(other.minSize) {
        // Don't copy blocks - each allocator manages its own memory
    }

    template<typename U>
    PoolAllocator(const PoolAllocator<U>& other) noexcept 
        : minSize(other.getMinSize()) {
    }

    ~PoolAllocator() {
        // Free all memory in destructor only
        for (auto& block : nonFullBlocks) {
            if (block.ptr) {
                std::free(block.ptr);
            }
        }
        for (auto& block : fullBlocks) {
            if (block.ptr) {
                std::free(block.ptr);
            }
        }
    }
    void flushExceptOne() {
        if (nonFullBlocks.empty()) {
            bool empty = true;
            for (FBlock & block : fullBlocks.end()) {
                if (empty && block.capacity >= minSize) {
                    nonFullBlocks.push_back(block.ptr, block.capacity);
                    empty = false;
                }
                else {
                    if (block.ptr) {
                        std::free(block.ptr);
                    }
                }
            }
        }
        fullBlocks.clear();
        while (nonFullBlocks.size() > 1) {
            Block & block = nonFullBlocks.back();
            nonFullBlocks.pop_back();
            if (block.ptr) {
                std::free(block.ptr);
            }
            nonFullBlocks.pop_back();
        }
    }
    pointer allocate(size_type n) {
        if (n == 0) {
            return nullptr;
        }

        size_type required_size = n * sizeof(T);
        
        if (required_size > minSize) {
            // Allocate new block for large requests
            void* ptr = allocateBlock(required_size);
            fullBlocks.emplace_back(ptr, required_size);
            return static_cast<pointer>(ptr);
        } else {
            // Try to find a non-full block
            for (auto it = nonFullBlocks.begin(); it != nonFullBlocks.end(); ++it) {
                if (it->capacity - it->used >= required_size) {
                    void* result = static_cast<char*>(it->ptr) + it->used;
                    it->used += required_size;
                    
                    // If block becomes full, move it to fullBlocks
                    if (it->used == it->capacity) {
                        fullBlocks.emplace_back(it->ptr, it->capacity);
                        nonFullBlocks.erase(it);
                    }
                    
                    return static_cast<pointer>(result);
                }
            }
            
            // No suitable non-full block found, allocate new one
            void* ptr = allocateBlock(minSize);
            void* result = ptr;
            size_t remaining = minSize - required_size;
            
            if (remaining > 0) {
                // Add as non-full block
                Block newBlock(ptr, minSize);
                newBlock.used = required_size;
                nonFullBlocks.push_back(newBlock);
            } else {
                // Block is exactly used up
                fullBlocks.emplace_back(ptr, minSize);
            }
            
            return static_cast<pointer>(result);
        }
    }

    void deallocate(pointer p, size_type n) noexcept {
        // Memory is only freed in destructor, so do nothing here
        (void)p; // Suppress unused parameter warning
        (void)n; // Suppress unused parameter warning
    }

    template<typename U, typename... Args>
    void construct(U* p, Args&&... args) {
        new (p) U(std::forward<Args>(args)...);
    }

    template<typename U>
    void destroy(U* p) {
        p->~U();
    }

    // Utility functions
    size_t getMinSize() const noexcept { return minSize; }
    
    size_t getNonFullBlockCount() const noexcept { 
        return nonFullBlocks.size(); 
    }
    
    size_t getFullBlockCount() const noexcept { 
        return fullBlocks.size(); 
    }
    
    size_t getTotalAllocatedMemory() const noexcept {
        size_t total = 0;
        for (const auto& block : nonFullBlocks) {
            total += block.capacity;
        }
        for (const auto& block : fullBlocks) {
            total += block.capacity;
        }
        return total;
    }

    // Comparison operators (required by allocator concept)
    bool operator==(const PoolAllocator& other) const noexcept {
        return this == &other;
    }

    bool operator!=(const PoolAllocator& other) const noexcept {
        return !(*this == other);
    }
};
       
}; // namespace evdm
