//
//  blocklist.h
//  DAG
//
//  Created by Troels Tang Karlsen on 02/09/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//

#ifndef blocklist_h
#define blocklist_h

#include <array>
#include <list>
#include <iterator>

using namespace std;

template <class T, size_t block_size>
class blockList
{
    // The data type we store data in, list of arrays of numbers(T) with specified block_size
    list<array<T, block_size>> data;
    
    // Define iterators of the same type as first array in data list
    using list_iter = decltype(data.begin());
    list_iter  cur_block; // Current block
    list_iter  last_block; // Last block in list
    
    // Define iterators of the same type as first element in last array in data
    using block_iter = decltype(data.back().begin());
    block_iter  next_space; // next available space
    block_iter  last_space; // last available space
    
    // Marks to rewind back to
    list_iter  marked_block;
    block_iter marked_space;
    
    // Add new array to list and update pointers
    void newBlock(){
        // Create new array at end of list
        data.emplace_back();
        
        // Set iterators to newly added block
        cur_block = last_block = prev(data.end()); // prev is iterator pointing to element before input
        next_space = cur_block->begin();
        last_space = cur_block->end();
    }
    
    void nextBlock(){
        if( cur_block == last_block){
            newBlock();
        }
        else {
            ++cur_block;
            next_space = cur_block->begin();
            last_space = cur_block->end();
        }
    }
    
public:
    // Constructor: create first array in list
    blockList() {
        newBlock();
    }
    
    // Empty list and construct new first array in list
    void clear(){
        data.clear();
        newBlock();
    }
    
    // Go back and overwrite without deleting contents of blocks
    void rewind(){
        cur_block = data.begin();
        next_space = cur_block->begin();
        last_space = cur_block->end();
    }
    
    // Save this point
    void setMark(){
        marked_block = cur_block;
        marked_space = next_space;
    }
    
    // Rewind to marked point
    void rewindToMark(){
        cur_block = marked_block;
        next_space = marked_space;
        last_space = cur_block->end();
    }
    
    // Overload emplace_back
    T* emplace_back() {
        // Create a new block if needed
        if (next_space == last_space) {
            nextBlock();
        }
        
        auto old_next_space = next_space;
        
        ++next_space;
        
        // Return pointer to what used to be last element in list
        return &*old_next_space;
    }
    
    // emplace_back for multiple blocks, if number of blocks (n) is known at compile time
    template <size_t n>
    T* emplace_back_multi() {
        
        // Create a new block if needed
        if (distance(next_space, last_space) < n) {
            nextBlock();
        }
        
        auto old_next_space = next_space;
        
        next_space += n;
        
        // Return pointer to what used to be last element in list
        return &*old_next_space;
    }
    
    // emplace_back for multiple blocks, if number of blocks (n) is not known at compile time
    T* emplace_back_multi(size_t n) {
        
        // Create a new block if needed
        if (distance(next_space, last_space) < n) {
            nextBlock();
        }
        
        auto old_next_space = next_space;
        
        next_space += n;
        
        // Return pointer to what used to be last element in list
        return &*old_next_space;
    }
    
    template<typename ...Args>
    T* emplace_back(Args&& ...args){
        // Create a new block if needed
        if (next_space == last_space) {
            nextBlock();
        }
        
        T* emplaced = new(&*next_space) T(forward<Args>(args)...);
        
        ++next_space;
        
        // Return pointer to what used to be last element in list
        return emplaced;
    }
    
    // Overload memset function, to set all elements of arrays in the data-list to value (default 0)
    void memset(unsigned char value = 0){
        for (auto& array : data){
            memset(&array[0], value, block_size * sizeof(T));
        }
    }
    
    // Defined inside blocklist class
    class iterator {
        list_iter   cur_block; // Current block
        block_iter  cur_space; // Current space
        
        // Define iterators of the same type as first element in last array in data
        block_iter  first_space; // first space in block
        block_iter  last_space;  // last available space in block
        
    public:
        using difference_type = ptrdiff_t;
        using reference = T&;
        using pointer = T*;
        using value_type = T;
        using iterator_category = bidirectional_iterator_tag;
        
        // Constructors
        iterator() {}
        
        iterator(list_iter cb, block_iter cs, block_iter fs, block_iter ls) :
        cur_block(cb), cur_space(cs), first_space(fs), last_space(ls) { }
        
        // Operator overloading
        iterator& operator++(){
            ++cur_space;
            
            if(cur_space == last_space){
                ++cur_block;
                
                first_space = cur_block->begin();
                last_space = cur_block->end();
                cur_space = first_space;
            }
            return *this;
        }
        iterator& operator--(){
            if(cur_space == last_space){
                --cur_block;
                
                first_space = cur_block->begin();
                last_space = cur_block->end();
                cur_space = last_space;
            }
            --cur_space;
            
            return *this;
        }
        
        T& operator*(){
            return *cur_space;
        }
        
        const T& operator*() const {
            return * cur_space;
        }
        T* operator->(){
            return &*cur_space;
        }
        const T* operator->() const {
            return &*cur_space;
        }
        bool operator == (const iterator& right) const
        {
            return (cur_block == right.cur_block && cur_space == right.cur_space);
        }
        bool operator != (const iterator& right) const
        {
            return (cur_block != right.cur_block || cur_space != right.cur_space);
        }
        
    };
    
    iterator begin()
    {
        return iterator(data.begin(), data.begin()->begin(), data.begin()->begin(), data.begin()->end() );
    }
    // Returns iterator to next available slot for storage of an element
    iterator end()
    {
        //auto last_block = prev(data.end()); // Why unused?
        return iterator(cur_block, next_space, cur_block->begin(), cur_block->end() );
    }
    
    iterator mark()
    {
        return iterator(marked_block, marked_space, marked_block->begin(), marked_block->end() );
    }
    
    // Find element by searching from end to beginning
    iterator find(const T* const element)
    {
        iterator it = end();
        iterator b = begin();
        
        while(it != b)
        {
            --it;
            if(&*it == element) return it;
        }
        
        if(&*it == element) return it;
        
        return end();
    }
};





#endif /* blocklist_h */
