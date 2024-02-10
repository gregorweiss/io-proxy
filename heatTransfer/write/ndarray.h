/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * ndarray.h
 *
 *  Created on: Feb 2024
 *      Author: Gregor Weiss
 */
#ifndef NDARRAY_H
#define NDARRAY_H

#include <cstdarg>     // variadic function
#include <numeric>     // std::accumulate
#include <functional>  // std::multiplies
#include <memory>      // std::unique_ptr
#include <utility>     // std::integer_sequence
#include <type_traits> // std::common_type

//template<typename T, typename... Types>
//concept is_all_same = (... && std::is_same<T, Types>::value);

std::vector<size_t> 
calculate_offsets( std::initializer_list<size_t> dimensions );

template<typename _Tp>
struct ndarray 
{
    using value_type             = _Tp ;
    using pointer                = value_type*;
    using const_pointer          = const value_type*;
    using reference              = value_type&;
    using const_reference        = const value_type&;
    using iterator               = value_type*;
    using const_iterator         = const value_type*;
    using size_type              = std::size_t;
    using difference_type        = std::ptrdiff_t;
    using reverse_iterator       = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    ndarray() = default;

    ndarray(std::initializer_list<size_type> dimensions, value_type value)
      : _size{ std::accumulate( dimensions.begin(), dimensions.end(), size_type(1),
                                std::multiplies<size_type>() ) }
      , _dims( dimensions )
      , _offsets( calculate_offsets(dimensions) )
      , _array( new value_type[ _size ] )
    { std::fill_n(_array.get(), _size, value); }

    ndarray(ndarray<value_type>& other) {
        other.swap(*this);
    }

    ndarray<value_type>&
    operator=(ndarray<value_type>& other) {
        other.swap(*this);
        return *this;
    }

    ~ndarray() = default;

    void swap(ndarray<value_type>& other) noexcept {
        using std::swap;
        swap(this->_size, other._size);
        swap(this->_dims, other._dims);
        swap(this->_offsets, other._offsets);
        swap(this->_array, other._array);
    }

    reference
    operator()( std::convertible_to<size_type> auto  && ... indices )
    { return _array.get()[ position(size_type(indices)...) ]; }

    const_reference
    operator()( std::convertible_to<size_type> auto  && ... indices ) const
    { return _array.get()[ position(size_type(indices)...) ]; }

  private:

    size_type _size{};
    std::vector<size_type> _dims{};
    std::vector<size_type> _offsets{};
    std::unique_ptr<value_type> _array{};

    template<typename ... IndexType>
    size_type position( IndexType ... indices ) const {
        size_t dim = 0;
        size_t pos = 0;
        for ( auto&& index : { indices... } ) {
            pos += index * _offsets[dim];
            ++dim;
        }
        return pos;
    }
};

template<typename _Tp>
void swap( ndarray<_Tp>& left, ndarray<_Tp>& right ) noexcept
{ left.swap(right); }

#endif /* NDARRAY_H_ */