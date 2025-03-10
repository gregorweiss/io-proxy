
#include "ndarray.h"

std::vector<size_t>
calculate_offsets( std::initializer_list<size_t> dimensions ) {
    std::vector<size_t> ret( dimensions.size(), 1 );
    auto iter = dimensions.end()-1;
    for ( size_t i = dimensions.size() - 1;
          i > 0; --iter, --i ) {
        ret[i-1] = ret[i] * (*iter);
    }
    return ret;
}

std::vector<size_t>
calculate_row_major_offsets( std::initializer_list<size_t> dimensions ) {
    std::vector<size_t> ret( dimensions.size(), 1 );
    size_t i = 1;
    for ( auto iter = dimensions.begin();
          iter < dimensions.end() - 1;
          ++iter, ++i ) {
        ret[i] = ret[i-1] * (*iter);
    }
    return ret;
}