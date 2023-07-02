#ifndef VEQF_FILTER_BITSUTIL_H_
#define VEQF_FILTER_BITSUTIL_H_

#include <cstdint>

namespace veqf
{
inline uint64_t upperpower2(uint64_t x)
{
    x--;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x |= x >> 32;
    x++;
    return x;
}
}

#endif