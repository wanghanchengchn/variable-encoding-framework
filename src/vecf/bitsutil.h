#ifndef VECF_BITSUTIL_H_
#define VECF_BITSUTIL_H_

#include <cstdint>

namespace vecf
{
#define haszero4(x) (((x)-0x1111ULL) & (~(x)) & 0x8888ULL)
#define hasvalue4(x, n) (haszero4((x) ^ (0x1111ULL * (n))))

#define haszero8(x) (((x)-0x01010101ULL) & (~(x)) & 0x80808080ULL)
#define hasvalue8(x, n) (haszero8((x) ^ (0x01010101ULL * (n))))

#define haszero12(x) (((x)-0x001001001001ULL) & (~(x)) & 0x800800800800ULL)
#define hasvalue12(x, n) (haszero12((x) ^ (0x001001001001ULL * (n))))

#define haszero16(x) \
    (((x)-0x0001000100010001ULL) & (~(x)) & 0x8000800080008000ULL)
#define hasvalue16(x, n) (haszero16((x) ^ (0x0001000100010001ULL * (n))))

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

template <uint32_t tag_length, typename T>
inline T MaskedTag(const T tag)
{
    // static_assert(std::is_same_v<T, uint32_t> || std::is_same_v<T, uint64_t>,
    // "unsupported types");
    return tag & ((1ULL << tag_length) - 1);
}

template <typename T>
inline T MaskedTag(const T tag, uint32_t tag_length)
{
    // static_assert(std::is_same_v<T, uint32_t> || std::is_same_v<T, uint64_t>,
    // "unsupported types");
    return tag & ((1ULL << tag_length) - 1);
}

}

#endif
