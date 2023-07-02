#ifndef VECF_SINGLE_TABLE_H_
#define VECF_SINGLE_TABLE_H_

#include <immintrin.h>

#include <cassert>
#include <cstring>
#include <memory>
#include <type_traits>

#include "vecf/bitsutil.h"

namespace vecf
{
// common interface of SingleTable
template <uint64_t kBitsPerItem>
class BaseSingleTable
{
  public:
    constexpr static uint32_t kFourSlotTagLen = kBitsPerItem;

    uint64_t NumBuckets() const
    {
        return num_buckets_;
    }

    size_t SizeInBytes() const
    {
        return kBytesPerBucket * num_buckets_;
    }

    size_t SizeInTags() const
    {
        return kTagsPerBucket * num_buckets_;
    }

  protected:
    // derived class is responsible to initialize `buckets_`
    explicit BaseSingleTable(const uint64_t num_buckets)
        : num_buckets_(num_buckets),
          buckets_(new Bucket[num_buckets + kPaddingBuckets])
    {
    }

    constexpr static uint64_t kTagsPerBucket = 4;
    constexpr static uint64_t kBytesPerBucket =
        (kBitsPerItem * kTagsPerBucket + 7) >> 3;
    constexpr static uint64_t kPaddingBuckets =
        ((((kBytesPerBucket + 7) / 8) * 8) - 1) / kBytesPerBucket;

    struct Bucket
    {
        char bits_[kBytesPerBucket];
    } __attribute__((__packed__));

    uint64_t num_buckets_;
    std::unique_ptr<Bucket[]> buckets_;

    template <uint64_t index, uint32_t tag_length, typename T>
    static inline T BucketTag(const T b)
    {
        static_assert(std::is_same_v<T, uint32_t> || std::is_same_v<T, uint64_t>,
                      "unsupported types");
        static_assert(index < kTagsPerBucket, "bucket index overflow");
        if constexpr (tag_length == 8)
        {
            static_assert(std::is_same_v<T, uint32_t>);
            return ((const uint8_t *)&b)[index]; // aligned byte access
        }
        else
        {
            return (b >> (index * tag_length)) & ((1ull << tag_length) - 1);
        }
    }

    template <uint32_t tag_length, typename T>
    static inline T BucketTag(const T b, uint64_t index)
    {
        static_assert(std::is_same_v<T, uint32_t> || std::is_same_v<T, uint64_t>,
                      "unsupported types");
        assert(index < kTagsPerBucket);
        if constexpr (tag_length == 8)
        {
            static_assert(std::is_same_v<T, uint32_t>);
            return ((const uint8_t *)&b)[index]; // aligned byte access
        }
        else
        {
            return (b >> (index * tag_length)) & ((1ull << tag_length) - 1);
        }
    }
};

template <size_t bits_per_tag>
class SingleTable
{
};

#define haszero29(x) (((x)-0x00000001ULL) & (~(x)) & 0x10000000ULL)
#define hasvalue29(x, n) (haszero29((x) ^ (0x00000001ULL * (n))))
#define haszero14(x) (((x)-0x0004001ULL) & (~(x)) & 0x8002000ULL)
#define hasvalue14(x, n) (haszero14((x) ^ (0x0004001ULL * (n))))
#define haszero9(x) (((x)-0x0040201ULL) & (~(x)) & 0x4020100ULL)
#define hasvalue9(x, n) (haszero9((x) ^ (0x0040201ULL * (n))))

template <>
class SingleTable<8> : public BaseSingleTable<8>
{
  public:
    explicit SingleTable(uint64_t num_buckets)
        : BaseSingleTable(num_buckets)
    {
        static_assert(kBytesPerBucket == 4, "ctor only work on this case");

        uint32_t *p = reinterpret_cast<uint32_t *>(buckets_.get());
        for (uint64_t i = 0; i < num_buckets_; ++i, ++p)
        {
            *p = kZeroSlotFlag;
        }
    }

    bool FindTagInBucket(const uint64_t i, const uint32_t unmasked_tag) const
    {
        const uint32_t bucket{*reinterpret_cast<uint32_t *>(buckets_[i].bits_)};
        const uint32_t flag{bucket & kFlagBitsMask};

        switch (flag)
        {
        case kZeroSlotFlag: {
            return false;
        }
        case kOneSlotFlag: {
            const uint32_t tags{_pext_u32(bucket, kTagBitsMask)};
            return hasvalue29(tags, MaskedTag<kOneSlotTagLen>(unmasked_tag)) ||
                   hasvalue29(tags, MaskedTag<kTwoSlotTagLen>(unmasked_tag)) ||
                   hasvalue29(tags, MaskedTag<kThreeSlotTagLen>(unmasked_tag)) ||
                   hasvalue29(tags, MaskedTag<kFourSlotTagLen>(unmasked_tag));
        }
        case kTwoSlotFlag: {
            const uint32_t tags{_pext_u32(bucket, kTagBitsMask)};
            return hasvalue14(tags, MaskedTag<kTwoSlotTagLen>(unmasked_tag)) ||
                   hasvalue14(tags, MaskedTag<kThreeSlotTagLen>(unmasked_tag)) ||
                   hasvalue14(tags, MaskedTag<kFourSlotTagLen>(unmasked_tag));
        }
        case kThreeSlotFlag: {
            const uint32_t tags{_pext_u32(bucket, kTagBitsMask)};
            return hasvalue9(tags, MaskedTag<kThreeSlotTagLen>(unmasked_tag)) ||
                   hasvalue9(tags, MaskedTag<kFourSlotTagLen>(unmasked_tag));
        }
        default: {
            return hasvalue8(bucket, MaskedTag<kFourSlotTagLen>(unmasked_tag));
        }
        }
    }

    inline bool InsertTagToBucket(const uint64_t i, const uint64_t tag,
                                  const bool kickout, uint64_t &oldtag)
    {
        uint32_t bucket{*reinterpret_cast<uint32_t *>(buckets_[i].bits_)};
        const uint32_t flag{bucket & kFlagBitsMask};

        switch (flag)
        {
        case kZeroSlotFlag: {
            bucket = _pdep_u32(tag, kTagBitsMask) | kOneSlotFlag;
            *reinterpret_cast<uint32_t *>(buckets_[i].bits_) = bucket;
            return true;
        }
        case kOneSlotFlag: {
            uint32_t tags{_pext_u32(bucket, 0x00003fff)}; // get one 14bit tags
            tags |= (MaskedTag<kTwoSlotTagLen>(tag) << kTwoSlotTagLen);
            *reinterpret_cast<uint32_t *>(buckets_[i].bits_) =
                _pdep_u32(tags, kTagBitsMask) | kTwoSlotFlag;
            return true;
        }
        case kTwoSlotFlag: {
            uint32_t tags{_pext_u32(bucket, 0x017f41ff)}; // get two 9bit tags
            tags |= MaskedTag<kThreeSlotTagLen>(tag) << (2 * kThreeSlotTagLen);
            *reinterpret_cast<uint32_t *>(buckets_[i].bits_) =
                _pdep_u32(tags, kTagBitsMask) | kThreeSlotFlag;
            return true;
        }
        case kThreeSlotFlag: {
            uint32_t tags{_pext_u32(bucket, 0x0f7b7eff)}; // get three 8bit tags
            tags |= (MaskedTag<kFourSlotTagLen>(tag) << (3 * kFourSlotTagLen));
            uint32_t tag3{BucketTag<3, kFourSlotTagLen>(tags)},
                tag2{BucketTag<2, kFourSlotTagLen>(tags)},
                tag1{BucketTag<1, kFourSlotTagLen>(tags)};
            if (tag3 > tag2 || tag2 > tag1)
            {
                if (tag3 > tag2)
                {
                    std::swap(tag3, tag2);
                }
                if (tag3 > tag1)
                {
                    std::swap(tag3, tag1);
                }
                if (tag2 > tag1)
                {
                    std::swap(tag2, tag1);
                }
                tags &= 0x000000ff;
                tags |= (tag3 << (3 * kFourSlotTagLen)) |
                        (tag2 << (2 * kFourSlotTagLen)) | (tag1 << kFourSlotTagLen);
            }
            *reinterpret_cast<uint32_t *>(buckets_[i].bits_) = tags;
            assert((BucketTag<3, kFourSlotTagLen>(tags) <=
                        BucketTag<2, kFourSlotTagLen>(tags) &&
                    BucketTag<2, kFourSlotTagLen>(tags) <=
                        BucketTag<1, kFourSlotTagLen>(tags)));
            return true;
        }
        default: {
            if (kickout)
            {
                static_assert(kFourSlotTagLen == 8,
                              "can not access bucket as uint8_t array");

                uint64_t r{rand() % kTagsPerBucket};
                oldtag = BucketTag<kFourSlotTagLen>(bucket, r);
                ((uint8_t *)(&bucket))[r] = MaskedTag<kFourSlotTagLen>(tag);
                uint32_t tag3{BucketTag<3, kFourSlotTagLen>(bucket)},
                    tag2{BucketTag<2, kFourSlotTagLen>(bucket)},
                    tag1{BucketTag<1, kFourSlotTagLen>(bucket)};
                if (tag3 > tag2 || tag2 > tag1)
                {
                    if (tag3 > tag2)
                    {
                        std::swap(tag3, tag2);
                    }
                    if (tag3 > tag1)
                    {
                        std::swap(tag3, tag1);
                    }
                    if (tag2 > tag1)
                    {
                        std::swap(tag2, tag1);
                    }
                    bucket &= 0x000000ff;
                    bucket |= (tag3 << (3 * kFourSlotTagLen)) |
                              (tag2 << (2 * kFourSlotTagLen)) |
                              (tag1 << kFourSlotTagLen);
                }
                *reinterpret_cast<uint32_t *>(buckets_[i].bits_) = bucket;
                assert((BucketTag<3, kFourSlotTagLen>(bucket) <=
                            BucketTag<2, kFourSlotTagLen>(bucket) &&
                        BucketTag<2, kFourSlotTagLen>(bucket) <=
                            BucketTag<1, kFourSlotTagLen>(bucket)));
            }
            return false;
        }
        }
    }

    void FindMaxMatchingTag(const uint64_t i, const uint32_t unmasked_tag,
                            uint64_t *max_bucket_idx,
                            uint32_t *max_tag_length) const
    {
        uint32_t bucket{*reinterpret_cast<uint32_t *>(buckets_[i].bits_)};
        const uint32_t flag{bucket & kFlagBitsMask};
        switch (flag)
        {
        case kZeroSlotFlag: {
            return;
        }
        case kOneSlotFlag: {
            const uint32_t tags{_pext_u32(bucket, kTagBitsMask)};
#define HELPER(tag_length)                                     \
    if (*max_tag_length >= tag_length)                         \
    {                                                          \
        return; /* Longer or equal matched tag found */        \
    }                                                          \
    if (hasvalue29(tags, MaskedTag<tag_length>(unmasked_tag))) \
    {                                                          \
        *max_bucket_idx = i;                                   \
        *max_tag_length = tag_length;                          \
        return;                                                \
    }

            HELPER(kOneSlotTagLen)
            HELPER(kTwoSlotTagLen)
            HELPER(kThreeSlotTagLen)
            HELPER(kFourSlotTagLen)
#undef HELPER
            return;
        }
        case kTwoSlotFlag: {
            const uint32_t tags{_pext_u32(bucket, kTagBitsMask)};
#define HELPER(tag_length)                                     \
    if (*max_tag_length >= tag_length)                         \
    {                                                          \
        return; /* Longer or equal matched tag found */        \
    }                                                          \
    if (hasvalue14(tags, MaskedTag<tag_length>(unmasked_tag))) \
    {                                                          \
        *max_bucket_idx = i;                                   \
        *max_tag_length = tag_length;                          \
        return;                                                \
    }

            HELPER(kTwoSlotTagLen)
            HELPER(kThreeSlotTagLen)
            HELPER(kFourSlotTagLen)
#undef HELPER
            return;
        }
        case kThreeSlotFlag: {
            const uint32_t tags{_pext_u32(bucket, kTagBitsMask)};
#define HELPER(tag_length)                                    \
    if (*max_tag_length >= tag_length)                        \
    {                                                         \
        return; /* Longer or equal matched tag found */       \
    }                                                         \
    if (hasvalue9(tags, MaskedTag<tag_length>(unmasked_tag))) \
    {                                                         \
        *max_bucket_idx = i;                                  \
        *max_tag_length = tag_length;                         \
        return;                                               \
    }

            HELPER(kThreeSlotTagLen)
            HELPER(kFourSlotTagLen)
#undef HELPER
            return;
        }
        default: {
            // encoding requirements
            assert((BucketTag<3, kFourSlotTagLen>(bucket) <=
                        BucketTag<2, kFourSlotTagLen>(bucket) &&
                    BucketTag<2, kFourSlotTagLen>(bucket) <=
                        BucketTag<1, kFourSlotTagLen>(bucket)));

            if (*max_tag_length >= kFourSlotTagLen)
            {
                return; // Longer or equal matched tag found
            }
            if (hasvalue8(bucket, MaskedTag<kFourSlotTagLen>(unmasked_tag)))
            {
                *max_bucket_idx = i;
                *max_tag_length = kFourSlotTagLen;
            }
            return;
        }
        }
    }

    void DeleteTagFromBucket(uint64_t bucket_idx, const uint32_t masked_tag)
    {
        uint32_t bucket{*reinterpret_cast<uint32_t *>(buckets_[bucket_idx].bits_)};
        const uint32_t flag{bucket & kFlagBitsMask};
        switch (flag)
        {
        case kZeroSlotFlag: {
            // never delete on empty bucket
            assert(false);
            __builtin_unreachable();
            return;
        }
        case kOneSlotFlag: {
            *reinterpret_cast<uint32_t *>(buckets_[bucket_idx].bits_) =
                kZeroSlotFlag;
            return;
        }
        case kTwoSlotFlag: {
            uint32_t masks[]{0x3f7f4000, 0x00003fff};
            const uint32_t tags{_pext_u32(bucket, kTagBitsMask)};
            for (uint32_t slot_idx = 0; slot_idx < 2; ++slot_idx)
            {
                if (BucketTag<kTwoSlotTagLen>(tags, slot_idx) == masked_tag)
                {
                    *reinterpret_cast<uint32_t *>(buckets_[bucket_idx].bits_) =
                        _pdep_u32(_pext_u32(bucket, masks[slot_idx]), kTagBitsMask) |
                        kOneSlotFlag;
                    return;
                }
            }
            assert(false);
            __builtin_unreachable();
        }
        case kThreeSlotFlag: {
            uint32_t masks[]{0x1f7f7e00, 0x1f7801ff, 0x000077fff};
            const uint32_t tags{_pext_u32(bucket, kTagBitsMask)};
            for (uint32_t slot_idx = 0; slot_idx < 3; ++slot_idx)
            {
                if (BucketTag<kThreeSlotTagLen>(tags, slot_idx) == masked_tag)
                {
                    *reinterpret_cast<uint32_t *>(buckets_[bucket_idx].bits_) =
                        _pdep_u32(_pext_u32(bucket, masks[slot_idx]), 0x017f41ff) |
                        kTwoSlotFlag;
                    return;
                }
            }
            assert(false);
            __builtin_unreachable();
        }
        default: {
            uint32_t masks[]{0xffffff00, 0xffff00ff, 0xff00ffff, 0x00ffffff};
            for (uint32_t slot_idx = 0; slot_idx < 4; ++slot_idx)
            {
                if (BucketTag<kFourSlotTagLen>(bucket, slot_idx) == masked_tag)
                {
                    *reinterpret_cast<uint32_t *>(buckets_[bucket_idx].bits_) =
                        _pdep_u32(_pext_u32(bucket, masks[slot_idx]), 0x0f7b7eff) |
                        kThreeSlotFlag;
                    return;
                }
            }
            assert(false);
            __builtin_unreachable();
        }
        }
    }

    void BucketCountStat(uint64_t *counter)
    {
        for (uint64_t i = 0; i < num_buckets_; ++i)
        {
            uint32_t bucket{*reinterpret_cast<uint32_t *>(buckets_[i].bits_)};
            const uint32_t flag{bucket & kFlagBitsMask};
            switch (flag)
            {
            case kZeroSlotFlag: {
                ++counter[0];
                break;
            }
            case kOneSlotFlag: {
                ++counter[1];
                break;
            }
            case kTwoSlotFlag: {
                ++counter[2];
                break;
            }
            case kThreeSlotFlag: {
                ++counter[3];
                break;
            }
            default: {
                ++counter[4];
                break;
            }
            }
        }
    }

  private:
    constexpr static uint32_t kFlagBitsMask = 0x80808000;
    constexpr static uint32_t kTagBitsMask = ~kFlagBitsMask;

    // other case is four slots
    constexpr static uint32_t kZeroSlotFlag = 0x00800000;
    constexpr static uint32_t kOneSlotFlag = 0x80000000;
    constexpr static uint32_t kTwoSlotFlag = 0x80008000;
    constexpr static uint32_t kThreeSlotFlag = 0x80800000;

    constexpr static uint32_t kOneSlotTagLen = 29;
    constexpr static uint32_t kTwoSlotTagLen = 14;
    constexpr static uint32_t kThreeSlotTagLen = 9;
};

#undef haszero29
#undef hasvalue29
#undef haszero14
#undef hasvalue14
#undef haszero9
#undef hasvalue9

#define haszero45(x) (((x)-0x000000000001ULL) & (~(x)) & 0x100000000000ULL)
#define hasvalue45(x, n) (haszero45((x) ^ (0x00000000001ULL * (n))))
#define haszero22(x) (((x)-0x00000400001ULL) & (~(x)) & 0x80000200000ULL)
#define hasvalue22(x, n) (haszero22((x) ^ (0x00000400001ULL * (n))))
#define haszero15(x) (((x)-0x000040008001ULL) & (~(x)) & 0x100020004000ULL)
#define hasvalue15(x, n) (haszero15((x) ^ (0x000040008001ULL * (n))))

template <>
class SingleTable<12> : public BaseSingleTable<12>
{
  public:
    explicit SingleTable(uint64_t num_buckets)
        : BaseSingleTable(num_buckets)
    {
        static_assert(kBytesPerBucket == 6, "ctor only work on this case");
        char *p{reinterpret_cast<char *>(buckets_.get())};
        for (uint64_t i = 0; i < num_buckets_; ++i)
        {
            *reinterpret_cast<uint32_t *>(p) = (kZeroSlotFlag & 0xffffffff);
            *reinterpret_cast<uint16_t *>(p + 4) = (kZeroSlotFlag >> 32) & 0xffff;
            p += kBytesPerBucket;
        }
    }

    bool FindTagInBucket(const uint64_t i, const uint64_t unmasked_tag) const
    {
        uint64_t bucket;
        std::memcpy(&bucket, buckets_[i].bits_, sizeof(uint64_t));
        const uint64_t flag{bucket & kFlagBitsMask};

        switch (flag)
        {
        case kZeroSlotFlag: {
            return false;
        }
        case kOneSlotFlag: {
            const uint64_t tags{_pext_u64(bucket, kTagBitsMask)};
            return hasvalue45(tags, MaskedTag<kOneSlotTagLen>(unmasked_tag)) ||
                   hasvalue45(tags, MaskedTag<kTwoSlotTagLen>(unmasked_tag)) ||
                   hasvalue45(tags, MaskedTag<kThreeSlotTagLen>(unmasked_tag)) ||
                   hasvalue45(tags, MaskedTag<kFourSlotTagLen>(unmasked_tag));
        }
        case kTwoSlotFlag: {
            const uint64_t tags{_pext_u64(bucket, kTagBitsMask)};
            return hasvalue22(tags, MaskedTag<kTwoSlotTagLen>(unmasked_tag)) ||
                   hasvalue22(tags, MaskedTag<kThreeSlotTagLen>(unmasked_tag)) ||
                   hasvalue22(tags, MaskedTag<kFourSlotTagLen>(unmasked_tag));
        }
        case kThreeSlotFlag: {
            const uint64_t tags{_pext_u64(bucket, kTagBitsMask)};
            return hasvalue15(tags, MaskedTag<kThreeSlotTagLen>(unmasked_tag)) ||
                   hasvalue15(tags, MaskedTag<kFourSlotTagLen>(unmasked_tag));
        }
        default: {
            return hasvalue12(bucket, MaskedTag<kFourSlotTagLen>(unmasked_tag));
        }
        }
    }

    inline bool InsertTagToBucket(const uint64_t i, const uint64_t tag,
                                  const bool kickout, uint64_t &oldtag)
    {
        uint64_t bucket;
        std::memcpy(&bucket, buckets_[i].bits_, sizeof(uint64_t));
        const uint64_t flag{bucket & kFlagBitsMask};

        switch (flag)
        {
        case kZeroSlotFlag: {
            uint64_t tags_and_next_bucket{
                _pext_u64(bucket, 0xffff000000000000)}; // get 16bit next bucket
            tags_and_next_bucket = (tags_and_next_bucket << kOneSlotTagLen) |
                                   MaskedTag<kOneSlotTagLen>(tag);
            tags_and_next_bucket =
                _pdep_u64(tags_and_next_bucket, 0xffff000000000000 | kTagBitsMask) |
                kOneSlotFlag;
            std::memcpy(buckets_[i].bits_, &tags_and_next_bucket, sizeof(uint64_t));
            return true;
        }
        case kOneSlotFlag: {
            uint64_t tags_and_next_bucket{_pext_u64(
                bucket,
                0xffff0000003fffff)}; // get one 22bit tag + 16bit next bucket
            tags_and_next_bucket = (tags_and_next_bucket << kTwoSlotTagLen) |
                                   MaskedTag<kTwoSlotTagLen>(tag);
            tags_and_next_bucket =
                _pdep_u64(tags_and_next_bucket, 0xffff3ff7ff7fffff) | kTwoSlotFlag;
            std::memcpy(buckets_[i].bits_, &tags_and_next_bucket, sizeof(uint64_t));
            return true;
        }
        case kTwoSlotFlag: {
            uint64_t tags_and_next_bucket{_pext_u64(
                bucket,
                0xffff0077ff407fff)}; // get two 15bit tags + 16bit next bucket
            tags_and_next_bucket = (tags_and_next_bucket << kThreeSlotTagLen) |
                                   MaskedTag<kThreeSlotTagLen>(tag);
            tags_and_next_bucket =
                _pdep_u64(tags_and_next_bucket, 0xffff000000000000 | kTagBitsMask) |
                kThreeSlotFlag;
            std::memcpy(buckets_[i].bits_, &tags_and_next_bucket, sizeof(uint64_t));
            return true;
        }
        case kThreeSlotFlag: {
            uint64_t tags_and_next_bucket{_pext_u64(
                bucket,
                0xffff0ff78f7f8fff)}; // get three 12bit tags + 16bit next bucket
            tags_and_next_bucket = (tags_and_next_bucket << kFourSlotTagLen) |
                                   MaskedTag<kFourSlotTagLen>(tag);
            uint64_t tag3{BucketTag<3, kFourSlotTagLen>(tags_and_next_bucket)},
                tag2{BucketTag<2, kFourSlotTagLen>(tags_and_next_bucket)},
                tag1{BucketTag<1, kFourSlotTagLen>(tags_and_next_bucket)};
            if (tag3 > tag2 || tag2 > tag1)
            {
                if (tag3 > tag2)
                {
                    std::swap(tag3, tag2);
                }
                if (tag3 > tag1)
                {
                    std::swap(tag3, tag1);
                }
                if (tag2 > tag1)
                {
                    std::swap(tag2, tag1);
                }
                tags_and_next_bucket &= 0xffff000000000fff;
                tags_and_next_bucket |= (tag3 << 3 * kFourSlotTagLen) |
                                        (tag2 << 2 * kFourSlotTagLen) |
                                        (tag1 << kFourSlotTagLen);
            }
            std::memcpy(buckets_[i].bits_, &tags_and_next_bucket, sizeof(uint64_t));
            assert((BucketTag<3, kFourSlotTagLen>(tags_and_next_bucket) <=
                        BucketTag<2, kFourSlotTagLen>(tags_and_next_bucket) &&
                    BucketTag<2, kFourSlotTagLen>(tags_and_next_bucket) <=
                        BucketTag<1, kFourSlotTagLen>(tags_and_next_bucket)));
            return true;
        }
        default: {
            if (kickout)
            {
                uint64_t r{rand() % kTagsPerBucket};
                oldtag = BucketTag<kFourSlotTagLen>(bucket, r);
                char *p = reinterpret_cast<char *>(&bucket) + (r + (r >> 1));

                uint16_t updated;
                std::memcpy(&updated, p, sizeof(updated));
                if ((r & 1) == 0)
                {
                    updated &= 0xf000;
                    updated |= MaskedTag<kFourSlotTagLen>(tag);
                }
                else
                {
                    updated &= 0x000f;
                    updated |= MaskedTag<kFourSlotTagLen>(tag) << 4;
                }
                std::memcpy(p, &updated, sizeof(updated));

                uint64_t tag3{BucketTag<3, kFourSlotTagLen>(bucket)},
                    tag2{BucketTag<2, kFourSlotTagLen>(bucket)},
                    tag1{BucketTag<1, kFourSlotTagLen>(bucket)};
                if (tag3 > tag2 || tag2 > tag1)
                {
                    if (tag3 > tag2)
                    {
                        std::swap(tag3, tag2);
                    }
                    if (tag3 > tag1)
                    {
                        std::swap(tag3, tag1);
                    }
                    if (tag2 > tag1)
                    {
                        std::swap(tag2, tag1);
                    }
                    bucket &= 0xffff000000000fff;
                    bucket |= (tag3 << 3 * kFourSlotTagLen) |
                              (tag2 << 2 * kFourSlotTagLen) | (tag1 << kFourSlotTagLen);
                }
                std::memcpy(buckets_[i].bits_, &bucket, sizeof(uint64_t));
            }
            return false;
        }
        }
    }

    void FindMaxMatchingTag(const uint64_t i, const uint64_t unmasked_tag,
                            uint64_t *max_bucket_idx,
                            uint32_t *max_tag_length) const
    {
        uint64_t bucket;
        std::memcpy(&bucket, buckets_[i].bits_, sizeof(uint64_t));
        const uint64_t flag{bucket & kFlagBitsMask};

        switch (flag)
        {
        case kZeroSlotFlag: {
            return;
        }
        case kOneSlotFlag: {
            const uint64_t tags{_pext_u64(bucket, kTagBitsMask)};
#define HELPER(tag_length)                                     \
    if (*max_tag_length >= tag_length)                         \
    {                                                          \
        return; /* Longer or equal matched tag found */        \
    }                                                          \
    if (hasvalue45(tags, MaskedTag<tag_length>(unmasked_tag))) \
    {                                                          \
        *max_bucket_idx = i;                                   \
        *max_tag_length = tag_length;                          \
        return;                                                \
    }

            HELPER(kOneSlotTagLen)
            HELPER(kTwoSlotTagLen)
            HELPER(kThreeSlotTagLen)
            HELPER(kFourSlotTagLen)
#undef HELPER
            return;
        }
        case kTwoSlotFlag: {
            const uint64_t tags{_pext_u64(bucket, kTagBitsMask)};
#define HELPER(tag_length)                                     \
    if (*max_tag_length >= tag_length)                         \
    {                                                          \
        return; /* Longer or equal matched tag found */        \
    }                                                          \
    if (hasvalue22(tags, MaskedTag<tag_length>(unmasked_tag))) \
    {                                                          \
        *max_bucket_idx = i;                                   \
        *max_tag_length = tag_length;                          \
        return;                                                \
    }

            HELPER(kTwoSlotTagLen)
            HELPER(kThreeSlotTagLen)
            HELPER(kFourSlotTagLen)
#undef HELPER
            return;
        }
        case kThreeSlotFlag: {
            const uint64_t tags{_pext_u64(bucket, kTagBitsMask)};
#define HELPER(tag_length)                                     \
    if (*max_tag_length >= tag_length)                         \
    {                                                          \
        return; /* Longer or equal matched tag found */        \
    }                                                          \
    if (hasvalue15(tags, MaskedTag<tag_length>(unmasked_tag))) \
    {                                                          \
        *max_bucket_idx = i;                                   \
        *max_tag_length = tag_length;                          \
        return;                                                \
    }

            HELPER(kThreeSlotTagLen)
            HELPER(kFourSlotTagLen)
#undef HELPER
            return;
        }
        default: {
            // encoding requirements
            assert((BucketTag<3, kFourSlotTagLen>(bucket) <=
                        BucketTag<2, kFourSlotTagLen>(bucket) &&
                    BucketTag<2, kFourSlotTagLen>(bucket) <=
                        BucketTag<1, kFourSlotTagLen>(bucket)));

            if (*max_tag_length >= kFourSlotTagLen)
            {
                return; // Longer or equal matched tag found
            }
            if (hasvalue12(bucket, MaskedTag<kFourSlotTagLen>(unmasked_tag)))
            {
                *max_bucket_idx = i;
                *max_tag_length = kFourSlotTagLen;
            }
            return;
        }
        }
    }

    void DeleteTagFromBucket(uint64_t bucket_idx, const uint64_t masked_tag)
    {
        uint64_t bucket;
        std::memcpy(&bucket, buckets_[bucket_idx].bits_, sizeof(uint64_t));
        const uint64_t flag{bucket & kFlagBitsMask};

        switch (flag)
        {
        case kZeroSlotFlag: {
            // never delete on empty bucket
            assert(false);
            __builtin_unreachable();
            return;
        }
        case kOneSlotFlag: {
            bucket &= 0xffff000000000000;
            bucket |= kZeroSlotFlag;
            std::memcpy(buckets_[bucket_idx].bits_, &bucket, sizeof(uint64_t));
            return;
        }
        case kTwoSlotFlag: {
            uint64_t masks[]{0xffff3ff7ff400000, 0xffff0000003fffff};
            const uint64_t tags{_pext_u64(bucket, kTagBitsMask)};
            for (uint32_t slot_idx = 0; slot_idx < 2; ++slot_idx)
            {
                if (BucketTag<kTwoSlotTagLen>(tags, slot_idx) == masked_tag)
                {
                    bucket = _pdep_u64(_pext_u64(bucket, masks[slot_idx]),
                                       0xffff0000003fffff) |
                             kOneSlotFlag;
                    std::memcpy(buckets_[bucket_idx].bits_, &bucket, sizeof(uint64_t));
                    return;
                }
            }
            assert(false);
            __builtin_unreachable();
        }
        case kThreeSlotFlag: {
            uint64_t masks[]{0xffff7ff7ff7f8000, 0xffff7ff780007fff,
                             0xffff00007f7fffff};
            const uint64_t tags{_pext_u64(bucket, kTagBitsMask)};
            for (uint32_t slot_idx = 0; slot_idx < 3; ++slot_idx)
            {
                if (BucketTag<kThreeSlotTagLen>(tags, slot_idx) == masked_tag)
                {
                    bucket = _pdep_u64(_pext_u64(bucket, masks[slot_idx]),
                                       0xffff0077ff407fff) |
                             kTwoSlotFlag;
                    std::memcpy(buckets_[bucket_idx].bits_, &bucket, sizeof(uint64_t));
                    return;
                }
            }
            assert(false);
            __builtin_unreachable();
        }
        default: {
            uint64_t masks[]{0xfffffffffffff000, 0xffffffffff000fff,
                             0xfffffff000ffffff, 0xffff000fffffffff};
            for (uint32_t slot_idx = 0; slot_idx < 4; ++slot_idx)
            {
                if (BucketTag<kFourSlotTagLen>(bucket, slot_idx) == masked_tag)
                {
                    bucket = _pdep_u64(_pext_u64(bucket, masks[slot_idx]),
                                       0xffff0ff78f7f8fff) |
                             kThreeSlotFlag;
                    std::memcpy(buckets_[bucket_idx].bits_, &bucket, sizeof(uint64_t));
                    return;
                }
            }
            assert(false);
            __builtin_unreachable();
        }
        }
    }

    void BucketCountStat(uint64_t *counter)
    {
        for (uint64_t i = 0; i < num_buckets_; ++i)
        {
            uint64_t bucket;
            std::memcpy(&bucket, buckets_[i].bits_, sizeof(uint64_t));
            const uint64_t flag{bucket & kFlagBitsMask};
            switch (flag)
            {
            case kZeroSlotFlag: {
                ++counter[0];
                break;
            }
            case kOneSlotFlag: {
                ++counter[1];
                break;
            }
            case kTwoSlotFlag: {
                ++counter[2];
                break;
            }
            case kThreeSlotFlag: {
                ++counter[3];
                break;
            }
            default: {
                ++counter[4];
                break;
            }
            }
        }
    }

  private:
    constexpr static uint64_t kFlagBitsMask = 0x0000800800800000;
    constexpr static uint64_t kTagBitsMask = 0x00007ff7ff7fffff;

    // other case is four slots
    constexpr static uint64_t kZeroSlotFlag = 0x0000000800000000;
    constexpr static uint64_t kOneSlotFlag = 0x0000800000000000;
    constexpr static uint64_t kTwoSlotFlag = 0x0000800000800000;
    constexpr static uint64_t kThreeSlotFlag = 0x0000800800000000;

    constexpr static uint32_t kOneSlotTagLen = 45;
    constexpr static uint32_t kTwoSlotTagLen = 22;
    constexpr static uint32_t kThreeSlotTagLen = 15;
};

#undef haszero45
#undef hasvalue45
#undef haszero22
#undef hasvalue22
#undef haszero15
#undef hasvalue15

#define haszero61(x) \
    (((x)-0x0000000000000001ULL) & (~(x)) & 0x1000000000000000ULL)
#define hasvalue61(x, n) (haszero61((x) ^ (0x0000000000000001ULL * (n))))
#define haszero30(x) \
    (((x)-0x0000000040000001ULL) & (~(x)) & 0x0800000020000000ULL)
#define hasvalue30(x, n) (haszero30((x) ^ (0x0000000040000001ULL * (n))))
#define haszero20(x) \
    (((x)-0x0000010000100001ULL) & (~(x)) & 0x0800008000080000ULL)
#define hasvalue20(x, n) (haszero20((x) ^ (0x0000010000100001ULL * (n))))

template <>
class SingleTable<16> : public BaseSingleTable<16>
{
  public:
    explicit SingleTable(uint64_t num_buckets)
        : BaseSingleTable(num_buckets)
    {
        static_assert(kBytesPerBucket == 8, "ctor only work on this case");
        char *p{reinterpret_cast<char *>(buckets_.get())};
        uint64_t src{kZeroSlotFlag};
        for (uint64_t i = 0; i < num_buckets_; ++i)
        {
            std::memcpy(p, &src, sizeof(uint64_t));
            p += kBytesPerBucket;
        }
    }

    bool FindTagInBucket(const uint64_t i, const uint64_t unmasked_tag) const
    {
        uint64_t bucket;
        std::memcpy(&bucket, buckets_[i].bits_, sizeof(uint64_t));
        const uint64_t flag{bucket & kFlagBitsMask};

        switch (flag)
        {
        case kZeroSlotFlag: {
            return false;
        }
        case kOneSlotFlag: {
            const uint64_t tags{_pext_u64(bucket, kTagBitsMask)};
            return hasvalue61(tags, MaskedTag<kOneSlotTagLen>(unmasked_tag)) ||
                   hasvalue61(tags, MaskedTag<kTwoSlotTagLen>(unmasked_tag)) ||
                   hasvalue61(tags, MaskedTag<kThreeSlotTagLen>(unmasked_tag)) ||
                   hasvalue61(tags, MaskedTag<kFourSlotTagLen>(unmasked_tag));
        }
        case kTwoSlotFlag: {
            const uint64_t tags{_pext_u64(bucket, kTagBitsMask)};
            return hasvalue30(tags, MaskedTag<kTwoSlotTagLen>(unmasked_tag)) ||
                   hasvalue30(tags, MaskedTag<kThreeSlotTagLen>(unmasked_tag)) ||
                   hasvalue30(tags, MaskedTag<kFourSlotTagLen>(unmasked_tag));
        }
        case kThreeSlotFlag: {
            const uint64_t tags{_pext_u64(bucket, kTagBitsMask)};
            return hasvalue20(tags, MaskedTag<kThreeSlotTagLen>(unmasked_tag)) ||
                   hasvalue20(tags, MaskedTag<kFourSlotTagLen>(unmasked_tag));
        }
        default: {
            // four slot

            // encoding requirements
            assert((BucketTag<3, kFourSlotTagLen>(bucket) <=
                        BucketTag<2, kFourSlotTagLen>(bucket) &&
                    BucketTag<2, kFourSlotTagLen>(bucket) <=
                        BucketTag<1, kFourSlotTagLen>(bucket)));
            return hasvalue16(bucket, MaskedTag<kFourSlotTagLen>(unmasked_tag));
        }
        }
    }

    inline bool InsertTagToBucket(const uint64_t i, const uint64_t tag,
                                  const bool kickout, uint64_t &oldtag)
    {
        uint64_t bucket;
        std::memcpy(&bucket, buckets_[i].bits_, sizeof(uint64_t));
        const uint64_t flag{bucket & kFlagBitsMask};

        switch (flag)
        {
        case kZeroSlotFlag: {
            bucket = _pdep_u64(tag, kTagBitsMask) | kOneSlotFlag;
            std::memcpy(buckets_[i].bits_, &bucket, sizeof(uint64_t));
            return true;
        }
        case kOneSlotFlag: {
            uint64_t tags{
                _pext_u64(bucket, 0x000000003fffffff)}; // get one 30bit tags
            tags |= MaskedTag<kTwoSlotTagLen>(tag) << kTwoSlotTagLen;
            tags = _pdep_u64(tags, kTagBitsMask) | kTwoSlotFlag;
            std::memcpy(buckets_[i].bits_, &tags, sizeof(uint64_t));
            return true;
        }
        case kTwoSlotFlag: {
            uint64_t tags{
                _pext_u64(bucket, 0x000f7fff400fffff)}; // get two 20bit tags
            tags |= MaskedTag<kThreeSlotTagLen>(tag) << (2 * kThreeSlotTagLen);
            tags = _pdep_u64(tags, kTagBitsMask) | kThreeSlotFlag;
            std::memcpy(buckets_[i].bits_, &tags, sizeof(uint64_t));
            return true;
        }
        case kThreeSlotFlag: {
            uint64_t tags{
                _pext_u64(bucket, 0x03ff7e1f7ff0ffff)}; // get three 16bit tags
            tags |= MaskedTag<kFourSlotTagLen>(tag) << (3 * kFourSlotTagLen);
            uint64_t tag3{BucketTag<3, kFourSlotTagLen>(tags)},
                tag2{BucketTag<2, kFourSlotTagLen>(tags)},
                tag1{BucketTag<1, kFourSlotTagLen>(tags)};
            if (tag3 > tag2 || tag2 > tag1)
            {
                if (tag3 > tag2)
                {
                    std::swap(tag3, tag2);
                }
                if (tag3 > tag1)
                {
                    std::swap(tag3, tag1);
                }
                if (tag2 > tag1)
                {
                    std::swap(tag2, tag1);
                }
                tags &= 0x000000000000ffff;
                tags |= (tag3 << 3 * kFourSlotTagLen) |
                        (tag2 << 2 * kFourSlotTagLen) | (tag1 << kFourSlotTagLen);
            }
            std::memcpy(buckets_[i].bits_, &tags, sizeof(uint64_t));
            assert((BucketTag<3, kFourSlotTagLen>(tags) <=
                        BucketTag<2, kFourSlotTagLen>(tags) &&
                    BucketTag<2, kFourSlotTagLen>(tags) <=
                        BucketTag<1, kFourSlotTagLen>(tags)));
            return true;
        }
        default: {
            if (kickout)
            {
                static_assert(kFourSlotTagLen == 16,
                              "can not access bucket as uint16_t array");
                uint64_t r{rand() % kTagsPerBucket};
                oldtag = BucketTag<kFourSlotTagLen>(bucket, r);
                char *p = reinterpret_cast<char *>(&bucket) + r * 2;

                uint16_t updated = MaskedTag<kFourSlotTagLen>(tag);
                std::memcpy(p, &updated, sizeof(updated));

                uint64_t tag3{BucketTag<3, kFourSlotTagLen>(bucket)},
                    tag2{BucketTag<2, kFourSlotTagLen>(bucket)},
                    tag1{BucketTag<1, kFourSlotTagLen>(bucket)};
                if (tag3 > tag2 || tag2 > tag1)
                {
                    if (tag3 > tag2)
                    {
                        std::swap(tag3, tag2);
                    }
                    if (tag3 > tag1)
                    {
                        std::swap(tag3, tag1);
                    }
                    if (tag2 > tag1)
                    {
                        std::swap(tag2, tag1);
                    }
                    bucket &= 0x000000000000ffff;
                    bucket |= (tag3 << 3 * kFourSlotTagLen) |
                              (tag2 << 2 * kFourSlotTagLen) | (tag1 << kFourSlotTagLen);
                }
                std::memcpy(buckets_[i].bits_, &bucket, sizeof(uint64_t));
                assert((BucketTag<3, kFourSlotTagLen>(bucket) <=
                            BucketTag<2, kFourSlotTagLen>(bucket) &&
                        BucketTag<2, kFourSlotTagLen>(bucket) <=
                            BucketTag<1, kFourSlotTagLen>(bucket)));
            }
            return false;
        }
        }
    }

    void FindMaxMatchingTag(const uint64_t i, const uint64_t unmasked_tag,
                            uint64_t *max_bucket_idx,
                            uint32_t *max_tag_length) const
    {
        uint64_t bucket;
        std::memcpy(&bucket, buckets_[i].bits_, sizeof(uint64_t));
        const uint64_t flag{bucket & kFlagBitsMask};

        switch (flag)
        {
        case kZeroSlotFlag: {
            return;
        }
        case kOneSlotFlag: {
            const uint64_t tags{_pext_u64(bucket, kTagBitsMask)};
#define HELPER(tag_length)                                     \
    if (*max_tag_length >= tag_length)                         \
    {                                                          \
        return; /* Longer or equal matched tag found */        \
    }                                                          \
    if (hasvalue61(tags, MaskedTag<tag_length>(unmasked_tag))) \
    {                                                          \
        *max_bucket_idx = i;                                   \
        *max_tag_length = tag_length;                          \
        return;                                                \
    }

            HELPER(kOneSlotTagLen)
            HELPER(kTwoSlotTagLen)
            HELPER(kThreeSlotTagLen)
            HELPER(kFourSlotTagLen)
#undef HELPER
            return;
        }
        case kTwoSlotFlag: {
            const uint64_t tags{_pext_u64(bucket, kTagBitsMask)};
#define HELPER(tag_length)                                     \
    if (*max_tag_length >= tag_length)                         \
    {                                                          \
        return; /* Longer or equal matched tag found */        \
    }                                                          \
    if (hasvalue30(tags, MaskedTag<tag_length>(unmasked_tag))) \
    {                                                          \
        *max_bucket_idx = i;                                   \
        *max_tag_length = tag_length;                          \
        return;                                                \
    }

            HELPER(kTwoSlotTagLen)
            HELPER(kThreeSlotTagLen)
            HELPER(kFourSlotTagLen)
#undef HELPER
            return;
        }
        case kThreeSlotFlag: {
            const uint64_t tags{_pext_u64(bucket, kTagBitsMask)};
#define HELPER(tag_length)                                     \
    if (*max_tag_length >= tag_length)                         \
    {                                                          \
        return; /* Longer or equal matched tag found */        \
    }                                                          \
    if (hasvalue20(tags, MaskedTag<tag_length>(unmasked_tag))) \
    {                                                          \
        *max_bucket_idx = i;                                   \
        *max_tag_length = tag_length;                          \
        return;                                                \
    }

            HELPER(kThreeSlotTagLen)
            HELPER(kFourSlotTagLen)
#undef HELPER
            return;
        }
        default: {
            // encoding requirements
            assert((BucketTag<3, kFourSlotTagLen>(bucket) <=
                        BucketTag<2, kFourSlotTagLen>(bucket) &&
                    BucketTag<2, kFourSlotTagLen>(bucket) <=
                        BucketTag<1, kFourSlotTagLen>(bucket)));

            if (*max_tag_length >= kFourSlotTagLen)
            {
                return; // Longer or equal matched tag found
            }
            if (hasvalue16(bucket, MaskedTag<kFourSlotTagLen>(unmasked_tag)))
            {
                *max_bucket_idx = i;
                *max_tag_length = kFourSlotTagLen;
            }
            return;
        }
        }
    }

    void DeleteTagFromBucket(uint64_t bucket_idx, const uint64_t masked_tag)
    {
        uint64_t bucket;
        std::memcpy(&bucket, buckets_[bucket_idx].bits_, sizeof(uint64_t));
        const uint64_t flag{bucket & kFlagBitsMask};

        switch (flag)
        {
        case kZeroSlotFlag: {
            // never delete on empty bucket
            assert(false);
            __builtin_unreachable();
            return;
        }
        case kOneSlotFlag: {
            bucket = kZeroSlotFlag;
            std::memcpy(buckets_[bucket_idx].bits_, &bucket, sizeof(uint64_t));
            return;
        }
        case kTwoSlotFlag: {
            uint64_t masks[]{0x3fff7fff40000000, 0x000000003fffffff};
            const uint64_t tags{_pext_u64(bucket, kTagBitsMask)};
            for (uint32_t slot_idx = 0; slot_idx < 2; ++slot_idx)
            {
                if (BucketTag<kTwoSlotTagLen>(tags, slot_idx) == masked_tag)
                {
                    bucket =
                        _pdep_u64(_pext_u64(bucket, masks[slot_idx]), kTagBitsMask) |
                        kOneSlotFlag;
                    std::memcpy(buckets_[bucket_idx].bits_, &bucket, sizeof(uint64_t));
                    return;
                }
            }
            assert(false);
            __builtin_unreachable();
        }
        case kThreeSlotFlag: {
            uint64_t masks[]{0x3fff7fff7ff00000, 0x3fff7e00000fffff,
                             0x000001ff7fffffff};
            const uint64_t tags{_pext_u64(bucket, kTagBitsMask)};
            for (uint32_t slot_idx = 0; slot_idx < 3; ++slot_idx)
            {
                if (BucketTag<kThreeSlotTagLen>(tags, slot_idx) == masked_tag)
                {
                    bucket = _pdep_u64(_pext_u64(bucket, masks[slot_idx]),
                                       0x000f7fff400fffff) |
                             kTwoSlotFlag;
                    std::memcpy(buckets_[bucket_idx].bits_, &bucket, sizeof(uint64_t));
                    return;
                }
            }
            assert(false);
            __builtin_unreachable();
        }
        default: {
            uint64_t masks[]{0xffffffffffff0000, 0xffffffff0000ffff,
                             0xffff0000ffffffff, 0x0000ffffffffffff};
            for (uint32_t slot_idx = 0; slot_idx < 4; ++slot_idx)
            {
                if (BucketTag<kFourSlotTagLen>(bucket, slot_idx) == masked_tag)
                {
                    bucket = _pdep_u64(_pext_u64(bucket, masks[slot_idx]),
                                       0x03ff7e1f7ff0ffff) |
                             kThreeSlotFlag;
                    std::memcpy(buckets_[bucket_idx].bits_, &bucket, sizeof(uint64_t));
                    return;
                }
            }
            assert(false);
            __builtin_unreachable();
        }
        }
    }

    bool AllZero()
    {
        char *p = reinterpret_cast<char *>(buckets_.get());
        for (uint64_t i = 0; i < num_buckets_; ++i, p += 8)
        {
            if (*reinterpret_cast<uint64_t *>(p) != kZeroSlotFlag)
            {
                return false;
            }
        }
        return true;
    }

  private:
    constexpr static uint64_t kFlagBitsMask = 0x8000800080000000;
    constexpr static uint64_t kTagBitsMask = ~kFlagBitsMask;

    // other case is four slots
    constexpr static uint64_t kZeroSlotFlag = 0x0000800000000000;
    constexpr static uint64_t kOneSlotFlag = 0x8000000000000000;
    constexpr static uint64_t kTwoSlotFlag = 0x8000000080000000;
    constexpr static uint64_t kThreeSlotFlag = 0x8000800000000000;

    constexpr static uint32_t kOneSlotTagLen = 61;
    constexpr static uint32_t kTwoSlotTagLen = 30;
    constexpr static uint32_t kThreeSlotTagLen = 20;
};

#undef haszero61
#undef hasvalue61
#undef haszero30
#undef hasvalue30
#undef haszero20
#undef hasvalue20

}

#endif
