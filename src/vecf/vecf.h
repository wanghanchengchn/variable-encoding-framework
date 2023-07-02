#ifndef VECF_H_
#define VECF_H_

#include "hashutil.h"
#include "vecf/singletable.h"

namespace vecf
{
// maximum number of cuckoo kicks before claiming failure
const size_t kMaxCuckooCount = 500;

template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType = SingleTable,
          typename HashFamily = hashutil::TwoIndependentMultiplyShift>
class VECF
{
  private:
    TableType<bits_per_item> *table_;

    size_t num_items_;

    typedef struct
    {
        uint64_t index;
        uint64_t tag;
        bool used;
    } VictimCache;

    VictimCache victim_;

    HashFamily hasher_one_, hasher_two_;

    inline uint64_t IndexHash(uint64_t hv) const
    {
        return hv & (table_->NumBuckets() - 1);
    }

    inline uint64_t TagHash(uint64_t hv) const
    {
        return hv;
    }

    inline uint64_t AltIndex(const uint64_t index,
                             const uint64_t unmasked_tag) const
    {
        return IndexHash(index ^
                         (MaskedTag<bits_per_item>(unmasked_tag) * 0x5bd1e995));
    }

    inline void GenerateIndexTagHash(const ItemType &item, uint64_t *index,
                                     uint64_t *unmasked_tag) const
    {
        *index = IndexHash(hasher_one_(item));
        *unmasked_tag = TagHash(hasher_two_(item));
    }

    bool InsertImpl(const uint64_t i, const uint64_t unmasked_tag);

  public:
    explicit VECF(const size_t max_num_keys)
        : num_items_(0), victim_(), hasher_one_(), hasher_two_()
    {
        size_t assoc = 4;
        size_t num_buckets =
            upperpower2(std::max<uint64_t>(1, max_num_keys / assoc));
        double frac = (double)max_num_keys / num_buckets / assoc;
        if (frac > 0.96)
        {
            num_buckets <<= 1;
        }
        victim_.used = false;
        table_ = new TableType<bits_per_item>(num_buckets);
    }

    ~VECF()
    {
        delete table_;
    }

    bool Insert(const ItemType &item);

    bool Lookup(const ItemType &item) const;

    bool Delete(const ItemType &item);

    size_t GetItemNum() const
    {
        return num_items_;
    }

    size_t SizeInBytes() const
    {
        return table_->SizeInBytes();
    }

    double LoadFactor() const
    {
        return 1.0 * GetItemNum() / table_->SizeInTags();
    }

    double BitsPerItem() const
    {
        return 8.0 * table_->SizeInBytes() / GetItemNum();
    }

    void BucketCountStat(uint64_t *counter)
    {
        table_->BucketCountStat(counter);
    }
};

template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType, typename HashFamily>
bool VECF<ItemType, bits_per_item, TableType, HashFamily>::InsertImpl(
    const uint64_t i, const uint64_t unmasked_tag)
{
    uint64_t curindex{i};
    uint64_t curtag{unmasked_tag}, oldtag;

    for (uint32_t count = 0; count < kMaxCuckooCount; count++)
    {
        bool kickout{count > 0}; // if count == 0 and insert failed, try another
                                 // bucket instead of kickout
        oldtag = 0;
        if (table_->InsertTagToBucket(curindex, curtag, kickout, oldtag))
        {
            ++num_items_;
            return true;
        }
        if (kickout)
        {
            curtag = oldtag;
        }

        curindex = AltIndex(curindex, curtag);
    }

    victim_.index = curindex;
    victim_.tag = curtag;
    victim_.used = true;
    return true;
}

template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType, typename HashFamily>
bool VECF<ItemType, bits_per_item, TableType, HashFamily>::Insert(
    const ItemType &item)
{
    size_t i;
    uint64_t tag;

    if (victim_.used)
    {
        return false;
    }

    GenerateIndexTagHash(item, &i, &tag);

    return InsertImpl(i, tag);
}

template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType, typename HashFamily>
bool VECF<ItemType, bits_per_item, TableType, HashFamily>::Lookup(
    const ItemType &item) const
{
    uint64_t i1, i2, unmasked_tag;

    GenerateIndexTagHash(item, &i1, &unmasked_tag);
    i2 = AltIndex(i1, unmasked_tag);

    if (victim_.used && (victim_.index == i1 || victim_.index == i2) &&
        victim_.tag == MaskedTag<bits_per_item>(unmasked_tag))
    {
        return true;
    }

    return table_->FindTagInBucket(i1, unmasked_tag) ||
           table_->FindTagInBucket(i2, unmasked_tag);
}

template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType, typename HashFamily>
bool VECF<ItemType, bits_per_item, TableType, HashFamily>::Delete(
    const ItemType &item)
{
    uint64_t i1, i2, unmasked_tag;
    GenerateIndexTagHash(item, &i1, &unmasked_tag);
    i2 = AltIndex(i1, unmasked_tag);

    uint64_t max_bucket_idx{};
    uint32_t max_tag_length{};
    table_->FindMaxMatchingTag(i1, unmasked_tag, &max_bucket_idx,
                               &max_tag_length);
    table_->FindMaxMatchingTag(i2, unmasked_tag, &max_bucket_idx,
                               &max_tag_length);

    if (max_tag_length == 0)
    {
        if (victim_.used && (victim_.index == i1 || victim_.index == i2) &&
            victim_.tag == MaskedTag<bits_per_item>(unmasked_tag))
        {
            victim_.used = false;
            --num_items_;
            return true;
        }
        return false;
    }

    table_->DeleteTagFromBucket(max_bucket_idx,
                                MaskedTag(unmasked_tag, max_tag_length));

    if (victim_.used)
    {
        victim_.used = false;
        uint64_t i{victim_.index}, tag{victim_.tag};
        InsertImpl(i, tag);
    }
    --num_items_;
    return true;
}

}

#endif