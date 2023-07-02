#ifndef VECBF_H_
#define VECBF_H_

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>

#include "hashutil.h"

namespace vecbf
{

template <typename ItemType, uint64_t kBitsPerCounter,
          typename HashFunction = hashutil::TwoIndependentMultiplyShift>
class VECBF
{
  private:
    static constexpr inline uint64_t LowMask(uint64_t n)
    {
        return (1ULL << n) - 1;
    }

    constexpr static uint64_t kPhase1LowerCounterMaxValue{
        LowMask(kBitsPerCounter / 2)};
    constexpr static uint64_t kPhase1UpperCounterBase{1 << (kBitsPerCounter / 2)};
    constexpr static uint64_t kCounterMask{LowMask(kBitsPerCounter)};

    bool is_overflow{false};
    uint64_t num_items_{0};

    const uint64_t max_num_keys_, counter_num_, hash_function_num_, table_size_;
    HashFunction hasher_;
    std::unique_ptr<uint64_t[]> table_;

    static uint64_t OptimalBitNum(uint64_t max_num_keys, double false_positive)
    {
        return (uint64_t)(max_num_keys * (-1.0 * log(false_positive)) /
                          (log(2) * log(2)));
    }

    static uint64_t OptimalHashFunctionNum(uint64_t max_num_keys,
                                           uint64_t counter_num)
    {
        uint64_t hash_function_num{
            (uint64_t)round(counter_num * log(2) / max_num_keys)};
        return hash_function_num > 1 ? hash_function_num : 1;
    }

    uint64_t GetCounter(uint64_t idx) const
    {
        uint64_t bit_idx{idx * kBitsPerCounter};
        uint64_t table_idx{bit_idx / 64}, slot_idx{bit_idx % 64};
        int64_t spillbits{static_cast<int64_t>(slot_idx + kBitsPerCounter) - 64};
        uint64_t ret{(table_[table_idx] >> slot_idx) & kCounterMask};
        if (spillbits > 0)
        {
            ++table_idx;
            uint64_t x{table_[table_idx] & LowMask(spillbits)};
            ret |= x << (kBitsPerCounter - spillbits);
        }
        return ret;
    }

    void SetCounter(uint64_t idx, uint64_t val)
    {
        uint64_t bit_idx{idx * kBitsPerCounter};
        uint64_t table_idx{bit_idx / 64}, slot_idx{bit_idx % 64};
        int64_t spillbits{static_cast<int64_t>(slot_idx + kBitsPerCounter) - 64};
        val &= kCounterMask;
        table_[table_idx] &= ~(kCounterMask << slot_idx);
        table_[table_idx] |= val << slot_idx;
        if (spillbits > 0)
        {
            ++table_idx;
            table_[table_idx] &= ~LowMask(spillbits);
            table_[table_idx] |= val >> (kBitsPerCounter - spillbits);
        }
    }

    static inline uint64_t Phase1LowerCounter(uint64_t counter)
    {
        return counter & LowMask(kBitsPerCounter / 2);
    }

    static inline uint64_t Phase1UpperCounter(uint64_t counter)
    {
        return (counter >> (kBitsPerCounter / 2));
    }

    void SwitchToPhase2()
    {
        for (uint64_t i = 0; i < counter_num_; ++i)
        {
            uint64_t counter{GetCounter(i)};
            SetCounter(i, counter & LowMask(kBitsPerCounter / 2));
        }
        is_overflow = true;
    }

  public:
    VECBF(const uint64_t max_num_keys, double false_positive = 0.04)
        : max_num_keys_(max_num_keys),
          counter_num_(OptimalBitNum(max_num_keys, false_positive)),
          hash_function_num_(OptimalHashFunctionNum(max_num_keys, counter_num_)),
          table_size_((counter_num_ * kBitsPerCounter + 63) / 64),
          hasher_(),
          table_(new uint64_t[table_size_])
    {
        memset(table_.get(), 0, table_size_ * sizeof(uint64_t));
    }

    bool Insert(const ItemType &item)
    {
        const uint64_t hash{hasher_(item)};
        const uint64_t hash1{hash & LowMask(32)}, hash2{hash >> 32};
        const uint64_t hash_function_num{
            is_overflow == true ? hash_function_num_ : hash_function_num_ * 2};

        if (is_overflow == false)
        {
            for (uint64_t i = 0; i < hash_function_num; ++i)
            {
                uint64_t combine_hash{hash1 + hash2 * i};
                uint64_t idx{combine_hash % counter_num_};
                uint64_t counter{GetCounter(idx)};

                // if (Phase1LowerCounter(counter) == LowMask(kBitsPerCounter / 2) ||
                // Phase1UpperCounter(counter) == LowMask(kBitsPerCounter -
                // kBitsPerCounter / 2))
                // {
                //     cout << "counter overflow" << endl;
                //     return false;
                // }

                counter += (i >= hash_function_num_) ? kPhase1UpperCounterBase : 1;

                SetCounter(idx, counter);
            }

            if (num_items_ >= int(max_num_keys_ * 0.5))
            {
                SwitchToPhase2();
            }
        }
        else
        {
            for (uint64_t i = 0; i < hash_function_num; ++i)
            {
                uint64_t combine_hash{hash1 + hash2 * i};
                uint64_t idx{combine_hash % counter_num_};
                uint64_t counter{GetCounter(idx)};

                // if (counter == LowMask(kBitsPerCounter))
                // {
                //     cout << "counter overflow" << endl;
                //     return false;
                // }

                SetCounter(idx, counter + 1);
            }
        }

        num_items_++;

        return true;
    }

    bool Lookup(const ItemType &key) const
    {
        const uint64_t hash{hasher_(key)};
        const uint64_t hash1{hash & LowMask(32)}, hash2{hash >> 32};
        const uint64_t hash_function_num{
            is_overflow == true ? hash_function_num_ : hash_function_num_ * 2};

        for (uint64_t i = 0; i < hash_function_num; ++i)
        {
            uint64_t combine_hash{hash1 + hash2 * i};
            uint64_t idx{combine_hash % counter_num_};
            if (GetCounter(idx) == 0)
            {
                return false;
            }
        }
        return true;
    }

    bool Delete(const ItemType &key)
    {
        const uint64_t hash{hasher_(key)};
        const uint64_t hash1{hash & LowMask(32)}, hash2{hash >> 32};
        const uint64_t hash_function_num{
            is_overflow == true ? hash_function_num_ : hash_function_num_ * 2};

        for (uint64_t i = 0; i < hash_function_num; ++i)
        {
            uint64_t combine_hash{hash1 + hash2 * i};
            uint64_t idx{combine_hash % counter_num_};
            uint64_t counter{GetCounter(idx)};

            if (counter == 0)
            {
                // cout << "counter is zero" << endl;
                return false;
            }

            counter -= (i >= hash_function_num_) ? kPhase1UpperCounterBase : 1;

            SetCounter(idx, counter);
        }

        num_items_--;
        return true;
    }

    size_t Size() const
    {
        return num_items_;
    }
    size_t SizeInBytes() const
    {
        return table_size_ * sizeof(uint64_t);
    }
    double LoadFactor() const
    {
        return 1.0 * Size() / max_num_keys_;
    }
    double BitsPerItem() const
    {
        return 8.0 * SizeInBytes() / Size();
    }

    bool CheckAllZero()
    {
        for (uint64_t i = 0; i < table_size_; i++)
        {
            if (table_[i] != 0)
            {
                // cout << "error: no zero" << endl;
                return false;
            }
        }

        return true;
    }
};

}

#endif