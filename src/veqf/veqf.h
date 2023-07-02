#ifndef VEQF_H_
#define VEQF_H_

#include <cassert>
#include <cstdint>
#include <cstring>
#include <memory>

#include "hashutil.h"
#include "veqf/bitsutil.h"

namespace veqf
{

// A remainder may use 1 or 2 slots.
template <typename ItemType, uint64_t kBitsPerItem,
          typename HashFunction = hashutil::TwoIndependentMultiplyShift>
class VEQF
{
  public:
    VEQF(uint64_t max_num_keys)
        : qbits_(__builtin_ctzl(upperpower2(max_num_keys))),
          index_mask_(LowMask(qbits_)),
          entries_(0),
          max_entries_(1ull << qbits_),
          items_(0),
          table_size_(CalcTableSize(qbits_)),
          table_(new uint64_t[table_size_])
    {
        memset(table_.get(), 0, sizeof(uint64_t) * table_size_);
    }

    bool Lookup(const ItemType &key) const
    {
        uint64_t quotient, remainder;
        GenerateQuotientRemainder(key, &quotient, &remainder);

        if (!IsOccupied(GetSlot(quotient)))
        {
            return false;
        }

        uint64_t run_idx{FindRunStart(quotient)}, cur_slot{GetSlot(run_idx)},
            one_slot_remainder{remainder & LowMask(kBitsPerItem)},
            two_slots_first_remainder{(remainder & LowMask(kBitsPerItem - 1)) |
                                      kRemainderHighestBit},
            max_remainder{std::max(one_slot_remainder, two_slots_first_remainder)};
        do
        {
            uint64_t partial_remainder{GetPartialRemainder(cur_slot)}, full_remainder;
            uint64_t step(GetRemainder(run_idx, cur_slot, &full_remainder));
            if ((step == 1 && partial_remainder == one_slot_remainder) ||
                (step == 2 && full_remainder == remainder))
            {
                return true;
            }
            else if (partial_remainder > max_remainder)
            {
                return false;
            }
            run_idx = IncrIdx(run_idx, step);
            cur_slot = GetSlot(run_idx);
        } while (IsContinuation(cur_slot));

        return false;
    }

    bool Insert(const ItemType &key)
    {
        if (items_ >= max_entries_)
        {
            return false;
        }

        uint64_t quotient, remainder,
            slot_count{IsInsertMultipleRemainder() ? kMaxOccupiedSlot : 1ul};
        GenerateQuotientRemainder(key, &quotient, &remainder);
        uint64_t quotient_entry{GetSlot(quotient)},
            to_insert_entry[]{
                (remainder & LowMask(kBitsPerItem)) << kMetadataBits,
                SetContinuation(
                    ((remainder >> kBitsPerItem) & LowMask(kBitsPerItem))
                    << kMetadataBits)};

        if (slot_count == 1)
        {
            to_insert_entry[0] = (remainder & LowMask(kBitsPerItem)) << kMetadataBits;
        }
        else
        {
            to_insert_entry[0] =
                ((remainder & LowMask(kBitsPerItem - 1)) | kRemainderHighestBit)
                << kMetadataBits;
            to_insert_entry[1] = SetContinuation(
                ((remainder >> (kBitsPerItem - 1)) & LowMask(kBitsPerItem - 1))
                << kMetadataBits);
        }

        if (entries_ == max_entries_)
        {
            bool is_quotient_occupied{IsOccupied(quotient_entry)};
            if (!is_quotient_occupied)
            {
                // set is_occupied to find run start
                SetSlot(quotient, SetOccupied(quotient_entry));
            }
            // Before inserting, We should delete the continued remainder of a
            // multiple remainder. Otherwise all remainders in filter are shifted, and
            // cluster start can not located.
            uint64_t multiple_remainder_idx{FindRunStart(quotient)},
                multiple_remainder_quotient{quotient};
            if (!is_quotient_occupied)
            {
                // Find the quotient of entry at run start
                do
                {
                    multiple_remainder_quotient = IncrIdx(multiple_remainder_quotient, 1);
                } while (!IsOccupied(GetSlot(multiple_remainder_quotient)));
                // clear is_occupied, for deleting the continued remainder
                SetSlot(quotient, quotient_entry);
            }
            uint64_t slot;
            do
            {
                multiple_remainder_idx = IncrIdx(multiple_remainder_idx, 1);
                slot = GetSlot(multiple_remainder_idx);
                if (IsRunStart(slot))
                {
                    do
                    {
                        multiple_remainder_quotient =
                            IncrIdx(multiple_remainder_quotient, 1);
                    } while (!IsOccupied(GetSlot(multiple_remainder_quotient)));
                }
            } while (!(IsContinuation(slot) && !IsShifted(slot)));
            AdjustTwoSlotsHighestBit(multiple_remainder_idx);
            DeleteFrom(multiple_remainder_idx, multiple_remainder_quotient,
                       IncrIdx(multiple_remainder_idx, 1));
            --entries_;
            quotient_entry = GetSlot(quotient);
        }

        // Special-case filling canonical slots
        if (IsEmpty(quotient_entry))
        {
            switch (slot_count)
            {
            case 1: {
                SetSlot(quotient, SetOccupied(to_insert_entry[0]));
                ++entries_;
                ++items_;
                return true;
            }
            case 2: {
                uint64_t next_quotient_idx{IncrIdx(quotient, 1)};
                if (!IsEmpty(GetSlot(next_quotient_idx)))
                {
                    // insert an empty entry here, so that 2 slot remainder can be
                    // inserted. ENSURE the empty slot is not less than 2, otherwise we
                    // can not get 2 continuous empty slot.
                    assert(max_entries_ - entries_ >= 2);
                    uint64_t to_insert_entry[]{0};
                    InsertTo(next_quotient_idx, to_insert_entry, 1, true);
                }

                assert(IsEmpty(GetSlot(quotient)) &&
                       (IsEmpty(GetSlot(next_quotient_idx)) ||
                        IsOccupied(GetSlot(next_quotient_idx))));
                if (IsOccupied(GetSlot(next_quotient_idx)))
                {
                    to_insert_entry[1] = SetOccupied(to_insert_entry[1]);
                }
                SetSlot(quotient, SetOccupied(to_insert_entry[0]));
                SetSlot(next_quotient_idx, to_insert_entry[1]);
                entries_ += 2;
                ++items_;
                return true;
            }
            default: {
                assert(false);
                return false;
            }
            }
        }

        bool is_quotient_occupied{IsOccupied(quotient_entry)};
        if (!is_quotient_occupied)
        {
            SetSlot(quotient, SetOccupied(quotient_entry));
        }

        uint64_t run_start{FindRunStart(quotient)};
        uint64_t insert_idx{run_start};

        if (is_quotient_occupied)
        {
            // Move the insert position in the run, if the run has already exists.
            uint64_t cur_slot{GetSlot(insert_idx)},
                one_slot_remainder{GetPartialRemainder(to_insert_entry[0])};
            do
            {
                uint64_t partial_remainder{GetPartialRemainder(cur_slot)},
                    full_remainder;
                uint64_t step(GetRemainder(insert_idx, cur_slot, &full_remainder));
                if (partial_remainder > one_slot_remainder)
                {
                    break;
                }
                insert_idx = IncrIdx(insert_idx, step);
                cur_slot = GetSlot(insert_idx);
            } while (IsContinuation(cur_slot));

            if (insert_idx == run_start)
            {
                // Current run start will be continuation. Also set is_shifted to
                // differentiate from multiple slot remainder.
                SetSlot(run_start, SetShifted(SetContinuation(GetSlot(run_start))));
            }
            else
            {
                // The inserted entry will be continuation.
                to_insert_entry[0] = SetContinuation(to_insert_entry[0]);
            }
        }

        if (insert_idx != quotient)
        {
            to_insert_entry[0] = SetShifted(to_insert_entry[0]);
        }
        entries_ += InsertTo(insert_idx, to_insert_entry, slot_count, false);
        ++items_;
        return true;
    }

    bool Delete(const ItemType &key)
    {
        uint64_t quotient, remainder;
        GenerateQuotientRemainder(key, &quotient, &remainder);
        uint64_t quotient_entry{GetSlot(quotient)};

        if (!IsOccupied(quotient_entry) || entries_ == 0)
        {
            return false;
        }

        uint64_t run_idx{FindRunStart(quotient)}, cur_slot{GetSlot(run_idx)},
            one_slot_remainder{remainder & LowMask(kBitsPerItem)},
            two_slots_first_remainder{(remainder & LowMask(kBitsPerItem - 1)) |
                                      kRemainderHighestBit},
            max_remainder{std::max(one_slot_remainder, two_slots_first_remainder)};
        uint64_t delete_idx{}, delete_step{}, remainder_match_len{};
        do
        {
            uint64_t partial_remainder{GetPartialRemainder(cur_slot)}, full_remainder;
            uint64_t step(GetRemainder(run_idx, cur_slot, &full_remainder));
            if ((step == 1 && partial_remainder == one_slot_remainder) ||
                (step == 2 && full_remainder == remainder))
            {
                if (step > remainder_match_len)
                {
                    delete_idx = run_idx;
                    delete_step = step;
                    remainder_match_len = step;
                    if (step == kMaxOccupiedSlot)
                    {
                        break;
                    }
                }
            }
            else if (partial_remainder > max_remainder)
            {
                break;
            }
            run_idx = IncrIdx(run_idx, step);
            cur_slot = GetSlot(run_idx);
        } while (IsContinuation(cur_slot));

        if (remainder_match_len == 0)
        {
            return false;
        }

        uint64_t delete_entry{delete_idx == quotient ? quotient_entry
                                                     : GetSlot(delete_idx)};
        bool is_run_start{IsRunStart(delete_entry)};
        uint64_t delete_next_idx{IncrIdx(delete_idx, delete_step)};
        uint64_t delete_next_entry{GetSlot(delete_next_idx)};
        if (is_run_start)
        {
            if (!IsContinuation(delete_next_entry))
            {
                // clear `is_occupied` if deleting the last entry in run
                SetSlot(quotient, ClearOccupied(quotient_entry));
            }
        }

        DeleteFrom(delete_idx, quotient, delete_next_idx);

        if (is_run_start)
        {
            uint64_t next{GetSlot(delete_idx)};
            uint64_t updated_next{next};
            if (IsContinuation(updated_next))
            {
                updated_next = ClearContinuation(updated_next);
            }
            if (delete_idx == quotient && IsShifted(updated_next))
            {
                updated_next = ClearShifted(updated_next);
            }
            if (updated_next != next)
            {
                SetSlot(delete_idx, updated_next);
            }
        }
        entries_ -= delete_step;
        --items_;
        return true;
    }

    void SetInsertLargeRemainderThreshold(double threshold)
    {
        insert_large_remainder_threshold_ = threshold;
    }

    size_t Size() const
    {
        return items_;
    }
    size_t SizeInBytes() const
    {
        return table_size_ * sizeof(uint64_t);
    }
    double LoadFactor() const
    {
        return 1.0 * entries_ / max_entries_;
    }
    double BitsPerItem() const
    {
        return 8.0 * SizeInBytes() / Size();
    }

  private:
    static constexpr inline uint64_t LowMask(uint64_t n)
    {
        return (1ull << n) - 1;
    }

    constexpr static uint8_t kIsOccupiedMask = 1;
    constexpr static uint8_t kIsContinuationMask = 2;
    constexpr static uint8_t kIsShiftedMask = 4;
    constexpr static uint8_t kMetadataMask =
        kIsOccupiedMask | kIsContinuationMask | kIsShiftedMask;
    constexpr static uint8_t kMetadataBits{3};
    constexpr static uint64_t kSlotBits{kBitsPerItem + kMetadataBits};
    constexpr static uint64_t kSlotMask{LowMask(kSlotBits)};
    constexpr static uint64_t kMaxOccupiedSlot{2};
    constexpr static uint64_t kRemainderHighestBit{1ull << (kBitsPerItem - 1)};

    // is_continuation, is_shifted = 1, 0 => remainder in multiple slots

    static inline uint64_t CalcTableSize(uint64_t qbits)
    {
        uint64_t total_bits{(1ull << qbits) * kSlotBits};
        return (total_bits + 63) / 64;
    }

    inline void GenerateQuotientRemainder(const ItemType &item,
                                          uint64_t *quotient,
                                          uint64_t *remainder) const
    {
        static_assert(kMaxOccupiedSlot * kBitsPerItem < 64, "no bits for quotient");
        const uint64_t hash = hasher_(item);
        *quotient = (hash >> (kMaxOccupiedSlot * kBitsPerItem - 2)) & index_mask_;
        *remainder = hash & LowMask(kMaxOccupiedSlot * kBitsPerItem - 2);
    }

    uint64_t GetSlot(uint64_t idx) const
    {
        uint64_t bitpos{idx * kSlotBits};
        uint64_t tabpos{bitpos / 64}, slotpos{bitpos % 64};
        int64_t spillbits{static_cast<int64_t>(slotpos + kSlotBits) - 64};
        uint64_t ret{(table_[tabpos] >> slotpos) & kSlotMask};
        if (spillbits > 0)
        {
            ++tabpos; // Move to next 64bit, won't overflow
            uint64_t x{table_[tabpos] & LowMask(spillbits)};
            ret |= x << (kSlotBits - spillbits);
        }
        return ret;
    }

    void SetSlot(uint64_t idx, uint64_t slot)
    {
        uint64_t bitpos{idx * kSlotBits};
        uint64_t tabpos{bitpos / 64}, slotpos{bitpos % 64};
        int64_t spillbits{static_cast<int64_t>(slotpos + kSlotBits) - 64};
        slot &= kSlotMask;
        table_[tabpos] &= ~(kSlotMask << slotpos);
        table_[tabpos] |= slot << slotpos;
        if (spillbits > 0)
        {
            ++tabpos; // Move to next 64bit, won't overflow
            table_[tabpos] &= ~LowMask(spillbits);
            table_[tabpos] |= slot >> (kSlotBits - spillbits);
        }
    }

    inline uint64_t IncrIdx(uint64_t idx, uint64_t step) const
    {
        return (idx + step) & index_mask_;
    }

    inline uint64_t DecrIdx(uint64_t idx) const
    {
        return (idx - 1) & index_mask_;
    }

    inline bool IsOccupied(uint64_t slot) const
    {
        return slot & kIsOccupiedMask;
    }

    inline bool IsContinuation(uint64_t slot) const
    {
        return slot & kIsContinuationMask;
    }

    inline bool IsShifted(uint64_t slot) const
    {
        return slot & kIsShiftedMask;
    }

    inline uint64_t SetOccupied(uint64_t slot) const
    {
        return slot | kIsOccupiedMask;
    }

    inline uint64_t SetContinuation(uint64_t slot) const
    {
        return slot | kIsContinuationMask;
    }

    inline uint64_t SetShifted(uint64_t slot) const
    {
        return slot | kIsShiftedMask;
    }

    inline uint64_t ClearOccupied(uint64_t slot) const
    {
        return slot & ~kIsOccupiedMask;
    }

    inline uint64_t ClearContinuation(uint64_t slot) const
    {
        return slot & ~kIsContinuationMask;
    }

    inline uint64_t ClearShifted(uint64_t slot) const
    {
        return slot & ~kIsShiftedMask;
    }

    inline uint64_t GetPartialRemainder(uint64_t slot) const
    {
        return slot >> kMetadataBits;
    }

    inline bool IsEmpty(uint64_t slot) const
    {
        return (slot & kMetadataMask) == 0;
    }

    inline bool IsClusterStart(uint64_t slot) const
    {
        return IsOccupied(slot) && !IsContinuation(slot) && !IsShifted(slot);
    }

    inline bool IsRunStart(uint64_t slot) const
    {
        return !IsContinuation(slot) && (IsOccupied(slot) || IsShifted(slot));
    }

    // Return how many slots are occupied by the remainder, slot can not be emoty
    uint64_t GetRemainder(uint64_t idx, uint64_t slot,
                          uint64_t *remainder) const
    {
        assert(!IsEmpty(slot));
        *remainder = GetPartialRemainder(slot);
        uint64_t counter{1}, next_idx{IncrIdx(idx, 1)};
        uint64_t next_slot{GetSlot(next_idx)};

        // slot and next_slot should in same run
        if (!IsEmpty(next_slot) && !IsRunStart(next_slot) &&
            (GetPartialRemainder(next_slot) < GetPartialRemainder(slot)))
        {
            assert(IsContinuation(next_slot) && !IsShifted(next_slot));
            *remainder &= LowMask(kBitsPerItem - 1);
            *remainder |= GetPartialRemainder(next_slot) << (kBitsPerItem - 1);
            ++counter;
        }
        assert(counter == 1 || counter == kMaxOccupiedSlot);
        return counter;
    }

    uint64_t FindRunStart(uint64_t quotient) const
    {
        uint64_t cluster_start{quotient};
        // Encounter cluster start before empty
        for (uint64_t cluster_content{GetSlot(cluster_start)};
             IsShifted(cluster_content) || IsContinuation(cluster_content);
             cluster_content = GetSlot(cluster_start))
        {
            cluster_start = DecrIdx(cluster_start);
        }

        uint64_t run_start{cluster_start};
        // IsOccupied(quotient) must be true, otherwise it will be infinite loop
        while (cluster_start != quotient)
        {
            do
            {
                run_start = IncrIdx(run_start, 1);
            } while (IsContinuation(GetSlot(run_start)));

            do
            {
                cluster_start = IncrIdx(cluster_start, 1);
            } while (!IsOccupied(GetSlot(cluster_start)));
        }
        return run_start;
    }

    // For entry insert
    class Queue
    {
      public:
        Queue() = default;

        bool IsEmpty() const
        {
            return first_ == last_;
        }
        bool IsFull() const
        {
            return last_ - first_ == kMaxOccupiedSlot || first_ - last_ == 1;
        }
        void Enqueue(uint64_t entry)
        {
            assert(!IsFull());
            slots_[last_] = entry;
            last_ = (last_ + 1) % (kMaxOccupiedSlot + 1);
        }
        uint64_t Dequeue()
        {
            assert(!IsEmpty());
            uint64_t ret{slots_[first_]};
            first_ = (first_ + 1) % (kMaxOccupiedSlot + 1);
            return ret;
        }

      private:
        uint64_t slots_[kMaxOccupiedSlot + 1];
        uint64_t first_{}, last_{};
    };

    // Return how many empty slots are used
    uint64_t InsertTo(uint64_t insert_idx, uint64_t *to_insert_entry,
                      uint64_t slot_count, bool force_disable_compaction)
    {
        Queue q;
        uint64_t prev;
        uint64_t curr;
        bool is_empty;
        uint64_t ret{slot_count};
        for (uint64_t i = 0; i < slot_count; ++i)
        {
            q.Enqueue(to_insert_entry[i]);
        }

        do
        {
            prev = GetSlot(insert_idx);
            curr = q.Dequeue();
            is_empty = IsEmpty(prev);
            bool need_move_backwards{false};

            if (!is_empty)
            {
                bool is_multiple_remainder{IsContinuation(prev) && !IsShifted(prev)};
                if (!is_multiple_remainder)
                {
                    // for remainders in multiple slots, don't set their is_shifted bit
                    prev = SetShifted(prev);
                }
                if (IsOccupied(prev))
                {
                    prev = ClearOccupied(prev);
                    curr = SetOccupied(curr);
                }
                bool reuse_multiple_remainder{!force_disable_compaction &&
                                              !IsInsertMultipleRemainder() &&
                                              is_multiple_remainder};
                if (!reuse_multiple_remainder)
                {
                    // otherwise reuse the remainder countinuation space
                    q.Enqueue(prev);
                }
                else
                {
                    if (GetPartialRemainder(prev) % 2 == 0)
                    {
                        curr &= ~(1 << (kSlotBits - 1));
                        need_move_backwards = true;
                    }
                    --ret;
                }
            }
            SetSlot(insert_idx, curr);
            if (need_move_backwards)
            {
                MoveCompactedSlot(insert_idx, curr);
            }
            insert_idx = IncrIdx(insert_idx, 1);
        } while (!q.IsEmpty());

        return ret;
    }

    void DeleteFrom(uint64_t delete_idx, uint64_t quotient,
                    uint64_t delete_next_idx)
    {
        uint64_t delete_curr_entry{GetSlot(delete_idx)};
        uint64_t delete_next_entry;
        const uint64_t orig_delete_idx{delete_idx};

        while (true)
        {
            delete_next_entry = GetSlot(delete_next_idx);
            bool curr_occupied{IsOccupied(delete_curr_entry)};

            if (IsEmpty(delete_next_entry) || IsClusterStart(delete_next_entry) ||
                delete_next_idx == orig_delete_idx)
            {
                for (uint64_t i = delete_idx; i != delete_next_idx; i = IncrIdx(i, 1))
                {
                    SetSlot(i, 0);
                }
                return;
            }

            uint64_t updated_next{delete_next_entry};
            // Fix entries which slide into canonical slots
            if (IsRunStart(updated_next))
            {
                do
                {
                    quotient = IncrIdx(quotient, 1);
                } while (!IsOccupied(GetSlot(quotient)));

                while (
                    IsDeleteIdxBeforeQuotient(delete_idx, quotient, delete_next_idx))
                {
                    // can not move remainder before quotient, padding zero
                    uint64_t padding_slot{};
                    if (curr_occupied)
                    {
                        padding_slot = SetOccupied(padding_slot);
                    }
                    SetSlot(delete_idx, padding_slot);
                    delete_idx = IncrIdx(delete_idx, 1);
                    delete_curr_entry = GetSlot(delete_idx);
                    curr_occupied = IsOccupied(delete_curr_entry);
                }

                if (curr_occupied && quotient == delete_idx)
                {
                    updated_next = ClearShifted(updated_next);
                }
            }

            SetSlot(delete_idx, curr_occupied ? SetOccupied(updated_next)
                                              : ClearOccupied(updated_next));
            delete_idx = IncrIdx(delete_idx, 1);
            delete_curr_entry = GetSlot(delete_idx);
            delete_next_idx = IncrIdx(delete_next_idx, 1);
            delete_next_entry = GetSlot(delete_next_idx);
        }
    }

    bool IsDeleteIdxBeforeQuotient(uint64_t delete_idx, uint64_t next_quotient,
                                   uint64_t delete_next_idx)
    {
        // If everyone is not overflow,
        // `delete_idx < next_quotient && next_quotient < delete_next_idx` is
        // enough. However, next_quotient and delete_next_idx may overflow, so
        // another two equation is added.
        return (delete_idx < next_quotient && next_quotient < delete_next_idx) ||
               (delete_next_idx < delete_idx && delete_idx < next_quotient) ||
               (next_quotient < delete_next_idx && delete_next_idx < delete_idx);
    }

    inline bool IsInsertMultipleRemainder()
    {
        return entries_ < max_entries_ * insert_large_remainder_threshold_;
    }

    void AdjustTwoSlotsHighestBit(uint64_t multiple_remainder_second_idx)
    {
        uint64_t multiple_remainder_first_idx{
            DecrIdx(multiple_remainder_second_idx)};
        uint64_t first_slot{GetSlot(multiple_remainder_first_idx)},
            second_slot{GetSlot(multiple_remainder_second_idx)};
        if (GetPartialRemainder(second_slot) % 2 == 0)
        {
            // two slots compacted to one slot, and value changed.
            // since slot value gets smaller, find new position backward.
            first_slot &= ~(1 << (kSlotBits - 1));
            MoveCompactedSlot(multiple_remainder_first_idx, first_slot);
        }
    }

    void MoveCompactedSlot(uint64_t multiple_remainder_first_idx,
                           uint64_t first_slot)
    {
        if (IsRunStart(first_slot))
        {
            // already run start, no need to find new position.
            SetSlot(multiple_remainder_first_idx, first_slot);
            return;
        }
        uint64_t new_position_idx{multiple_remainder_first_idx},
            new_position_slot{first_slot};
        {
            // find new position
            uint64_t curr_idx{DecrIdx(multiple_remainder_first_idx)}, curr_slot,
                first_slot_remainder{GetPartialRemainder(first_slot)};
            do
            {
                curr_slot = GetSlot(curr_idx);
                if (!(IsContinuation(curr_slot) && !IsShifted(curr_slot)))
                {
                    if (GetPartialRemainder(curr_slot) >= first_slot_remainder)
                    {
                        new_position_idx = curr_idx;
                        new_position_slot = curr_slot;
                    }
                    else
                    {
                        break;
                    }
                }
                curr_idx = DecrIdx(curr_idx);
            } while (!IsRunStart(curr_slot));
        }
        if (new_position_idx != multiple_remainder_first_idx)
        {
            bool multiple_remainder_first_slot_is_occupied{IsOccupied(first_slot)};
            // Set metadata bits for first_slot
            if (IsOccupied(new_position_slot))
            {
                first_slot = SetOccupied(first_slot);
            }
            else
            {
                first_slot = ClearOccupied(first_slot);
            }
            if (IsClusterStart(new_position_slot))
            {
                first_slot = ClearContinuation(first_slot);
                first_slot = ClearShifted(first_slot);
                // Current run start will be continuation. Also set is_shifted to
                // differentiate from multiple slot remainder.
                SetSlot(new_position_idx,
                        SetShifted(SetContinuation(new_position_slot)));
            }
            else if (IsRunStart(new_position_slot))
            {
                first_slot = ClearContinuation(first_slot);
                first_slot = SetShifted(first_slot);
                // Current run start will be continuation. Also set is_shifted to
                // differentiate from multiple slot remainder.
                SetSlot(new_position_idx,
                        SetShifted(SetContinuation(new_position_slot)));
            }
            else
            {
                first_slot = SetContinuation(first_slot);
                first_slot = SetShifted(first_slot);
            }
            // Make a empty slot, as insert end
            SetSlot(multiple_remainder_first_idx, 0);
            InsertTo(new_position_idx, &first_slot, 1, true);
            if (multiple_remainder_first_slot_is_occupied)
            {
                SetSlot(multiple_remainder_first_idx,
                        SetOccupied(GetSlot(multiple_remainder_first_idx)));
            }
        }
        else
        {
            SetSlot(new_position_idx, first_slot);
        }
    }

    uint8_t qbits_;
    uint64_t index_mask_, entries_, // count of occupied slots
        max_entries_, items_,       // count of inserted items
        table_size_;
    HashFunction hasher_;
    std::unique_ptr<uint64_t[]> table_;
    double insert_large_remainder_threshold_{0.2};
};

}

#endif