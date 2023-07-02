#include <gtest/gtest.h>

#include <cstdint>

#include "vecbf/vecbf.h"
#include "vecf/vecf.h"
#include "veqf/veqf.h"

template <typename T>
class VEFrameworkTest : public testing::Test
{
  protected:
    VEFrameworkTest()
        : filter_(total_items)
    {
    }
    ~VEFrameworkTest() = default;

    constexpr static uint64_t total_items = 1024 * 1024 * 7;

    T filter_;
};

using Implementations =
    testing::Types<vecf::VECF<uint64_t, 8>, vecf::VECF<uint64_t, 12>,
                   vecf::VECF<uint64_t, 16>, veqf::VEQF<uint64_t, 8>,
                   veqf::VEQF<uint64_t, 10>, veqf::VEQF<uint64_t, 12>,
                   veqf::VEQF<uint64_t, 14>, veqf::VEQF<uint64_t, 16>,
                   vecbf::VECBF<uint64_t, 8>>;
TYPED_TEST_SUITE(VEFrameworkTest, Implementations);

TYPED_TEST(VEFrameworkTest, Correctness)
{
    // Insert items to this cuckoo filter
    uint64_t num_inserted = 0;
    for (uint64_t i = 0; i < this->total_items; i++, num_inserted++)
    {
        if (!this->filter_.Insert(i))
        {
            break;
        }
    }
    printf("Total items count: %lu, inserted items count: %lu\n",
           this->total_items, num_inserted);

    // Check if previously inserted items are in the filter, expected true for all
    // items
    for (uint64_t i = 0; i < num_inserted; i++)
    {
        ASSERT_TRUE(this->filter_.Lookup(i));
    }

    // Check non-existing items, a few false positives expected
    uint64_t total_queries = 0, false_queries = 0;
    for (uint64_t i = this->total_items; i < 2 * this->total_items; i++)
    {
        if (this->filter_.Lookup(i))
        {
            false_queries++;
        }
        total_queries++;
    }

    // Output the measured false positive rate
    printf("false positive rate is %f %%\n",
           100.0 * false_queries / total_queries);

    // Check if previously inserted items can be removed, expected true for all
    // items
    for (uint64_t i = 0; i < num_inserted; i++)
    {
        ASSERT_TRUE(this->filter_.Delete(i));
    }
}