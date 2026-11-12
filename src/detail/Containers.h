#pragma once

#include <algorithm>
#include <bitset>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace mesh_arrange::detail
{

/// @brief Multi-word dynamic bitset supporting fast `.set()`, `.test()`, `.count()`, and `operator<`.
class DynamicBitset
{
  public:
    DynamicBitset() = default;

    explicit DynamicBitset(std::size_t n_bits)
    {
        resize(n_bits);
    }

    void resize(std::size_t n_bits)
    {
        const std::size_t n_blocks = (n_bits + 63) / 64;
        blocks_.assign(n_blocks, 0ULL);
    }

    void set(std::size_t bit, bool value = true)
    {
        const std::size_t blk = bit / 64;
        if (blk >= blocks_.size())
        {
            blocks_.resize(blk + 1, 0ULL);
        }
        const uint64_t mask = uint64_t{ 1 } << (bit % 64);
        if (value)
        {
            blocks_[blk] |= mask;
        }
        else
        {
            blocks_[blk] &= ~mask;
        }
    }

    bool test(std::size_t bit) const noexcept
    {
        const std::size_t blk = bit / 64;
        if (blk >= blocks_.size())
        {
            return false;
        }
        return (blocks_[blk] & (uint64_t{ 1 } << (bit % 64))) != 0ULL;
    }

    std::size_t count() const noexcept
    {
        std::size_t c = 0;
        for (uint64_t w : blocks_)
        {
            c += std::bitset<64>(w).count();
        }
        return c;
    }

    bool operator==(const DynamicBitset& other) const noexcept
    {
        return blocks_ == other.blocks_;
    }

    bool operator<(const DynamicBitset& other) const noexcept
    {
        if (blocks_.size() != other.blocks_.size())
        {
            return blocks_.size() < other.blocks_.size();
        }
        for (std::size_t i = blocks_.size(); i > 0; --i)
        {
            const std::size_t k = i - 1;
            if (blocks_[k] != other.blocks_[k])
            {
                return blocks_[k] < other.blocks_[k];
            }
        }
        return false;
    }

  private:
    std::vector<uint64_t> blocks_;
};

/// @brief Small flat deduplicated vector set.
template<typename T>
class VectorSet
{
  public:
    bool insert(const T& val)
    {
        if (std::find(v_.begin(), v_.end(), val) != v_.end())
        {
            return false;
        }
        v_.push_back(val);
        return true;
    }

    std::size_t size() const noexcept { return v_.size(); }
    auto begin() const noexcept { return v_.begin(); }
    auto end() const noexcept { return v_.end(); }

  private:
    std::vector<T> v_;
};

} // namespace mesh_arrange::detail
