#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>
#include <unordered_map>
#include <unordered_set>

template<class T>
struct TMin {
    const T& operator()(const T& a, const T& b) const noexcept {
        return std::min(a, b);
    }
};

template<class T>
struct TMax {
    const T& operator()(const T& a, const T& b) const noexcept {
        return std::max(a, b);
    }
};

template<class T, class TOp>
class TRmq {
public:
    TRmq(const std::vector<T>& data)
        : Shifts(1, 0)
    {
        const auto n = data.size();

        Modifications.reserve(n);

        {
            auto sz = n;
            for (auto cut = 2u; cut < n + 1u; cut <<= 1) {
                sz += n + 1 - cut;
            }

            SparseTable.reserve(sz);
            SparseTable.resize(n);
            Sizes.push_back(n);
            std::copy_n(data.begin(), n, SparseTable.begin());
        }

        for (auto i = 1u, cut = 2u; cut < n + 1u; ++i, cut <<= 1) {
            Shifts.push_back(SparseTable.size());

            const auto size = n + 1 - cut;
            Sizes.push_back(size);
            SparseTable.insert(SparseTable.end(), size, 0);

            const auto cur = SparseTable.begin() + Shifts.back();
            const auto prev = SparseTable.begin() + Shifts[Shifts.size() - 2];

            const auto shift = (cut >> 1);

            for (auto j = 0u; j < size; ++j) {
                cur[j] = TOp()(prev[j], prev[j + shift]);
            }
        }

        Log2.resize(n + 1);
        for (auto i = 0u, power = 2u, exponent = 0u; i < Log2.size(); ++i) {
            if (i == power) {
                ++exponent;
                power <<= 1;
            }
            Log2[i] = exponent;
        }
    }

    const auto& operator[](const size_t i) const {
        return SparseTable.front()[i];
    }

    auto Get(const size_t i, const size_t j) {
        assert(i <= j);

        FlushModifications();

        const auto k = Log2[j + 1 - i];
        const auto pos = SparseTable.begin() + Shifts[k];
        const auto& a = pos[i];
        const auto& b = pos[j + 1 - (1 << k)];

        return TOp()(a, b);
    }

    void Modify(const size_t pos, T val) {
        assert(pos < Sizes.at(0));

        Modifications[pos] = std::move(val);
    }

private:
    void FlushModifications() {
        if (Modifications.empty()) {
            return;
        }

        auto& prevLevelIdxs = Sets[0] = {};
        for (auto& data : Modifications) {
            SparseTable[data.first] = std::move(data.second);
            prevLevelIdxs.insert(data.first);
        }

        auto& curIdxs = Sets[1];

        for (auto i = 1u, shift = 1u; i < Shifts.size() && !prevLevelIdxs.empty(); ++i, shift <<= 1) {
            curIdxs.clear();
            for (const auto idx : prevLevelIdxs) {
                if (idx < Sizes[i]) {
                    curIdxs.insert(idx);
                }

                if (idx >= shift) {
                    curIdxs.insert(idx - shift);
                }
            }

            const auto cur = SparseTable.begin() + Shifts[i];
            const auto prev = SparseTable.begin() + Shifts[i - 1];

            for (const auto j : curIdxs) {
                cur[j] = TOp()(prev[j], prev[j + shift]);
            }
            std::swap(prevLevelIdxs, curIdxs);
        }

        Modifications.clear();
    }

private:
    std::unordered_map<size_t, T> Modifications;
    std::vector<T> SparseTable;
    std::array<std::unordered_set<size_t>, 2> Sets;
    std::vector<size_t> Shifts;
    std::vector<size_t> Sizes;
    std::vector<size_t> Log2;
};

int main() {
    using namespace std;

    vector<int> data(50);
    srand(777);
    for (auto i = 0u; i < data.size(); ++i) {
        data[i] = -100 + (rand() % 200);
    }

    using TOp = TMin<int>;
    TRmq<int, TOp> st(data);

    for (auto i = 0u; i + 1 < data.size(); ++i) {
        for (auto j = i + 1; j < data.size(); ++j) {

            for (auto k = 0u; k < 1 + rand() % 9; ++k) {
                const auto modIdx = i + (rand() % (j - i));

                data[modIdx] = -100 + (rand() % 200);
                st.Modify(modIdx, data[modIdx]);
            }

            const auto trueSum = std::accumulate(data.begin() + i, data.begin() + j + 1, data[i], TOp());
            const auto gotSum = st.Get(i, j);
            assert(trueSum == gotSum);
        }
    }

    return EXIT_SUCCESS;
}
