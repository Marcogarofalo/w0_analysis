#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <string>
#include <stdexcept>

std::pair<std::vector<size_t>, std::vector<size_t>>
matching_indices(const std::vector<std::string>& a,
    const std::vector<std::string>& b) {
    
    // ---- uniqueness check for a ----
    {
        std::unordered_set<std::string> seen;
        seen.reserve(a.size());
        for (const auto& s : a) {
            if (!seen.insert(s).second) {
                throw std::runtime_error("Duplicate entry in vector a: " + s);
            }
        }
    }

    // ---- uniqueness check for b ----
    {
        std::unordered_set<std::string> seen;
        seen.reserve(b.size());
        for (const auto& s : b) {
            if (!seen.insert(s).second) {
                throw std::runtime_error("Duplicate entry in vector b: " + s);
            }
        }
    }

    // ---- map a -> index ----
    std::unordered_map<std::string, size_t> pos;
    pos.reserve(a.size());

    for (size_t i = 0; i < a.size(); ++i) {
        pos.emplace(a[i], i);
    }

    // ---- match ----
    std::vector<size_t> idx_a;
    std::vector<size_t> idx_b;
    idx_a.reserve(std::min(a.size(), b.size()));
    idx_b.reserve(std::min(a.size(), b.size()));

    for (size_t j = 0; j < b.size(); ++j) {
        auto it = pos.find(b[j]);
        if (it != pos.end()) {
            idx_a.push_back(it->second);
            idx_b.push_back(j);
        }
    }

    return { idx_a, idx_b };
}

// int main() {
//     std::vector<std::string> a = { "a", "b", "c", "d", "e" };
//     std::vector<std::string> b = { "e", "c", "b" };

//     auto [ia, ib] = matching_indices(a, b);
//     for (size_t k = 0; k < ia.size(); ++k) {
//         printf("a[%zu] = %s matches b[%zu] = %s\n", ia[k], a[ia[k]].c_str(), ib[k], b[ib[k]].c_str());
//     }
//     return 0;

// }