/**
 * @file maximal_clique_hashmap.cpp
 * @brief Implementation: MCE Degeneracy (Eppstein, Löffler, Strash 2011)
 *
 * Triển khai đúng theo thuật toán BronKerboschDegeneracy được mô tả trong paper:
 * "Listing All Maximal Cliques in Sparse Graphs in Near-Optimal Time"
 * David Eppstein, Maarten Löffler, Darren Strash — Dagstuhl 2011.
 *
 * Cấu trúc thuật toán (theo Figure 4 trong paper):
 *
 *   proc BronKerboschDegeneracy(V, E)
 *     for each vertex v_i in a degeneracy ordering v0, v1, v2, ... of (V, E) do
 *       P ← Γ(v_i) ∩ { v_{i+1}, ..., v_{n-1} }   // later neighbors
 *       X ← Γ(v_i) ∩ { v_0,    ..., v_{i-1} }     // earlier neighbors
 *       BronKerboschPivot(P, {v_i}, X)
 *     end for
 *
 *   proc BronKerboschPivot(P, R, X)                  // Tomita et al. pivot rule
 *     if P ∪ X = ∅ then report R as a maximal clique
 *     choose pivot u ∈ P ∪ X maximizing |P ∩ Γ(u)|
 *     for each v ∈ P \ Γ(u) do
 *       BronKerboschPivot(P ∩ Γ(v), R ∪ {v}, X ∩ Γ(v))
 *       P ← P \ {v}
 *       X ← X ∪ {v}
 *     end for
 *
 * Độ phức tạp thời gian: O(d · n · 3^{d/3})  (Theorem 2 trong paper)
 * với d = degeneracy của đồ thị, n = số đỉnh.
 */

#include "maximal_clique_hashmap.h"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>

namespace {

    using Node = const SpatialInstance*;
    using CliqueVec = std::vector<Node>;
    using AdjMap = std::unordered_map<Node, CliqueVec>;

    // Kiểu dữ liệu kết quả
    using ResultMap = std::map<
        Colocation,
        std::unordered_map<FeatureType, std::set<const SpatialInstance*>>
    >;

    // =========================================================================
    // HELPER: Set Operations (tất cả CliqueVec đều được giữ sorted theo pointer)
    // =========================================================================

    /**
     * Đếm |A ∩ B| — dùng merge scan vì cả hai đều sorted.
     */
    int count_intersection(const CliqueVec& A, const CliqueVec& B) {
        int count = 0;
        auto it1 = A.begin();
        auto it2 = B.begin();
        while (it1 != A.end() && it2 != B.end()) {
            if (*it1 < *it2) ++it1;
            else if (*it2 < *it1) ++it2;
            else { ++count; ++it1; ++it2; }
        }
        return count;
    }

    /**
     * Trả về A \ B  (các phần tử trong A mà không có trong B).
     * Yêu cầu A và B đều đã sorted.
     */
    CliqueVec set_difference_helper(const CliqueVec& A, const CliqueVec& B) {
        CliqueVec result;
        result.reserve(A.size());
        std::set_difference(
            A.begin(), A.end(),
            B.begin(), B.end(),
            std::back_inserter(result));
        return result;
    }

    /**
     * Trả về A ∩ B.
     * Yêu cầu A và B đều đã sorted.
     */
    CliqueVec set_intersection_helper(const CliqueVec& A, const CliqueVec& B) {
        CliqueVec result;
        result.reserve(std::min(A.size(), B.size()));
        std::set_intersection(
            A.begin(), A.end(),
            B.begin(), B.end(),
            std::back_inserter(result));
        return result;
    }

    // =========================================================================
    // OUTPUT: Lưu clique vào ResultMap
    // =========================================================================

    /**
     * Ghi nhận một maximal clique R vào hashMap.
     * Chỉ lưu nếu |R| >= 2 (co-location pattern yêu cầu ít nhất 2 feature type).
     */
    void report_clique(const CliqueVec& R, ResultMap& hashMap) {
        if (R.size() < 2) return;

        // Tạo khoá Colocation từ danh sách feature-type, đã sắp xếp
        Colocation colocationKey;
        colocationKey.reserve(R.size());
        for (const auto& instancePtr : R) {
            colocationKey.push_back(instancePtr->type);
        }
        std::sort(colocationKey.begin(), colocationKey.end());

        // Ghi instance vào inner map theo feature type
        auto& innerMap = hashMap[colocationKey];
        for (const auto& instancePtr : R) {
            innerMap[instancePtr->type].insert(instancePtr);
        }
    }

    // =========================================================================
    // ALGORITHM 1: BronKerboschPivot  (Fig. 2 — phiên bản Tomita et al.)
    //
    // Tham số:
    //   R   — partial clique hiện tại (passed by value để backtracking tự nhiên)
    //   P   — ứng viên có thể mở rộng clique (passed by value)
    //   X   — các đỉnh đã xử lý / cần loại trừ (passed by value)
    //   adj — adjacency map toàn cục (by const ref)
    //   hashMap — nơi lưu kết quả (by ref)
    // =========================================================================
    void runBKPivot(
        CliqueVec       R,          // current clique
        CliqueVec       P,          // candidates
        CliqueVec       X,          // excluded
        const AdjMap& adj,
        ResultMap& hashMap)
    {
        // Base case: P ∪ X = ∅  →  R là maximal clique
        if (P.empty() && X.empty()) {
            report_clique(R, hashMap);
            return;
        }

        // Nếu P rỗng nhưng X khác rỗng: R không tối đại (bị chặn bởi X), dừng
        if (P.empty()) return;

        // ------------------------------------------------------------------
        // Bước 1: Chọn pivot u ∈ P ∪ X sao cho |P ∩ Γ(u)| lớn nhất
        //         (Tomita et al. pivot rule — đảm bảo O(3^{n/3}) worst-case)
        // ------------------------------------------------------------------
        Node u_pivot = nullptr;
        int  max_inter = -1;

        auto try_pivot = [&](Node candidate) {
            auto it = adj.find(candidate);
            if (it != adj.end()) {
                int inter = count_intersection(P, it->second);
                if (inter > max_inter) {
                    max_inter = inter;
                    u_pivot = candidate;
                }
            }
            };

        for (Node node : P) try_pivot(node);
        for (Node node : X) try_pivot(node);

        // ------------------------------------------------------------------
        // Bước 2: Tập ứng viên = P \ Γ(u_pivot)
        //         Các đỉnh trong P ∩ Γ(u) bị trì hoãn sang các lần gọi đệ quy sau
        // ------------------------------------------------------------------
        CliqueVec candidates;
        if (u_pivot != nullptr) {
            auto it = adj.find(u_pivot);
            candidates = (it != adj.end())
                ? set_difference_helper(P, it->second)
                : P;
        }
        else {
            candidates = P;
        }

        // ------------------------------------------------------------------
        // Bước 3: Duyệt qua từng v ∈ candidates và gọi đệ quy
        // ------------------------------------------------------------------
        for (Node v : candidates) {

            // Lấy hàng xóm của v (sorted)
            const CliqueVec* neighbors_v_ptr = nullptr;
            static const CliqueVec empty_vec;
            auto it = adj.find(v);
            neighbors_v_ptr = (it != adj.end()) ? &it->second : &empty_vec;
            const CliqueVec& neighbors_v = *neighbors_v_ptr;

            // Gọi đệ quy:  R ∪ {v},  P ∩ Γ(v),  X ∩ Γ(v)
            CliqueVec newR = R;
            newR.push_back(v);

            runBKPivot(
                std::move(newR),
                set_intersection_helper(P, neighbors_v),
                set_intersection_helper(X, neighbors_v),
                adj,
                hashMap);

            // Backtrack: chuyển v từ P sang X
            auto itP = std::lower_bound(P.begin(), P.end(), v);
            if (itP != P.end() && *itP == v) P.erase(itP);

            auto itX = std::lower_bound(X.begin(), X.end(), v);
            X.insert(itX, v);
        }
    }

    // =========================================================================
    // ALGORITHM 2: Degeneracy Ordering  (Lick & White 1970; O(n + m) linear)
    //
    // Trả về thứ tự suy biến: mỗi đỉnh có tối đa d hàng xóm xuất hiện SAU nó.
    // Degeneracy d = max degree tại thời điểm loại bỏ.
    //
    // Triển khai dùng set<pair<degree, Node>> để luôn lấy đỉnh có bậc nhỏ nhất.
    // Độ phức tạp: O((n + m) log n) — đủ đúng đắn cho mọi đồ thị thực tế.
    // =========================================================================
    std::vector<Node> getDegeneracyOrdering(const AdjMap& adj) {
        std::unordered_map<Node, int> degrees;
        degrees.reserve(adj.size());

        // Min-heap giả lập bằng ordered set: (degree, Node)
        std::set<std::pair<int, Node>> sortedNodes;

        for (const auto& entry : adj) {
            int d = static_cast<int>(entry.second.size());
            degrees[entry.first] = d;
            sortedNodes.insert({ d, entry.first });
        }

        std::vector<Node> ordering;
        ordering.reserve(adj.size());

        std::unordered_map<Node, bool> removed;
        removed.reserve(adj.size());

        while (!sortedNodes.empty()) {
            // Lấy đỉnh có bậc hiện tại nhỏ nhất
            auto it = sortedNodes.begin();
            Node u = it->second;
            sortedNodes.erase(it);

            ordering.push_back(u);
            removed[u] = true;

            // Giảm bậc của các hàng xóm chưa bị loại bỏ
            auto itAdj = adj.find(u);
            if (itAdj == adj.end()) continue;

            for (Node v : itAdj->second) {
                if (removed.count(v)) continue;

                int old_deg = degrees[v];
                sortedNodes.erase({ old_deg, v });
                degrees[v] = old_deg - 1;
                sortedNodes.insert({ old_deg - 1, v });
            }
        }

        return ordering;
    }

} // anonymous namespace

// =============================================================================
// PUBLIC METHOD: executeBK
//
// Triển khai proc BronKerboschDegeneracy(V, E)  (Figure 4 trong paper):
//
//   for each vertex v_i in degeneracy ordering do
//       P ← Γ(v_i) ∩ { later neighbors }
//       X ← Γ(v_i) ∩ { earlier neighbors }
//       BronKerboschPivot(P, {v_i}, X)
//   end for
// =============================================================================
std::map<Colocation, std::unordered_map<FeatureType, std::set<const SpatialInstance*>>>
MaximalCliqueHashmap::executeBK(const std::vector<NeighborSet>& neighborSets)
{
    // ------------------------------------------------------------------
    // Bước 1: Xây dựng adjacency map
    //         Mỗi neighbor list được sort theo pointer để set ops hoạt động.
    // ------------------------------------------------------------------
    AdjMap adj;
    adj.reserve(neighborSets.size());

    for (const auto& ns : neighborSets) {
        Node u = ns.center;
        CliqueVec sorted_neighbors = ns.neighbors;
        std::sort(sorted_neighbors.begin(), sorted_neighbors.end());
        adj[u] = std::move(sorted_neighbors);
    }

    // ------------------------------------------------------------------
    // Bước 2: Tính degeneracy ordering
    // ------------------------------------------------------------------
    std::vector<Node> ordering = getDegeneracyOrdering(adj);

    // Map từ Node → vị trí index trong ordering (để phân loại earlier/later)
    std::unordered_map<Node, int> orderIndex;
    orderIndex.reserve(ordering.size());
    for (int i = 0; i < static_cast<int>(ordering.size()); ++i) {
        orderIndex[ordering[i]] = i;
    }

    // ------------------------------------------------------------------
    // Bước 3: Vòng lặp ngoài — duyệt theo degeneracy ordering
    //         Đây là outer loop của BronKerboschDegeneracy (Figure 4).
    //
    //   P ← Γ(v_i) ∩ { v_{i+1}, ..., v_{n-1} }  (later neighbors)
    //   X ← Γ(v_i) ∩ { v_0,    ..., v_{i-1} }   (earlier neighbors)
    //   BronKerboschPivot(P, {v_i}, X)
    //
    //   Nhờ degeneracy ordering, |P| ≤ d tại mọi outer call,
    //   đảm bảo tổng thời gian O(d · n · 3^{d/3}).
    // ------------------------------------------------------------------
    ResultMap hashMap;

    for (int i = 0; i < static_cast<int>(ordering.size()); ++i) {
        Node v = ordering[i];

        auto itAdj = adj.find(v);
        if (itAdj == adj.end()) continue;
        const CliqueVec& neighbors = itAdj->second;

        CliqueVec P, X;
        P.reserve(neighbors.size());
        X.reserve(neighbors.size());

        for (Node nb : neighbors) {
            auto idxIt = orderIndex.find(nb);
            if (idxIt == orderIndex.end()) continue;

            if (idxIt->second > i)
                P.push_back(nb);   // later neighbor  → P
            else
                X.push_back(nb);   // earlier neighbor → X
        }

        // P và X phải sorted (set ops trong BronKerboschPivot yêu cầu điều này)
        std::sort(P.begin(), P.end());
        std::sort(X.begin(), X.end());

        // Gọi BronKerboschPivot với R = {v}
        runBKPivot(
            CliqueVec{ v },     // R = {v_i}
            std::move(P),
            std::move(X),
            adj,
            hashMap);
    }

    return hashMap;
}

// =============================================================================
// PUBLIC METHOD: extractInitialCandidates
// Trích xuất tất cả co-location patterns từ hashMap vào priority queue.
// =============================================================================
std::priority_queue<Colocation, std::vector<Colocation>, ColocationPriorityComp>
MaximalCliqueHashmap::extractInitialCandidates(
    const std::map<Colocation,
    std::unordered_map<FeatureType, std::set<const SpatialInstance*>>>& hashMap)
{
    std::priority_queue<Colocation, std::vector<Colocation>, ColocationPriorityComp> candidateQueue;
    for (const auto& entry : hashMap) {
        candidateQueue.push(entry.first);
    }
    return candidateQueue;
}