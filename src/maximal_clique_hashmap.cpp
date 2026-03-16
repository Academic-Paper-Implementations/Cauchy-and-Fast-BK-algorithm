/**
 * @file maximal_clique_hashmap.cpp
 * @brief Implementation: Hybrid MCE (Degeneracy + RCD + Pivot)
 * Based on "Fast Maximal Clique Enumeration for Real-World Graphs"
 */

#include "maximal_clique_hashmap.h"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <cmath> // For floor/ceil if needed

namespace {

    using Node = const SpatialInstance*;
    using CliqueVec = std::vector<Node>;
    using AdjMap = std::unordered_map<Node, CliqueVec>;

    // Type definition for the result map structure
    using ResultMap = std::map<Colocation, std::unordered_map<FeatureType, std::set<const SpatialInstance*>>>;

    // --- HELPER FUNCTIONS (Set Operations) ---

    // Đếm số phần tử chung (Intersection Size)
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

    // P \ N(u)
    CliqueVec set_difference_helper(const CliqueVec& A, const CliqueVec& B) {
        CliqueVec result;
        result.reserve(A.size());
        std::set_difference(A.begin(), A.end(), B.begin(), B.end(), std::back_inserter(result));
        return result;
    }

    // P intersection N(u)
    CliqueVec set_intersection_helper(const CliqueVec& A, const CliqueVec& B) {
        CliqueVec result;
        result.reserve(std::min(A.size(), B.size()));
        std::set_intersection(A.begin(), A.end(), B.begin(), B.end(), std::back_inserter(result));
        return result;
    }

    // Hàm lưu kết quả vào Hashmap
    void report_clique(const CliqueVec& R, ResultMap& hashMap) {
        if (R.size() < 2) return;

        Colocation colocationKey;
        colocationKey.reserve(R.size());
        for (const auto& instancePtr : R) {
            colocationKey.push_back(instancePtr->type);
        }
        std::sort(colocationKey.begin(), colocationKey.end());

        auto& innerMap = hashMap[colocationKey];
        for (const auto& instancePtr : R) {
            innerMap[instancePtr->type].insert(instancePtr);
        }
    }

    // --- ALGORITHM 1: BK PIVOT (Standard) ---
    void runBKPivot(
        CliqueVec R,
        CliqueVec P,
        CliqueVec X,
        const AdjMap& adj,
        ResultMap& hashMap)
    {
        if (P.empty() && X.empty()) {
            report_clique(R, hashMap);
            return;
        }
        if (P.empty()) return;

        // 1. Select Pivot u in P U X maximizing |P n N(u)|
        Node u_pivot = nullptr;
        int max_inter = -1;

        auto check_pivot = [&](Node candidate) {
            auto it = adj.find(candidate);
            if (it != adj.end()) {
                int inter_size = count_intersection(P, it->second);
                if (inter_size > max_inter) {
                    max_inter = inter_size;
                    u_pivot = candidate;
                }
            }
            };

        for (Node node : P) check_pivot(node);
        for (Node node : X) check_pivot(node);

        // 2. Candidates = P \ N(pivot)
        CliqueVec candidates;
        if (u_pivot != nullptr) {
            auto it = adj.find(u_pivot);
            if (it != adj.end()) {
                candidates = set_difference_helper(P, it->second);
            }
            else {
                candidates = P;
            }
        }
        else {
            candidates = P;
        }

        // 3. Recurse
        for (Node v : candidates) {
            CliqueVec newR = R;
            newR.push_back(v);

            CliqueVec neighbors_v;
            auto it = adj.find(v);
            if (it != adj.end()) neighbors_v = it->second;

            runBKPivot(
                newR,
                set_intersection_helper(P, neighbors_v),
                set_intersection_helper(X, neighbors_v),
                adj, hashMap
            );

            // Backtrack: Move v from P to X
            auto itP = std::lower_bound(P.begin(), P.end(), v);
            if (itP != P.end() && *itP == v) P.erase(itP);

            auto itX = std::lower_bound(X.begin(), X.end(), v);
            X.insert(itX, v);
        }
    }

    // --- ALGORITHM 2: BK RCD (Recursive Core Decomposition) ---
    // Được sử dụng cho các vùng "đặc" (dense neighborhoods)
    void runBKRcd(
        CliqueVec R,
        CliqueVec P,
        CliqueVec X,
        const AdjMap& adj,
        ResultMap& hashMap)
    {
        if (P.empty() && X.empty()) {
            report_clique(R, hashMap);
            return;
        }

        // Loop Decomposition: Tiếp tục loại bỏ đỉnh cho đến khi P là Clique
        while (true) {
            // Kiểm tra xem P có phải là Clique hay không, đồng thời tìm đỉnh có bậc thấp nhất trong P
            // Trong ngữ cảnh này: Bậc thấp nhất trong P <=> Có nhiều non-neighbor nhất trong P

            bool isClique = true;
            Node u_worst = nullptr;
            int min_degree_in_P = 2147483647; // INT_MAX

            // Duyệt qua tất cả đỉnh trong P để tính bậc nội bộ
            for (Node u : P) {
                // Tính bậc của u trong subgraph P (giao của N(u) và P)
                auto itAdj = adj.find(u);
                int deg_in_P = 0;
                if (itAdj != adj.end()) {
                    deg_in_P = count_intersection(P, itAdj->second);
                }

                // Nếu có bất kỳ đỉnh nào không nối với tất cả đỉnh còn lại (bậc < |P| - 1)
                // thì P chưa phải là Clique.
                if (deg_in_P < (int)P.size() - 1) {
                    isClique = false;
                }

                // Tìm đỉnh tệ nhất (bậc nhỏ nhất) để loại bỏ
                if (deg_in_P < min_degree_in_P) {
                    min_degree_in_P = deg_in_P;
                    u_worst = u;
                }
            }

            // CASE 1: P đã là Clique (Remaining Clique)
            if (isClique) {
                // Kiểm tra tính tối đại với X
                // Một clique P hợp R là tối đại nếu không có node x nào trong X nối với TẤT CẢ node trong P
                bool isMaximal = true;
                if (!P.empty()) {
                    for (Node x : X) {
                        bool connectedToAll = true;
                        auto itAdjX = adj.find(x);
                        const CliqueVec& neighbors_x = (itAdjX != adj.end()) ? itAdjX->second : CliqueVec();

                        // Check if neighbors_x contains all of P
                        // Vì cả 2 đều sorted, có thể check nhanh, nhưng ở đây dùng count_intersection cho đơn giản
                        // Nếu intersection(P, N(x)) == |P| -> x nối hết với P
                        if (count_intersection(P, neighbors_x) != (int)P.size()) {
                            connectedToAll = false;
                        }

                        if (connectedToAll) {
                            isMaximal = false;
                            break; // P bị chặn bởi x
                        }
                    }
                }

                if (isMaximal) {
                    // Output R U P
                    CliqueVec resultClique = R;
                    resultClique.insert(resultClique.end(), P.begin(), P.end());
                    report_clique(resultClique, hashMap);
                }
                return; // Kết thúc nhánh này
            }

            // CASE 2: P chưa là Clique -> Gọt vỏ (Decomposition)
            // Chọn u_worst (đỉnh có nhiều non-neighbor nhất)

            // a. Gọi đệ quy với u_worst được bao gồm
            CliqueVec newR = R;
            newR.push_back(u_worst);

            auto itAdj = adj.find(u_worst);
            const CliqueVec& neighbors_u = (itAdj != adj.end()) ? itAdj->second : CliqueVec();

            runBKRcd(
                newR,
                set_intersection_helper(P, neighbors_u),
                set_intersection_helper(X, neighbors_u),
                adj, hashMap
            );

            // b. Loại bỏ u_worst khỏi P và thêm vào X cho vòng lặp while tiếp theo
            // (Tương đương P = P \ {u}, X = X U {u})
            auto itP = std::lower_bound(P.begin(), P.end(), u_worst);
            if (itP != P.end() && *itP == u_worst) P.erase(itP);

            auto itX = std::lower_bound(X.begin(), X.end(), u_worst);
            X.insert(itX, u_worst);

            // Nếu P rỗng thì dừng
            if (P.empty()) return;
        }
    }

    // --- STRUCTURAL ANALYSIS (s, k-graph) ---
    // Tính s và k để quyết định dùng thuật toán nào
    struct StructureInfo {
        int s; // Kernel size (kết nối với tất cả)
        int k; // Shell size
    };

    StructureInfo analyzeStructure(const CliqueVec& P, const AdjMap& adj) {
        int n_sub = (int)P.size();
        if (n_sub == 0) return { 0, 0 };

        int s = 0;
        int k = 0;

        for (Node u : P) {
            // Tính bậc trong P
            auto it = adj.find(u);
            int deg_in_P = 0;
            if (it != adj.end()) {
                deg_in_P = count_intersection(P, it->second);
            }

            // Nếu nối với tất cả (trừ chính nó) -> thuộc S
            if (deg_in_P == n_sub - 1) {
                s++;
            }
            else {
                k++;
            }
        }
        return { s, k };
    }

    // --- DEGENERACY ORDERING ---
    // Tính thứ tự suy biến của đồ thị
    std::vector<Node> getDegeneracyOrdering(const AdjMap& adj) {
        // 1. Tính bậc ban đầu
        std::unordered_map<Node, int> degrees;
        // Bucket sort: max degree < N. Dùng vector<vector<Node>> làm bucket
        // Tuy nhiên vì Node là pointer, max degree có thể lớn bằng số instance.
        // Để đơn giản và an toàn bộ nhớ, ta dùng set<pair<degree, Node>> để luôn lấy min.
        // Độ phức tạp O(M log N), đủ tốt. (Linear O(N+M) bucket sort phức tạp hơn chút để implement)

        std::set<std::pair<int, Node>> sortedNodes;

        for (const auto& entry : adj) {
            int d = (int)entry.second.size();
            degrees[entry.first] = d;
            sortedNodes.insert({ d, entry.first });
        }

        std::vector<Node> ordering;
        ordering.reserve(adj.size());

        // Theo dõi các node đã bị xóa
        std::set<Node> removed;

        // 2. Core Decomposition
        while (!sortedNodes.empty()) {
            // Lấy đỉnh có bậc nhỏ nhất
            auto it = sortedNodes.begin();
            Node u = it->second;
            sortedNodes.erase(it);

            ordering.push_back(u);
            removed.insert(u);

            // Giảm bậc các lân cận
            auto itAdj = adj.find(u);
            if (itAdj != adj.end()) {
                for (Node v : itAdj->second) {
                    if (removed.find(v) != removed.end()) continue; // Đã xóa v rồi thì bỏ qua

                    // Cập nhật v trong sortedNodes
                    int old_deg = degrees[v];
                    auto searchPair = sortedNodes.find({ old_deg, v });
                    if (searchPair != sortedNodes.end()) {
                        sortedNodes.erase(searchPair);
                        degrees[v] = old_deg - 1;
                        sortedNodes.insert({ old_deg - 1, v });
                    }
                }
            }
        }
        return ordering;
    }
}

// ============================================================================
// PUBLIC METHODS IMPLEMENTATION
// ============================================================================

std::map<Colocation, std::unordered_map<FeatureType, std::set<const SpatialInstance*>>> MaximalCliqueHashmap::executeBK(
    const std::vector<NeighborSet>& neighborSets) {

    // --- Step 1: Build Adjacency Map Directly ---
    AdjMap adj;
    adj.reserve(neighborSets.size());

    // Thu thập tất cả các nodes
    for (const auto& ns : neighborSets) {
        Node u = ns.center;
        CliqueVec sorted_neighbors = ns.neighbors;
        std::sort(sorted_neighbors.begin(), sorted_neighbors.end());
        adj[u] = std::move(sorted_neighbors);
    }

    // --- Step 2: Compute Degeneracy Ordering ---
    std::vector<Node> ordering = getDegeneracyOrdering(adj);

    // --- Step 3: Iterate in Degeneracy Order ---
    // MCE Degeneracy Logic:
    // Với mỗi đỉnh v trong thứ tự suy biến:
    // P = N(v) giao {các đỉnh đứng SAU v trong thứ tự}
    // X = N(v) giao {các đỉnh đứng TRƯỚC v trong thứ tự}

    // Để tra cứu nhanh "đứng sau/trước", ta map Node -> index trong ordering
    std::unordered_map<Node, int> orderIndex;
    for (int i = 0; i < (int)ordering.size(); ++i) {
        orderIndex[ordering[i]] = i;
    }

    ResultMap hashMap;

    for (int i = 0; i < (int)ordering.size(); ++i) {
        Node v = ordering[i];

        // Lấy hàng xóm của v
        auto itAdj = adj.find(v);
        if (itAdj == adj.end()) continue;
        const CliqueVec& neighbors = itAdj->second;

        // Phân loại hàng xóm vào P (sau) và X (trước)
        CliqueVec P, X;
        P.reserve(neighbors.size());
        X.reserve(neighbors.size());

        for (Node neighbor : neighbors) {
            if (orderIndex[neighbor] > i) {
                P.push_back(neighbor);
            }
            else {
                X.push_back(neighbor);
            }
        }
        // P và X cần được sort để dùng cho set intersection trong các bước sau
        std::sort(P.begin(), P.end());
        std::sort(X.begin(), X.end());

        // --- HYBRID SWITCH ---
        // Phân tích cấu trúc của đồ thị con P
        StructureInfo info = analyzeStructure(P, adj);

        // Điều kiện chọn thuật toán (từ paper: s >= 2.8k - 11)
        // RCD tốt cho vùng đặc (s lớn, k nhỏ)
        // Pivot tốt cho vùng thưa
        double threshold = 2.8 * info.k - 11.0;

        // R khởi tạo chứa {v}
        CliqueVec R_init = { v };

        if (info.s >= threshold) {
            // Gọi BK RCD
            runBKRcd(R_init, P, X, adj, hashMap);
        }
        else {
            // Gọi BK Pivot
            runBKPivot(R_init, P, X, adj, hashMap);
        }
    }

    return hashMap;
}

std::priority_queue<Colocation, std::vector<Colocation>, ColocationPriorityComp> MaximalCliqueHashmap::extractInitialCandidates(
    const std::map<Colocation, std::unordered_map<FeatureType, std::set<const SpatialInstance*>>>& hashMap) {

    std::priority_queue<Colocation, std::vector<Colocation>, ColocationPriorityComp> candidateQueue;
    for (const auto& entry : hashMap) {
        const Colocation& maximalClique = entry.first;
        candidateQueue.push(maximalClique);
    }
    return candidateQueue;
}