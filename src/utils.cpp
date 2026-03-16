/**
 * @file utils.cpp
 * @brief Implementation: Feature counting, dispersion calculation, and rare intensity
 */

#include "utils.h"
#include <unordered_map>
#include <set>
#include <chrono>
#include <windows.h>
#include <psapi.h>
#include <iostream> 
#include <iomanip>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <vector>
#include <numeric>

// Count instances per feature type and sort by frequency (ascending)
std::map<FeatureType, int> countFeatures(
	const std::vector<SpatialInstance>& instances) {
	//////// TODO: Implement (1)//////////
	std::map<FeatureType, int> counts;
	for (const auto& instance : instances) {
		counts[instance.type]++;
	}
	return counts;
};

// Calculate dispersion (delta) from feature distribution
double calculateDispersion(const std::map<FeatureType, int>& featureCount) {
    //////// TODO: Implement (2)//////////
    if (featureCount.empty()) return 0.0;
    size_t m = featureCount.size();
    if (m <= 1) return 0.0;

    // --- Bước 1: Trích xuất và Sắp xếp (Sort) ---
    // Chuyển tần suất từ Map sang Vector để sắp xếp
    std::vector<double> frequencies;
    frequencies.reserve(m);

    for (const auto& pair : featureCount) {
        frequencies.push_back(static_cast<double>(pair.second));
    }

    // Sắp xếp tăng dần (Ascending) để đảm bảo num(fi) <= num(fj) khi i < j
    std::sort(frequencies.begin(), frequencies.end());

    // --- Bước 2: Tính tổng tỷ lệ ---
    // Công thức: sum_{i<j} (num(fj) / num(fi))
    double sumRatio = 0.0;

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = i + 1; j < m; ++j) {
            if (frequencies[i] > 0) {
                sumRatio += frequencies[j] / frequencies[i];
            }
        }
    }

    // --- Bước 3: Áp dụng hệ số ---
    // Hệ số: 2 / (m * (m - 1))
    double factor = 2.0 / (static_cast<double>(m) * (static_cast<double>(m) - 1.0));

    double delta = factor * sumRatio;

    return delta;
};


// Calculate rare intensity for each instance in a colocation
std::unordered_map<FeatureType, double> calcRareIntensity(
    Colocation c,
    const std::map<FeatureType, int>& featureCounts,
    double delta) {
        //////// TODO: Implement (11)//////////
    std::unordered_map<FeatureType, double> intensityMap;
    if (c.empty()) return intensityMap;

    // 1. Find N(f_min)
    // Initialize with first feature's count
    int minCount = -1;
    
    for (const auto& f : c) {
        if (featureCounts.find(f) != featureCounts.end()) {
            int count = featureCounts.at(f);
            if (minCount == -1 || count < minCount) {
                minCount = count;
            }
        }
    }

    if (minCount <= 0) return intensityMap; 

    // Tránh lỗi chia cho 0 nếu delta = 0
    if (delta == 0) delta = 1e-9; 

    for (const auto& f : c) {
        if (featureCounts.find(f) != featureCounts.end()) {
            int count = featureCounts.at(f);
            if (count > 0) {
                // B. Tính v
                double v = static_cast<double>(count) / minCount;

                // C. Tính trọng số theo CAUCHY KERNEL
                // w = 1 + ((v - 1) / delta)^2
                double w = 1.0 + std::pow((v - 1.0) / delta, 2);

                // Rare Intensity (RI) là nghịch đảo của trọng số (RI = 1 / w)
                // Để khi sang miner.cpp, w_log = 1 / RI sẽ trả lại đúng giá trị w này.
                double ri = 1.0 / w;

                // Store result (Mapping FeatureType -> RI)
                intensityMap[f] = ri;
            }
        }
    }

    return intensityMap;
};