/**
 * @file miner.cpp
 * @brief Implementation: Mining prevalent colocation patterns using weighted PI
 */

#include "miner.h"
#include "utils.h"
#include <vector>
#include <queue>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <iostream>


// Main mining algorithm: find all prevalent colocation patterns
std::set<Colocation> Miner::minePCPs(
	std::priority_queue<Colocation, std::vector<Colocation>, ColocationPriorityComp>& candidateColocations,
	const std::map<Colocation, std::unordered_map<FeatureType, std::set<const SpatialInstance*>>>& hashMap,
	const std::map<FeatureType, int>& featureCounts,
	double delta,
	double min_prev) {

	std::set<Colocation> prevalentPCs;
	std::set<Colocation> nonPrevalentPCs;
	std::set<Colocation> visited;

	while (!candidateColocations.empty()) {
		Colocation c = candidateColocations.top();
		candidateColocations.pop();

		if (visited.count(c)) continue;
		visited.insert(c);

		std::set<Colocation> newCs;

		auto partInstances = queryInstances(c, hashMap);
		auto rareIntensityMap = calcRareIntensity(c, featureCounts, delta);

		double weightedPI = computeWeightedPI(partInstances, c, rareIntensityMap, featureCounts);
		newCs = generateSubsets(c);

		if (weightedPI >= min_prev) {
			prevalentPCs.insert(c);

			auto prevalentSubsets = deducePrevalentSubsets(newCs, c, featureCounts);
			for (const auto& subset : prevalentSubsets) {
				prevalentPCs.insert(subset);
			}

			std::set<Colocation> filteredSubsets;
			for (const auto& subset : newCs) {
				if (!prevalentSubsets.count(subset)) {
					filteredSubsets.insert(subset);
				}
			}
			newCs = filteredSubsets;
		}
		else {
			nonPrevalentPCs.insert(c);
		}

		for (const auto& subset : newCs) {
			if (!visited.count(subset)) {
				candidateColocations.push(subset);
			}
		}
	}

	return prevalentPCs;
}


// Query instances of a colocation from hashmap
std::map<FeatureType, std::set<const SpatialInstance*>> Miner::queryInstances(
	Colocation c,
	const std::map<Colocation, std::unordered_map<FeatureType, std::set<const SpatialInstance*>>>& hashMap) {
		//////// TODO: Implement (10)/////////

	std::map<FeatureType, std::set<const SpatialInstance*>> instancesMap;

	for (const auto& entry : hashMap) {
		const Colocation& maximalClique = entry.first;
		const auto& cliqueInstances = entry.second;

		// Check if c is a subset of maximalClique
		bool isSubset = std::includes(
			maximalClique.begin(), maximalClique.end(),
			c.begin(), c.end()
		);

		// If c is a subset, merge instances
		if (isSubset) {
			for (const auto& f : c) {
				if (cliqueInstances.count(f)) {
					const auto& insts = cliqueInstances.at(f);
					instancesMap[f].insert(insts.begin(), insts.end());
				}
			}
		}
	}

	return instancesMap;
};

// Compute weighted participation index for a colocation
double Miner::computeWeightedPI(
	const std::map<FeatureType, std::set<const SpatialInstance*>>& partInstances,
	Colocation c,
	const std::unordered_map<FeatureType, double>& rareIntensityMap,
	const std::map<FeatureType, int>& featureCounts) {
		//////// TODO: Implement (12)/////////
	if (c.empty()) return 0.0;
	
	double minWPR = -1.0;

	for (const auto& f : c) {
		// Calculate PR = count / N
		int count = 0;
		if (partInstances.count(f)) {
			count = partInstances.at(f).size();
		}

		int totalCount = 0;
		if (featureCounts.count(f)) {
			totalCount = featureCounts.at(f);
		}

		if (totalCount == 0) continue;

		double pr = static_cast<double>(count) / totalCount;

		// Retrieve RI
		double ri = 1.0; // Default if missing? Or should skip? 
		// If RI is missing, it implies error or skip. 
		// According to logic, RI depends on featureCounts. If >0 count, RI exists.
		if (rareIntensityMap.count(f)) {
			ri = rareIntensityMap.at(f);
		}

		// Calculate W_log = 1 / RI
		if (ri < 1e-9) {
			std::cout<<"Warning!!";
			ri = 1e-9;
		}
		double w_log = 1.0 / ri;

		// WPR = PR * W_log
		double wpr = pr * w_log;

		if (minWPR < 0 || wpr < minWPR) {
			minWPR = wpr;
		}
	}

	return (minWPR < 0) ? 0.0 : minWPR;
};

// Generate all size-1 subsets (remove one feature at a time)
std::set<Colocation> Miner::generateSubsets(const Colocation& c) {
	std::set<Colocation> subsets;
	if (c.size() <= 2) return subsets;

	for (size_t i = 0; i < c.size(); ++i) {
		Colocation sub = c;
		sub.erase(sub.begin() + i);
		subsets.insert(sub);
	}
	return subsets;
};

// Deduce prevalent subsets using downward closure property
std::set<Colocation> Miner::deducePrevalentSubsets(
	std::set<Colocation>& subsets,
	const Colocation& c,
	const std::map<FeatureType, int>& featureCounts) {
	
	std::set<Colocation> provenPrevalent;
	
	// 1. Find f_min (feature with min frequency) in parent c
	FeatureType f_min;
	int minCount = -1;

	for (const auto& f : c) {
		int count = 0;
		if (featureCounts.count(f)) {
			count = featureCounts.at(f);
		}
		
		if (minCount == -1 || count < minCount) {
			minCount = count;
			f_min = f;
		} else if (count == minCount) {
			// Optional tie-breaker
			if (f < f_min) f_min = f;
		}
	}

	// 2. Lemma 2: If C is prevalent, any subset C' containing f_min is also prevalent.
	if (!f_min.empty()) {
		for (const auto& subset : subsets) {
			bool hasMinFeature = false;
			for (const auto& f : subset) {
				if (f == f_min) {
					hasMinFeature = true;
					break;
				}
			}

			if (hasMinFeature) {
				provenPrevalent.insert(subset);
			}
		}
	}
	return provenPrevalent;
};