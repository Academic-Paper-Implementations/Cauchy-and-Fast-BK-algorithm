/**
 * @file main.cpp
 * @brief Simplified Entry point for Co-location Mining
 */

#include "config.h"
#include "data_loader.h"
#include "neighbor_graph.h"
#include "maximal_clique_hashmap.h"
#include "miner.h"
#include "types.h"
#include "utils.h"
#include <iostream>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <fstream>

//Show memmory usage
#include <windows.h>
#include <psapi.h>
#include <stdio.h>
#pragma comment(lib, "psapi.lib")

int main(int argc, char* argv[]) {
    auto programStart = std::chrono::high_resolution_clock::now();
    // --- Step 1: Config & Load Data ---
    std::cout << "Running... (Results will be saved to result.txt)\n";
    std::string config_path = (argc > 1) ? argv[1] : "./config/config.txt";
    AppConfig config = ConfigLoader::load(config_path);

    auto instances = DataLoader::load_csv(config.datasetPath, config.percentageData);

    // --- Step 2: Pre-processing (Indexing & Structures) ---
    // 1. Feature Counting & Sorting
    auto featureCount = countFeatures(instances);

	// 2. Delta Calculation
	double delta = calculateDispersion(featureCount);

	// 3. Neighbor Graph Building
    NeighborGraph neighborGraph;
    auto graph = neighborGraph.buildNeighborGraph(instances, config.neighborDistance);

	// 4. Build Instance Hashmap from Maximal Cliques
	MaximalCliqueHashmap mcHashmap;
    auto hashMap = mcHashmap.executeBK(graph);

	// 5. Get Candidate Colocations
	auto candidateQueue = mcHashmap.extractInitialCandidates(hashMap);

    // --- Step 3: Mining Prevalent Co-location Patterns ---
    Miner miner;
    auto colocations = miner.minePCPs(
        candidateQueue,
        hashMap,
        featureCount,
        delta,
        config.minPrev
    );

    // --- END OF PROCESSING ---
    auto programEnd = std::chrono::high_resolution_clock::now();
    double totalTime = std::chrono::duration<double>(programEnd - programStart).count();

    // --- REPORT GENERATION (FILE ONLY) ---
    // 1. Get Memory Info (Peak)
    HANDLE handle = GetCurrentProcess();
    PROCESS_MEMORY_COUNTERS memCounter;
    SIZE_T peakMemMB = 0;

    if (GetProcessMemoryInfo(handle, &memCounter, sizeof(memCounter))) {
        peakMemMB = memCounter.PeakWorkingSetSize / 1024 / 1024; // Convert to MB
    }

    // 2. Write to File
    std::ofstream outFile("../results.txt");
    if (!outFile.is_open()) {
        std::cerr << "Cannot open results.txt for writing.\n";
        return 1;
    }
    // (A) Thông tin Dataset & Config
    outFile << "=== FINAL REPORT ===\n";
    outFile << "Dataset Path:      " << config.datasetPath << "\n";
    outFile << "Total Instances:   " << instances.size() << "\n";
    outFile << "Neighbor Distance: " << config.neighborDistance << "\n";
    outFile << "Min Prevalence:    " << config.minPrev << "\n";
    outFile << "Percentage Data:    " << (config.percentageData * 100) << "%\n";
    outFile << "----------------------------------------\n";

	// (B) Execution Time
    outFile << "Execution Time: " << std::fixed << std::setprecision(3) << totalTime << " s\n";

    // (C) Peak Memory Usage
    outFile << "Peak Memory Usage: " << peakMemMB << " MB\n";

	// (D) Number of Patterns Found
    outFile << "Patterns Found: " << colocations.size() << "\n";
    outFile << "----------------------------------------\n";

	// (E) List of Patterns
    if (!colocations.empty()) {
        int idx = 1;
        for (const auto& col : colocations) {
            outFile << "[" << idx++ << "] {";
            for (size_t i = 0; i < col.size(); ++i) {
                outFile << (i > 0 ? ", " : "") << col[i];
            }
            outFile << "}\n";
        }
    }
    else {
        outFile << "No patterns found.\n";
    }

    outFile.close();

    std::cout << "Done! Please check 'result.txt'.\n";
    return 0;
}