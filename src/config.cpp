/**
 * @file config.cpp
 * @brief Implementation of configuration file loading
 */

#include "config.h"
#include <iostream>
#include <fstream>
#include <sstream>


 /**
  * @brief Load configuration from a file
  * @param configPath Path to the configuration file
  * @return AppConfig Configuration object with loaded or default values
  *
  * Parses key=value pairs from the configuration file. Lines starting with '#'
  * are treated as comments. If the file cannot be opened, returns default configuration.
  */
AppConfig ConfigLoader::load(const std::string& configPath) {
    AppConfig config;
    std::ifstream file(configPath);

    if (!file.is_open()) {
        std::cerr << "Warning: Config file not found, using defaults.\n";
        return config;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream is_line(line);
        std::string key;

        if (std::getline(is_line, key, '=')) {
            std::string value;
            if (std::getline(is_line, value)) {
                if (key == "dataset_path") config.datasetPath = value;
                else if (key == "neighbor_distance") config.neighborDistance = std::stod(value);
                else if (key == "min_prevalence") config.minPrev = std::stod(value);
                else if (key == "min_cond_prob") config.minCondProb = std::stod(value);
                else if (key == "percentage_instances") config.percentageData = std::stod(value);
                else if (key == "debug_mode") config.debugMode = (value == "true" || value == "1");
            }
        }
    }

    return config;
}