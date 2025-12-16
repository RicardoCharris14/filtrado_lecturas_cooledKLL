/** uhr: generic time performance tester
 * Author: LELE
 *
 * Things to set up:
 * 0. Includes: include all files to be tested,
 * 1. Time unit: in elapsed_time,
 * 2. What to write on time_data,
 * 3. Data type and distribution of RNG,
 * 4. Additive or multiplicative stepping,
 * 5. The experiments: in outer for loop. */

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// Include to be tested files here
#include "../include/procesarKmers.hpp"
#include "cooled-kll.cpp"

inline void validate_input(int argc, char *argv[], std::int64_t& runs, int& method, int& k)
{
    if (argc != 5) {
        std::cerr << "Usage: <filename> <RUNS> <METHOD> <k>" << std::endl;
        std::cerr << "<filename> is the name of the file where performance data will be written." << std::endl;
        std::cerr << "It is recommended for <filename> to have .csv extension and it should not previously exist." << std::endl;
        std::cerr << "<RUNS>: numbers of runs per test case: should be >= 32." << std::endl;
        std::cerr << "<METHOD>: 1 = plain vector | 2 = compressed vector | 3 = sketch" << std::endl;
        std::cerr << "<k>: length of kmers." << std::endl;
        std::cerr << "These should all be positive." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Read command line arguments
    try {
        runs = std::stoll(argv[2]);
        method = std::stoi(argv[3]);
        k = std::stoi(argv[4]);
    } catch (std::invalid_argument const& ex) {
        std::cerr << "std::invalid_argument::what(): " << ex.what() << std::endl;
        std::exit(EXIT_FAILURE);
    } catch (std::out_of_range const& ex) {
        std::cerr << "std::out_of_range::what(): " << ex.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Validate arguments
    if (runs < 4) {
        std::cerr << "<RUNS> must be at least 4." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (method != 1 and method != 2 and method != 3){
        std::cerr << "<METHOD>: must be 1, 2 or 3." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (k < 1) {
        std::cerr << "<k> must be at least 1." << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

inline void display_progress(std::int64_t u, std::int64_t v)
{
    const double progress = u / double(v);
    const std::int64_t width = 70;
    const std::int64_t p = width * progress;
    std::int64_t i;

    std::cout << "\033[1m[";
    for (i = 0; i < width; i++) {
        if (i < p)
            std::cout << "=";
        else if (i == p)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << std::int64_t(progress * 100.0) << "%\r\033[0m";
    std::cout.flush();
}

int main(int argc, char *argv[])
{
    // Validate and sanitize input
    std::int64_t runs, lower = 1, upper = 1000, step = 1;
    int method, k;
    validate_input(argc, argv, runs, method, k);

    // Set up clock variables
    std::int64_t n, i, executed_runs;
    std::int64_t total_runs_additive = runs * (((upper - lower) / step) + 1);
    std::int64_t total_runs_multiplicative = runs * (floor(log(upper / double(lower)) / log(step)) + 1);
    std::vector<double> times(runs);
    std::vector<double> times2(runs);
    std::vector<double> q;
    double mean_time, time_stdev, dev;
    double mean_time2, time_stdev2, dev2;
    auto begin_time = std::chrono::high_resolution_clock::now();
    auto end_time = std::chrono::high_resolution_clock::now();
    auto begin_time2 = std::chrono::high_resolution_clock::now();
    auto end_time2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::nano> elapsed_time = end_time - begin_time;
    std::chrono::duration<double, std::nano> elapsed_time2 = end_time - begin_time;

    // Set up random number generation
    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::uniform_int_distribution<std::int64_t> u_distr; // change depending on app

    // File to write time data
    std::ofstream time_data;
    time_data.open(argv[1]);
    time_data << "quantile,quantile_t_mean,quantile_t_stdev,rank,rank_t_mean,rank_t_stdev" << std::endl;

    // Begin testing
    float quantile = 0.001;
    float increment = 0.001;

    std::vector<std::pair<uint64_t, size_t>> kmers = procesarKMers("Genomas", k);

    std::vector<std::pair<uint64_t, size_t>> compressed_vector;

    size_t n_buckets = 100, buckets_capacity = 10;
    int compactor_size = 100, eviction_threshold = 16;
    float compression_factor = 0.7;
    size_t total_kmers = kmers.size();
    CooledKLL sketch(n_buckets, buckets_capacity, eviction_threshold, compactor_size, compression_factor);

    size_t quantile_val;

    switch(method){
        // plain vector
        case 1:
        {
            std::sort(kmers.begin(), kmers.end(), [](const auto& a, const auto& b){
                return a.second < b.second;
            });
            break;
        }

        // compressed vector
        case 2:
        {
            std::sort(kmers.begin(), kmers.end(), [](const auto& a, const auto& b){
                return a.second < b.second;
            });
            size_t total_kmers = kmers.size();
            uint64_t prev_item = kmers[0].second, frequency = 0;
            for (size_t i=0 ; i<total_kmers ; i++){
                if (kmers[i].second == prev_item){
                    frequency += 1;
                } else {
                    compressed_vector.push_back({prev_item, frequency});
                    prev_item = kmers[i].second;
                    frequency = 1;
                }

                if (i == total_kmers-1){
                    compressed_vector.push_back({prev_item, frequency});
                }
            }
            break;
        }

        // sketch
        case 3:
        {
            for (size_t i=0 ; i<total_kmers ; i++){
                sketch.insert(kmers[i].second);
            }
            break;
        }

        default:
            break;
    }

    std::cout << "\033[0;36mRunning tests...\033[0m" << std::endl << std::endl;
    executed_runs = 0;
    for (n = lower; n <= upper; n += step) {
        mean_time = 0;
        time_stdev = 0;

        // Test configuration goes here
        
        
        // Run to compute elapsed time
        for (i = 0; i < runs; i++) {
            // Remember to change total depending on step type
            display_progress(++executed_runs, total_runs_additive);

            
            // Function to test goes here

            switch(method){
                // plain vector
                case 1:
                {
                    // quantile experiment
                    begin_time = std::chrono::high_resolution_clock::now();

                    size_t idx = static_cast<size_t>(std::ceil(total_kmers * quantile));
                    if (idx >= total_kmers) idx = total_kmers - 1;
                    quantile_val = kmers[idx].second;

                    end_time = std::chrono::high_resolution_clock::now();

                    // rank experiment
                    begin_time2 = std::chrono::high_resolution_clock::now();

                    size_t rank = 0;
                    bool found = false;
                    for (size_t i=0 ; i<total_kmers ; i++){
                        if (found and kmers[i].second != quantile_val){
                            break;
                        }
                        if (not found and kmers[i].second == quantile_val) found = true;

                        rank += 1;
                    }

                    end_time2 = std::chrono::high_resolution_clock::now();
                    break;
                }

                // compressed vector
                case 2:
                {
                    // quantile experiment
                    begin_time = std::chrono::high_resolution_clock::now();

                    size_t idx = static_cast<size_t>(std::ceil(total_kmers * quantile));
                    size_t cumulative_frequency = 0;
                    size_t different_kmers = compressed_vector.size();
                    for (size_t i=0 ; i<different_kmers ; i++){
                        cumulative_frequency += compressed_vector[i].second;
                        if (cumulative_frequency-1 >= idx){
                            quantile_val = compressed_vector[i].first;
                            break;
                        }
                        if (i == different_kmers-1) quantile_val = compressed_vector[i].first;
                    }

                    end_time = std::chrono::high_resolution_clock::now();

                    // rank experiment
                    begin_time2 = std::chrono::high_resolution_clock::now();

                    size_t rank = 0;
                    for (size_t i=0 ; i<different_kmers ; i++){
                        rank += compressed_vector[i].second;
                        if (compressed_vector[i].first == quantile_val){
                            break;
                        }
                    }

                    end_time2 = std::chrono::high_resolution_clock::now();
                    break;
                }

                // Sketch
                case 3:
                {
                    // quantile experiment
                    begin_time = std::chrono::high_resolution_clock::now();

                    quantile_val = sketch.quantile(quantile);

                    end_time = std::chrono::high_resolution_clock::now();

                    // rank experiment
                    begin_time2 = std::chrono::high_resolution_clock::now();

                    size_t rank = sketch.rank(quantile_val);

                    end_time2 = std::chrono::high_resolution_clock::now();
                    break;
                }
            }


            elapsed_time = end_time - begin_time;
            elapsed_time2 = end_time2 - begin_time2;
            times[i] = elapsed_time.count();
            times2[i] = elapsed_time2.count();

            mean_time += times[i];
            mean_time2 += times2[i];
        }

        quantile += increment;

        // Compute statistics
        mean_time /= runs;
        mean_time2 /= runs;

        for (i = 0; i < runs; i++) {
            dev = times[i] - mean_time;
            time_stdev += dev * dev;

            dev2 = times2[i] - mean_time2;
            time_stdev2 += dev2 * dev2;
        }

        time_stdev /= runs - 1; // Subtract 1 to get unbiased estimator
        time_stdev = std::sqrt(time_stdev);

        time_stdev2 /= runs - 1;
        time_stdev2 = std::sqrt(time_stdev2);


        time_data << quantile << "," << mean_time << "," << time_stdev << "," << quantile_val << "," << mean_time2 << "," << time_stdev2 << std::endl;
    }

    // This is to keep loading bar after testing
    std::cout << std::endl << std::endl;
    std::cout << "\033[1;32mDone!\033[0m" << std::endl;

    time_data.close();

    return 0;
}