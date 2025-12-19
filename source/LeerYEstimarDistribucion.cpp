#include <iostream>
#include <string>
#include <vector>
#include "../include/procesarKmers.hpp"
#include "../include/experiments.hpp"


int main(int argc, char* argv[]){
    // Verificacion de correctitud en la ejecucion del programa
    if (argc != 7){
        std::cerr << "correct usage: ./exe <folder_file> <k-mers_length> <distribution> <N_buckets> <B_capacity> <C_size>" << std::endl;
        std::cerr << "<folder_file>: path to the folder with genomic data." << std::endl;
        std::cerr << "<k-mers_length>: length of kmers." << std::endl;
        std::cerr << "<distribution>: define type of distribution: 0 = kmers distribution | 1 = frequency distribution." << std::endl;
        std::cerr << "<N_buckets>: number of buckets in the hot filter part of the sketch." << std::endl;
        std::cerr << "<B_capacity>: number of entries a bucket have." << std::endl;
        std::cerr << "<C_size>: number of elements of the largest compactor in the classic kll part." << std::endl;
        return 1;
    }

    std::string folder_path;
    int k_;
    bool frequency_distribution;
    double quantile_ratio = 0.001;
    // KLL settings
    size_t n_buckets, buckets_capacity; 
    int compactor_size; 
    float compression_factor = 0.7;

    // Verificacion de pertinencia de los argumentos
    try{
        folder_path = argv[1];
        k_ = std::stoi(argv[2]);
        frequency_distribution = std::stoi(argv[3]);
        n_buckets = std::stoll(argv[4]);
        buckets_capacity = std::stoll(argv[5]);
        compactor_size = std::stoll(argv[6]); 

        if (k_ <= 0 or k_ > 31){
            std::cerr << "<k-mer length> must be a number belonging to [1, 31]." << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (frequency_distribution != 0 and frequency_distribution != 1){
            std::cerr << "<distribution> must be a 0 or 1." << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (n_buckets <= 0){
            std::cerr << "<N_buckets> must be greater than 0." << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (buckets_capacity <= 0){
            std::cerr << "<B_capacity> must be greater than 0." << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (compactor_size <= 0){
            std::cerr << "<C_size> must be greater than 0." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    } catch (std::exception e){
        std::cerr << "Error: " << e.what() << std::endl;
        std::exit(1);
    }

    std::cout << "!Leyendo kmers!" << std::endl;

    std::vector<std::pair<uint64_t, size_t>> kmers_distribution = procesarKMers(folder_path, k_);

    std::sort(kmers_distribution.begin(), kmers_distribution.end(), [](const auto& a, const auto &b){
        return a.second < b.second;
    });

    if (frequency_distribution){
        frequencyExperiments(kmers_distribution, k_, quantile_ratio, n_buckets, buckets_capacity, compactor_size, compression_factor);
    } else {
        kmersExperiments(kmers_distribution, k_, quantile_ratio, n_buckets, buckets_capacity, compactor_size, compression_factor);
    }
    
}