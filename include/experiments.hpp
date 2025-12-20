#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H

#include <iostream>
#include <vector>
#include <utility>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include "../source/cooled-kll.cpp"

void frequencyExperiments(std::vector<std::pair<uint64_t, uint64_t>>& kmers_dist, int k_, float quantile_ratio,
     size_t n_buckets=100, size_t buckets_capacity = 10, int compactor_size = 100, float compression_factor = 0.7){
    std::cout << "!ESTIMACION DE DISTRIBUCION DE FRECUENCIAS!" << std::endl;
    std::cout << "!Ordenando la distribucion de frecuencias!" << std::endl;

    // Ordena los datos en funcion de la frecuencia de los kmers
    std::sort(kmers_dist.begin(), kmers_dist.end(), [](const auto& a, const auto& b){
        return a.second < b.second;
    });

    std::cout << "!Calculando quantiles y ranks reales!" << std::endl;

    // Obtiene los cuantiles y ranks reales
    size_t total_kmers = kmers_dist.size(), idx, last_item_rank = total_kmers;
    std::vector<uint64_t> real_quantiles;
    std::vector<size_t> real_ranks;
    uint64_t last_item = UINT64_MAX;
    for (double i=0.0 ; i<=1.00000 ; i+=quantile_ratio){
        idx = static_cast<size_t>(std::ceil(total_kmers * i));
        if (idx >= total_kmers) idx = total_kmers - 1;
        real_quantiles.push_back(kmers_dist[idx].second);
        int j = 1;

        if (kmers_dist[idx].second == last_item){
            real_ranks.push_back(last_item_rank);
            continue;
        }

        while (true){
            last_item = kmers_dist[idx].second;
            if (idx + j == total_kmers){
                real_ranks.push_back(total_kmers);
                last_item_rank = total_kmers;
                break;
            }
            if (kmers_dist[idx].second != kmers_dist[idx + j].second){
                real_ranks.push_back(idx + j);
                last_item_rank = idx + j;
                break;
            }
            j++;
        }
    }

    std::cout << "!Calculando memoria real usada!" << std::endl;

    // Se calcula memoria del vector de frecuencias
    size_t vector_memory = 0;
    vector_memory += sizeof(std::vector<uint64_t>);
    vector_memory += kmers_dist.size() * sizeof(uint64_t);

    // Se calcula memoria usada por vector en formato comprimido
    std::vector<std::pair<uint64_t, uint64_t>> compressed_vector;
    size_t cv_memory = 0;
    uint64_t prev_item = kmers_dist[0].second, frequency = 0;
    for (size_t i=0 ; i<total_kmers ; i++){
        if (kmers_dist[i].second == prev_item){
            frequency += 1;
        } else {
            compressed_vector.push_back({prev_item, frequency});
            prev_item = kmers_dist[i].second;
            frequency = 1;
        }

        if (i == total_kmers-1){
            compressed_vector.push_back({prev_item, frequency});
        }
    }

    cv_memory += sizeof(compressed_vector);
    cv_memory += compressed_vector.size() * sizeof(std::pair<uint64_t, uint64_t>);
    cv_memory += compressed_vector.size() * 2 * sizeof(uint64_t);

    size_t unique_elements = compressed_vector.size();

    std::vector<std::pair<uint64_t, uint64_t>>().swap(compressed_vector);

    // Inserta en el sketch todas las frecuencias de los kmers
    int eviction_threshold = 16;

    CooledKLL sketch(n_buckets, buckets_capacity, eviction_threshold, compactor_size, compression_factor);

    std::cout << "Insertando datos en el sketch" << std::endl;

    for (size_t i=0 ; i<kmers_dist.size() ; i++){
        sketch.insert(kmers_dist[i].second);
    }

    std::cout << "!Estimando y guardando distribución de los datos!" << std::endl;

    std::filesystem::path folder_route = "data/frequency_distribution/NB_"+
                    std::to_string(n_buckets)+"_BC_"+std::to_string(buckets_capacity)+"_CS_"
                    +std::to_string(compactor_size);

    try{
        if (not std::filesystem::exists(folder_route)){
            std::filesystem::create_directories(folder_route);
        }
    } catch (const std::filesystem::filesystem_error& e){
        std::cerr << "Error creating folder: " << e.what() << std::endl;
    }

    std::string csvFilename = folder_route.string()+"/"+ std::to_string(k_) + "mers_distribution.csv";
    std::ofstream csvFile(csvFilename);
    if (csvFile.is_open()){
        csvFile << "quantile,real_quantile,estimated_quantile,rank,real_rank,estimated_rank\n";

        // Obtiene los cuantiles estimados y guarda los datos en un csv
        size_t estimated_quantile, estimated_rank, j = 0;
        for (double i=0.0 ; i <= 1.00000 ; i+=quantile_ratio){
            estimated_quantile = sketch.quantile(i);
            estimated_rank = sketch.rank(kmers_dist[real_ranks[j] - 1].second);

            csvFile << i << "," << real_quantiles[j] << "," << estimated_quantile << "," << 
            kmers_dist[real_ranks[j] - 1].second << "," << real_ranks[j] << "," << estimated_rank << std::endl;

            j++;
        }

        csvFile.close();
        std::cout << "!Datos guardados exitosamente!" << std::endl;
    }else{
        std::cerr << "Error creando CSV: " << csvFilename << std::endl;
    }

    std::cout << "!Almacenando memoria usada!" << std::endl;

    csvFilename = folder_route.string()+"/"+std::to_string(k_) + "mers_memory.csv";
    std::ofstream csvFile1(csvFilename);
    if (csvFile1.is_open()){
        csvFile1 << "elements,unique_elements,sketch_memory,vector_memory,compressed_vector_memory,n_buckets,b_capacity,comp_size,comp_factor\n";
        csvFile1 << total_kmers << "," << unique_elements << "," << sketch.memory() << "," << vector_memory << "," << cv_memory
         << "," << n_buckets << "," << buckets_capacity << "," << compactor_size << "," << compression_factor << std::endl;

        csvFile1.close();
    }else{
        std::cerr << "Error creando CSV: " << csvFilename << std::endl;
    }
}

void kmersExperiments(std::vector<std::pair<uint64_t, uint64_t>> kmers_dist, int k_, float quantile_ratio,
     size_t n_buckets=100, size_t buckets_capacity = 10, int compactor_size = 100, float compression_factor = 0.7){
    std::cout << "!ESTIMACION DE DISTRIBUCION DE KMERS!" << std::endl;
    std::cout << "!Calculando distribucion real de los datos!" << std::endl;

    // se define la distribucion de los kmers
    size_t different_kmers = kmers_dist.size();
    std::vector<uint64_t> kmers;
    for (size_t i=0 ; i<different_kmers ; i++){
        for (uint64_t j=0 ; j<kmers_dist[i].second ; j++){
            kmers.push_back(kmers_dist[i].first);
        }
    }

    std::cout << "!Calculando quantiles y ranks reales!" << std::endl;

    // Obtiene los cuantiles reales
    std::vector<uint64_t> real_quantiles;
    std::vector<size_t> real_ranks;
    size_t quantile_idx, total_kmers = kmers.size();
    uint64_t last_item = UINT64_MAX, last_item_rank = total_kmers;
    for (double i=0.0 ; i<=1.00000 ; i+=quantile_ratio){
        quantile_idx = static_cast<size_t>(std::ceil(total_kmers * i));
        if (quantile_idx == total_kmers){
            quantile_idx = total_kmers - 1;
        }
        real_quantiles.push_back(kmers[quantile_idx]);
        int k=1;

        if (kmers[quantile_idx] == last_item){
            real_ranks.push_back(last_item_rank);
            continue;
        }

        while (true){
            last_item = kmers[quantile_idx];
            if (quantile_idx + k == total_kmers){
                real_ranks.push_back(total_kmers);
                last_item_rank = total_kmers;
                break;
            }
            if (kmers[quantile_idx] != kmers[quantile_idx + k]){
                real_ranks.push_back(quantile_idx + k);
                last_item_rank = quantile_idx + k;
                break;
            }
            k++;
        }
    }

    

    std::cout << "!Calculando memoria real usada!" << std::endl;

    // se calcula memoria usada por vector de kmers
    size_t vector_memory = 0;
    vector_memory += sizeof(kmers);
    vector_memory += total_kmers * sizeof(uint64_t);

    // Se calcula memoria usada por vector en formato comprimido
    std::vector<std::pair<uint64_t, uint64_t>> compressed_vector;
    size_t cv_memory = 0;
    uint64_t prev_item = kmers[0], frequency = 0;
    for (size_t i=0 ; i<total_kmers ; i++){
        if (kmers[i] == prev_item){
            frequency += 1;
        } else {
            compressed_vector.push_back({prev_item, frequency});
            prev_item = kmers[i];
            frequency = 1;
        }

        if (i == total_kmers-1){
            compressed_vector.push_back({prev_item, frequency});
        }
    }

    cv_memory += sizeof(compressed_vector);
    cv_memory += compressed_vector.size() * sizeof(std::pair<uint64_t, uint64_t>);
    cv_memory += compressed_vector.size() * 2 * sizeof(uint64_t);

    size_t unique_elements = compressed_vector.size();

    std::vector<std::pair<uint64_t, uint64_t>>().swap(compressed_vector);

    std::cout << "!Insertando datos en el sketch!" << std::endl;

    // Inserta en el sketch todos los kmers
    int eviction_threshold = 16;
    CooledKLL sketch(n_buckets, buckets_capacity, eviction_threshold, compactor_size, compression_factor);
    for (size_t i=0 ; i<different_kmers ; i++){
        sketch.insert(kmers_dist[i].first, kmers_dist[i].second);
    }

    std::cout << "!Estimando y guardando distribución de los datos!" << std::endl;

    std::filesystem::path folder_route = "data/kmers_distribution/NB_"+
                    std::to_string(n_buckets)+"_BC_"+std::to_string(buckets_capacity)+"_CS_"
                    +std::to_string(compactor_size);

    try{
        if (not std::filesystem::exists(folder_route)){
            std::filesystem::create_directories(folder_route);
        }
    } catch (const std::filesystem::filesystem_error& e){
        std::cerr << "Error creating folder: " << e.what() << std::endl;
    }
    
    std::string csvFilename = folder_route.string()+"/"+std::to_string(k_) + "mers_distribution.csv";
                
    std::ofstream csvFile(csvFilename);
    if (csvFile.is_open()){
        csvFile << "quantile,real_quantile,estimated_quantile,rank,real_rank,estimated_rank\n";

        // Obtiene los cuantiles estimados y guarda los datos en un csv
        uint64_t estimated_quantile, j = 0;
        size_t estimated_rank;
        for (double i=0.0 ; i <= 1.00000 ; i+=quantile_ratio){
            estimated_quantile = sketch.quantile(i);
            estimated_rank = sketch.rank(kmers[real_ranks[j] - 1]);

            csvFile << i << "," << real_quantiles[j] << "," << estimated_quantile << 
            "," << kmers[real_ranks[j] - 1] << "," << real_ranks[j] << "," << estimated_rank << std::endl;

            j++;
        }

        csvFile.close();
        std::cout << "!Datos guardados exitosamente!" << std::endl;
    }else{
        std::cerr << "Error creando CSV: " << csvFilename << std::endl;
    }
    
    std::cout << "!Almacenando memoria usada!" << std::endl;

    csvFilename = folder_route.string()+"/"+std::to_string(k_) + "mers_memory.csv";
    std::ofstream csvFile1(csvFilename);
    if (csvFile1.is_open()){
        csvFile1 << "elements,unique_elements,sketch_memory,vector_memory,compressed_vector_memory,n_buckets,b_capacity,comp_size,comp_factor\n";
        csvFile1 << total_kmers << "," << unique_elements << "," << sketch.memory() << "," << vector_memory << "," << cv_memory
         << "," << n_buckets << "," << buckets_capacity << "," << compactor_size << "," << compression_factor << std::endl;

        csvFile1.close();
    }else{
        std::cerr << "Error creando CSV: " << csvFilename << std::endl;
    }

}

#endif