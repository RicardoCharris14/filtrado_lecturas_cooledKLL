#include "../include/lectorGenomas.hpp"
#include "cooled-kll.cpp"
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <math.h>
#include <filesystem>

// Mapa de caracteres a 2 bits: A=00, C=01, G=10, T=11
uint64_t charToBits(char c) {
    switch(c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4; // Error
    }
}

// Convierte una string en bits
uint64_t stringToBits(std::string kmer) {
    uint64_t bits = 0;
    for (int i=0 ; i<kmer.size() ; i++){
        bits = (bits << 2) | charToBits(kmer[i]);
    }
    return bits;
}

// Decodificar bits a string (solo para guardar en CSV)
std::string bitsToString(uint64_t kmer, int k) {
    std::string s(k, ' ');
    for (int i = 0; i < k; ++i) {
        uint64_t base = (kmer >> (2 * (k - 1 - i))) & 3;
        const char bases[] = "ACGT";
        s[i] = bases[base];
    }
    return s;
}

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
        csvFile1 << "total_elements,sketch_memory,vector_memory,compressed_vector_memory,n_buckets,b_capacity,comp_size,comp_factor\n";
        csvFile1 << total_kmers << "," << sketch.memory() << "," << vector_memory << "," << cv_memory
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
        csvFile1 << "total_elements,sketch_memory,vector_memory,compressed_vector_memory,n_buckets,b_capacity,comp_size,comp_factor\n";
        csvFile1 << total_kmers << "," << sketch.memory() << "," << vector_memory << "," << cv_memory
         << "," << n_buckets << "," << buckets_capacity << "," << compactor_size << "," << compression_factor << std::endl;

        csvFile1.close();
    }else{
        std::cerr << "Error creando CSV: " << csvFilename << std::endl;
    }

}

std::vector<std::pair<uint64_t, uint64_t>> leerKmers(const std::string& kmers_path) {
    std::vector<std::pair<uint64_t, uint64_t>> kmers_dist;
    std::ifstream file(kmers_path);

    // Verificación de error
    if (!file.is_open()) {
        std::cerr << "Error: No se pudo abrir el archivo en la ruta: " << kmers_path << std::endl;
        std::exit(1);
    } 

    std::string line;

    // Omitir la primera línea (cabecera "kmer,frequency")
    if (file.good()) {
        std::getline(file, line);
    }

    // Leer línea por línea
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string kmer;
        std::string kmer_frequency;

        // leemos el kmer y su frecuencia
        if (std::getline(ss, kmer, ',')) {
            if (std::getline(ss, kmer_frequency)) {
                uint64_t frequency = std::stoul(kmer_frequency);
                kmers_dist.push_back({std::stoull(kmer), frequency});
            }
        }
    }

    file.close();
    
    return kmers_dist;
}

int main(int argc, char* argv[]){
    // Verificacion de correctitud en la ejecucion del programa
    if (argc != 7){
        std::cerr << "correct usage: ./exe <kmers_file> <k-mers_length> <distribution> <N_buckets> <B_capacity> <C_size>" << std::endl;
        std::cerr << "<kmers_file>: path to the file with the kmers." << std::endl;
        std::cerr << "<k-mers_length>: length of kmers." << std::endl;
        std::cerr << "<distribution>: define type of distribution: 0 = kmers distribution | 1 = frequency distribution." << std::endl;
        std::cerr << "<N_buckets>: number of buckets in the hot filter part of the sketch." << std::endl;
        std::cerr << "<B_capacity>: number of entries a bucket have." << std::endl;
        std::cerr << "<C_size>: number of elements of the largest compactor in the classic kll part." << std::endl;
        return 1;
    }

    std::string kmers_path;
    int k_;
    bool frequency_distribution;
    double quantile_ratio = 0.001;
    // KLL settings
    size_t n_buckets, buckets_capacity; 
    int compactor_size; 
    float compression_factor = 0.7;

    // Verificacion de pertinencia de los argumentos
    try{
        kmers_path = argv[1];
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

    std::vector<std::pair<uint64_t, uint64_t>> kmers_dist = leerKmers(kmers_path);

    if (frequency_distribution){
        frequencyExperiments(kmers_dist, k_, quantile_ratio, n_buckets, buckets_capacity, compactor_size, compression_factor);
    } else {
        kmersExperiments(kmers_dist, k_, quantile_ratio, n_buckets, buckets_capacity, compactor_size, compression_factor);
    }
    
}