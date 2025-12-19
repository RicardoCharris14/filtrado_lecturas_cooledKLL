#include <filesystem>
#include <fstream>

#include "../include/procesarKmers.hpp"
#include "cooled-kll.cpp"

int main(int argc, char* argv[]){
    // Verificacion de correctitud en la ejecucion del programa
    if (argc != 9){
        std::cerr << "correct usage: ./exe <folder_file> <save_file> <k-mers_length> <N_buckets> <B_capacity> <C_size> <l_quantile> <u_quantile>" << std::endl;
        std::cerr << "<folder_file>: path to the folder with genomic lectures of FASTA type." << std::endl;
        std::cerr << "<save_file>: path to the file where statistics of filtering will be saved." << std::endl;
        std::cerr << "<k-mers_length>: length of kmers." << std::endl;
        std::cerr << "<N_buckets>: number of buckets in the hot filter part of the sketch." << std::endl;
        std::cerr << "<B_capacity>: number of entries a bucket have." << std::endl;
        std::cerr << "<C_size>: number of elements of the largest compactor in the classic kll part." << std::endl;
        std::cerr << "<l_quantile>: lower quantile to filter data." << std::endl;
        std::cerr << "<u_quantile>: upper quantile to filter data." << std::endl;
        return 1;
    }

    std::string folder_path;
    std::string file_path;
    int k;
    // KLL settings
    size_t n_buckets, buckets_capacity; 
    int compactor_size; 
    float compression_factor = 0.7;
    int eviction_threshold = 16;

    // Filter settings
    float lower_quantile;
    float upper_quantile;


    // Verificacion de pertinencia de los argumentos
    try{
        folder_path = argv[1];
        file_path = argv[2];
        k = std::stoi(argv[3]);
        n_buckets = std::stoll(argv[4]);
        buckets_capacity = std::stoll(argv[5]);
        compactor_size = std::stoll(argv[6]);
        lower_quantile = std::stof(argv[7]);
        upper_quantile = std::stof(argv[8]);

        if (k <= 0 or k > 31){
            std::cerr << "<k-mer length> must be a number belonging to [1, 31]." << std::endl;
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
        if (lower_quantile <= 0 or lower_quantile >= 1 or upper_quantile <= 0 or upper_quantile >= 1){
            std::cerr << "<l_quantile> and <r_quantile> must belong to ]0,1[" << std::endl;
        }
    } catch (std::exception e){
        std::cerr << "Error: " << e.what() << std::endl;
        std::exit(1);
    }

    std::cout << "!Leyendo kmers!" << std::endl;

    std::vector<std::pair<uint64_t, size_t>> kmers = procesarKMers("Genomas", k);

    std::cout << "!Creando el sketch!" << std::endl;
    size_t total_kmers = kmers.size();
    CooledKLL sketch(n_buckets, buckets_capacity, eviction_threshold, compactor_size, compression_factor);
    for (size_t i=0 ; i<total_kmers ; i++){
        sketch.insert(kmers[i].second);
    }

    size_t lower_bound = sketch.quantile(lower_quantile);
    size_t upper_bound = sketch.quantile(upper_quantile);
    
    std::cout << "Se eliminaran los K-mers con abundancia menor a " << lower_bound << " y mayor a " << upper_bound << "." << std::endl;
    sketch.~CooledKLL();

    std::sort(kmers.begin(), kmers.end(), [](const auto& a, const auto& b){
        return a.second < b.second;
    });

    size_t kmers_eliminados = 0;
    size_t kmers_diferentes_eliminados = 0;
    size_t total_elements = 0;
    for (size_t i=0 ; i<total_kmers ; i++){
        if (kmers[i].second < lower_bound or kmers[i].second > upper_bound){
            kmers_eliminados += kmers[i].second;
            kmers_diferentes_eliminados += 1;
        }
        total_elements += kmers[i].second;
    }

    // 1. Manejo de directorios (Crea la carpeta si no existe)
    std::filesystem::path path_obj(file_path);
    if (path_obj.has_parent_path()) {
        std::filesystem::create_directories(path_obj.parent_path());
    }

    // 2. Verificar si el archivo ya existe para saber si poner cabecera
    bool file_exists = std::filesystem::exists(path_obj);

    // 3. Abrir archivo en modo append (std::ios::app)
    std::ofstream outfile(file_path, std::ios::app);

    if (outfile.is_open()) {
        // Si el archivo no existía, escribimos la cabecera primero
        if (!file_exists) {
            outfile << "k,lower_quantile,upper_quantile,lower_bound,upper_bound,elements,unique_elim_e,elim_e\n";
        }

        // Escribimos los datos
        outfile << k << ","
                << lower_quantile << ","
                << upper_quantile << ","
                << lower_bound << ","
                << upper_bound << ","
                << total_elements << ","
                << kmers_diferentes_eliminados << ","
                << kmers_eliminados << "\n";
        
        outfile.close();
        std::cout << "Resultados guardados exitosamente en: " << file_path << std::endl;
    } else {
        std::cerr << "Error: No se pudo abrir o crear el archivo en " << file_path << std::endl;
    }

}