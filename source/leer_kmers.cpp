#include <algorithm>
#include "../include/procesarKmers.hpp"


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

int main(int argc, char* argv[]){
    if (argc != 2){
        std::cerr << "correct usage: ./exec <folder_url> <k>" << std::endl;
        std::cerr << "<folder_url>: path to the folder where FASTA files are located." << std::endl;
        std::cerr << "<k>: length of the kmer." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    int k = std::stoi(argv[2]);
    std::string folder_url = argv[1];
    if (k <= 0 or k > 31){
        std::cerr << "<k> must be lower or equal than 31 and greater than 0" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Obtiene los kmers a partir de las lecturas presentes en la carpeta indicada
    std::vector<std::pair<uint64_t, size_t>> kmers_distribution = procesarKMers(folder_url, k);

    // Ordena los kmers en base a su representacion binaria.
    std::sort(kmers_distribution.begin(), kmers_distribution.end());

    std::cout << "Guardando resultados en CSV" << std::endl;
    // Guardar en CSV

    // Crea el directorio si no existe.
    std::filesystem::path folder_route = "data/kmers";
    try{
        if (not std::filesystem::exists(folder_route)){
            std::filesystem::create_directories(folder_route);
        }
    } catch (const std::filesystem::filesystem_error& e){
        std::cerr << "Error creating folder: " << e.what() << std::endl;
    }


    std::string csvFilename = folder_route.string()+"/" + std::to_string(k) + "mers_frequency.csv";
    std::ofstream csvFile(csvFilename);
    
    if (csvFile.is_open()) {
        csvFile << "kmer,frequency\n";
        
        for (size_t i = 0; i < kmers_distribution.size(); i++) {
            csvFile << kmers_distribution[i].first << "," 
                   << kmers_distribution[i].second << std::endl;
        }
        
        csvFile.close();
        std::cout << "Guardado en: " << csvFilename << std::endl;
        std::cout << "Registros: " << kmers_distribution.size() << std::endl;
    } else {
        std::cerr << "Error creando CSV: " << csvFilename << std::endl;
    }
}