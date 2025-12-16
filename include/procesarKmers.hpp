#ifndef PROCESAR_K_MERS_HPP
#define PROCESAR_K_MERS_HPP

#include <unordered_map>
#include <iostream>
#include "../include/lectorGenomas.hpp"

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

// Obtener canónico usando operaciones de bits (mucho más rápido)
uint64_t getCanonicalKmerBits(uint64_t kmer, int k) {
    uint64_t revComp = 0;
    uint64_t temp = kmer;
    for (int i = 0; i < k; ++i) {
        uint64_t base = temp & 3;
        uint64_t comp = base ^ 3; // 00<->11 (A-T), 01<->10 (C-G)
        revComp = (revComp << 2) | comp;
        temp >>= 2;
    }
    return std::min(kmer, revComp);
}

// Función para obtener los k-mers a partir de las lecturas en la carpeta indicada (sin ordenar)
std::vector<std::pair<uint64_t, size_t>> procesarKMers(std::string folder, int k) {
    if (k > 31) throw std::runtime_error("K demasiado grande para uint64_t (max 31)");
    std::cout << "\n=== Lectura de archivos iniciada ===" << std::endl;

    LectorGenomas reader(folder);
    std::unordered_map<uint64_t, uint64_t> kmers_frequency;
    
    uint64_t currentKmer = 0;
    uint64_t mask = (k == 32) ? ~0ULL : ((1ULL << (2 * k)) - 1); // Máscara para mantener solo k bases
    int basesInWindow = 0;

    // --- VARIABLES PARA EL PROGRESO ---
    size_t total_processed = 0;
    const size_t PRINT_INTERVAL = 1000000; // Imprimir cada 1 millón de k-mers
    // ----------------------------------

    std::cout << "=== Iniciando procesamiento de k-mers (k=" << k << ") ===" << std::endl;
    std::cout << "Leyendo archivos del directorio 'Genomas'..." << std::endl;
    
    // Iteramos usando lógica de sliding window (ventana deslizante)
    // Esto evita substr() y re-leer strings.
    while (reader.hasMoreKmers(1)) { // Chequeamos si queda al menos 1 base
        char base = reader.getBaseAt(reader.getCurrentPosition());
        uint64_t val = charToBits(base);

        reader.advancePosition(); // Avanzar manualmente

        if (val > 3) { 
            // Base inválida (N o similar), reiniciar ventana
            basesInWindow = 0;
            currentKmer = 0;
            continue;
        }

        // Shift a la izquierda y añadir nueva base
        currentKmer = ((currentKmer << 2) | val) & mask;
        basesInWindow++;

        if (basesInWindow >= k) {
            // Tenemos un k-mer válido en currentKmer
            uint64_t canonical = getCanonicalKmerBits(currentKmer, k);
            kmers_frequency[canonical]++;
            // --- BLOQUE DE IMPRESIÓN DE PROGRESO ---
            total_processed++;
            if (total_processed % PRINT_INTERVAL == 0) {
                std::cout << "\r[Progreso] Procesados: " << (total_processed / 1000000) << "M"
                          << " | Unicos: " << kmers_frequency.size()
                          << " | Archivo: " << reader.getCurrentFilename() 
                          << "          " << std::flush; // Espacios extra para limpiar residuos visuales
            }
        }
        
        // Si cambiamos de archivo, reader.advancePosition() lo maneja, 
        // pero debemos reiniciar la ventana porque los k-mers no cruzan archivos
        if (reader.getCurrentPosition() == 0) {
            basesInWindow = 0;
            currentKmer = 0;
        }
    }

    std::cout << "\n\n=== Procesamiento Finalizado ===" << std::endl;
    std::cout << "Total k-mers procesados: " << total_processed << std::endl;
    std::cout << "Total k-mers unicos: " << kmers_frequency.size() << std::endl;
    std::cout << "Generando vector de resultados..." << std::endl;

    // Convertir mapa a vector para ordenar y devolver
    std::vector<std::pair<uint64_t, size_t>> kmers_distribution;
    kmers_distribution.reserve(kmers_frequency.size());

    for (const auto& pair : kmers_frequency) {
        kmers_distribution.push_back({pair.first, pair.second});
    }

    return kmers_distribution;
}

#endif