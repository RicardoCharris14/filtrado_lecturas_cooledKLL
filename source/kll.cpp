#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <ctime>
#include <set>

typedef uint64_t int_t;
// 192 bits
class KLL{
private:
    std::vector<std::vector<int_t>> compactors;
    size_t height, sketch_size;
    float k, c;
public:
    /**
     * @brief Construct a new KLL sketch and initialize de first compactor
     * 
     * @param k Capacity of the biggest compactor.
     * @param c Constant for reducing the next compactor size when creating a new one.
     */
    KLL(int k, float c) : k(k), c(c){
        if (not (0.5 < c and c < 1)){
            throw std::invalid_argument("c must belong to (0.5, 1)");
        }
        if (k <= 0){
            throw std::invalid_argument("k must be greater than 0");
        }
        height = 0;
        compactors.emplace_back();
    }

    /**
     * @brief Insert an element into the sketch.
     * 
     * @param element Element to insert.
     */
    void insert(int_t element){
        compactors[height].push_back(element);
        compaction();
    }

    /**
     * @brief Insert multiple times an element into the sketch on an efficient way.
     * 
     * @param element Element to insert.
     * @param frequency Times the element must be inserted.
     */
    void insert(int_t element, size_t frequency){
        size_t rest = frequency, w;
        size_t exponent = static_cast<int>(log2(frequency));
        if (exponent > height) exponent = height;
        
        // Insert the element in the respective compactor following the binary representation of the number.
        while (rest > 0){
            w = weight(exponent);
            while (rest >= w){
                compactors[height - exponent].push_back(element);
                rest -= w;
            }
            exponent = static_cast<int>(log2(rest));
        }

        compaction();
    }

    /**
     * @brief Estimates the rank of the given element.
     * 
     * @param element Element to calculate the rank.
     * @return int_t Rank of the element.
     */
    size_t rank(int_t element){
        size_t rank = 0, compactor_rank = 0;
        size_t compactor_size;
        for (size_t i=0 ; i<=height ; i++){
            compactor_size = compactors[i].size();
            for (size_t j=0 ; j<compactor_size ; j++){
                if (compactors[i][j] <= element){
                    compactor_rank++;
                }
            }
            rank += compactor_rank * weight(height-i);
            compactor_rank = 0;
        }
        return rank;
    }
    

    /**
     * @brief Estimates the delta-quantile of the data in the sketch.
     * 
     * @param delta Quantile to estimate.
     * @return int_t Element that is the delta-quantile of the data.
     */
    int_t quantile(float delta){
        if (delta < 0 or 1 < delta){
            throw std::invalid_argument("delta must belong to [0, 1]");
        }
        size_t data_size, counter = 0, quantile_pos;
        std::vector<std::pair<int_t, size_t>> elements;
        elements = data();

        quantile_pos = static_cast<size_t>(std::round(delta * sketch_size));
        data_size = elements.size();
        for (size_t i=0 ; i<data_size ; i++){
            counter+=elements[i].second;
            if (counter > quantile_pos){
                if (i != 0) return elements[i-1].first;
                else return elements[i].first;
            }
        }
        return elements[data_size-1].first;
    }

    /**
     * @brief Creates a sorted vector with the elements in the sketch.
     * 
     * @return std::vector<std::pair<int_t, size_t>> Ordered vector with all the elements in the sketch.
     */
    std::vector<std::pair<int_t, size_t>> data(){
        size_t compactor_size, w;
        sketch_size = 0;

        std::vector<std::pair<int_t, size_t>> data;
        for (size_t i=0; i<=height ; i++){
            compactor_size = compactors[i].size();
            for (size_t j=0 ; j<compactor_size ; j++){
                w = weight(height-i);
                sketch_size += w;
                data.push_back(std::make_pair(compactors[i][j], w));
            }
        }

        std::sort(data.begin(), data.end(), [](const std::pair<int_t, size_t>& a, const std::pair<int_t, size_t>& b){
            return a.first < b.first;
        });

        return data;
    }

    /**
     * @brief Determines used memory by the object.
     * 
     * @return size_t: used memory in bytes.
     */
    size_t memory(){
        size_t memory_used = 0;

        // Calculates used memory by the elements themselves.
        uint32_t size = compactors.size();
        for (uint32_t i=0 ; i<size ; i++){
            memory_used += compactors[i].size() * sizeof(size_t);
        }

        // Calculates used memory by the internal elements of 'compactors' vector.
        memory_used += compactors.size() * sizeof(std::vector<size_t>);
        
        // static used memory
        memory_used += sizeof(*this);

        return memory_used;
    }

    size_t getSketch_size(){
        return sketch_size;
    }

private:
    /**
     * @brief compacts each compactor from bottom to up if it has reached its maximum capacity.
     * 
     */
    void compaction(){
        for (size_t i=0 ; i<=height ; i++){
            if (compactors[height-i].size() >= compactorCapacity(i)){
                sortCompactor(compactors[height-i]);
                if (i == height){
                    height++;
                    compactors.emplace_back();
                    compactLastLevel();
                } else{
                    compactLevel(compactors[height-i], height-i);
                }
            }
        }
    }

    /**
     * @brief Calculates compactor's capacity of the given level.
     * 
     * @param level Level of the compactor.
     * @return uint64_t Capacity of the compactor.
     */
    int_t compactorCapacity(size_t level){
        return std::max(static_cast<int>(std::round(k * std::pow(c, height - level))), 2);
    }

    /**
     * @brief Sort the given compactor.
     * 
     * @param compactor Compactor to sort.
     */
    void sortCompactor(std::vector<int_t>& compactor){
        std::sort(compactor.begin(), compactor.end());
    }

    /**
     * @brief Determinates randomly even or odd.
     * 
     * @return true If even.
     * @return false If odd.
     */
    bool even(){
        std::mt19937 gen(static_cast<unsigned int>(std::time(0)));
        std::uniform_int_distribution<int> dist(1, 2);
        int random_number = dist(gen);
        return random_number == 2;
    }

    /**
     * @brief Compacts the highest compactor maintaining half of the elements and deleting the other half.
     * 
     */
    void compactLastLevel(){
        size_t size = compactors[0].size(), i, idx = 0;

        // selects which elements will remain
        if (even()) i=0; else i=1;
        
        // moves selected elements to arrays begin
        for (; i<size ; i+=2){
            compactors[0][idx] = compactors[0][i];
            idx++;
        }

        // Delete unselected elements
        compactors[0].resize(idx);
    }

    /**
     * @brief Compacts the given compactor transferring half of the elements to the next level and deleting the other half.
     * 
     * @param compactor Compactor to compact.
     * @param idx Index of the compactor in the vector of compactors.
     */
    void compactLevel(std::vector<int_t>& compactor, size_t idx){
        size_t size = compactor.size(), i;
        std::vector<int_t> buffer;

        // selects which elements will remain
        if (even()) i=0; else i=1;
        
        // temporarily stores selected elements in the buffer
        for (; i<size ; i+=2){
            buffer.push_back(compactor[i]);
        }

        // clears the compactor and transfers selected elements to the next level
        compactor.clear();
        size = buffer.size();
        for(size_t i=0 ; i<size ; i++){
            compactors[idx-1].push_back(buffer[i]);
        }
    }

    /**
     * @brief Calculates the weight of a given level.
     * 
     * @param level Level to obtain the weight.
     * @return size_t Weight of the level.
     */
    size_t weight(size_t level){
        return std::pow(2, level);
    }
};

// int main(){
//     int num = 10;
//     KLL sketch(10, 0.6);
//     for(int i=0 ; i<50 ; i++){
//         sketch.insert(i);
//     }

//     std::cout << "Rank de " << num << ": "<< sketch.rank(10) << std::endl;
// }