#include "kll.cpp"


class CooledKLL{
    enum Component{
        HOT_FILTER = 1,
        CLASIC_KLL = 2
    };

    class Bucket{
    public:
        std::vector<int_t> items;
        std::vector<size_t> frequencys;
        size_t vote, capacity;
        
        Bucket(size_t capacity) :
            items(std::vector<int_t>()),
            frequencys(std::vector<size_t>())
        {
            this->capacity = capacity;
            this->vote = 0;        
        }

        /**
         * @brief Find the position of element in items's vector if it exists.
         * 
         * @param element Element to find.
         * @return size_t Position of the element in items or -1 if element is not present.
         */
        std::pair<size_t,bool> find(int_t element){
            size_t items_size = items.size();
            for (size_t i=0 ; i<items_size ; i++){
                if (items[i] == element) return std::make_pair(i, true);
            }
            return std::make_pair(0, false);
        }

        /**
         * @brief Find the element with the minimum frequency.
         * 
         * @return size_t Index of the element.
         */
        size_t find_minimum(){
            size_t lowest_idx, min = UINT64_MAX;
            size_t entry_size = items.size();
            for (size_t i=0 ; i<entry_size ; i++){
                if (frequencys[i] < min){
                    min = frequencys[i];
                    lowest_idx = i;
                }
            }
            return lowest_idx;
        }

        /**
         * @brief Determines used memory by the object.
         * 
         * @return size_t: used memory in bytes.
         */
        size_t memory(){
            size_t memory_used = 0;

            // used memory by the items in containers.
            memory_used += items.size() * sizeof(int_t);
            memory_used += frequencys.size() * sizeof(size_t);

            // static used memory
            memory_used += sizeof(*this);

            return memory_used;
        }
        
    };

private:
    KLL kll;
    std::vector<Bucket> buckets;
    size_t n_buckets;
    int eviction_threshold;
public:
    /**
     * @brief Construct a new Cooled-KLL object
     * 
     * @param n_buckets Number of buckets in the hot filter.
     * @param buckets_capacity Number of entries per bucket.
     * @param eviction_threshold Threshold to evict an element to the kll component.
     * @param k Capacity of the largest compactor in kll.
     * @param c Compression factor for creating a new compactor in the KLL.
     */
    CooledKLL(size_t n_buckets, size_t buckets_capacity, int eviction_threshold, int k, float c) : 
        kll(KLL(k, c)),
        buckets(std::vector<Bucket>(n_buckets, Bucket(buckets_capacity)))
    {
        this->eviction_threshold = eviction_threshold;
        this->n_buckets = n_buckets;
    }

    /**
     * @brief Inserts an element into the Hot filter side (buckets) or into the KLL as appropriate.
     * 
     * @param element Element to insert.
     */
    void insert(int_t element, size_t frequency = 1){
        size_t bucket_idx = hash(element);
        std::pair<size_t, bool> pair = buckets[bucket_idx].find(element);
        size_t entry_idx = pair.first;
        // If the element exits in the bucket, increase its frequency by 1
        if (pair.second){
            buckets[bucket_idx].frequencys[entry_idx] += frequency;
            return;
        }
        // If element is not in the bucket, but there is an empty entry, push the pair (element, 1) into an entry.
        if (buckets[bucket_idx].items.size() < buckets[bucket_idx].capacity){
            buckets[bucket_idx].items.push_back(element);
            buckets[bucket_idx].frequencys.push_back(frequency);
            return;
        }
        // If element is not in the bucket and there isn't any entry available, increase vote by 1 and continue
        buckets[bucket_idx].vote += 1;
        
        // Find the minimum element
        size_t lowest_idx = buckets[bucket_idx].find_minimum();
        size_t min = buckets[bucket_idx].frequencys[lowest_idx];

        // If vote/min_frequency < eviction_threshold, insert element into KLL.
        int condition = static_cast<int>(std::round(buckets[bucket_idx].vote/min));
        if (condition < eviction_threshold){
            kll.insert(element, frequency);
            return;
        // If not, replace minimum frequency element with incoming element and insert the evicted element into KLL
        } else {
            buckets[bucket_idx].vote = 0;
            int_t min_element = buckets[bucket_idx].items[lowest_idx];
            kll.insert(min_element, min);
            buckets[bucket_idx].items[lowest_idx] = element;
            buckets[bucket_idx].frequencys[lowest_idx] = frequency;
            return;
        }
    }

    /**
     * @brief Estimates the rank of the given element.
     * 
     * @param element Element to estimate the rank.
     * @return size_t Amount of elements that are less or equal to element.
     */
    size_t rank(int_t element){
        int_t rank = 0;
        size_t entries_number;

        // Go through each element of the Hot filter to estimate de Hot filter rank
        size_t buckets_number = buckets.size();
        for (size_t i=0 ; i<buckets_number ; i++){
            entries_number = buckets[i].items.size();
            for (size_t j=0 ; j<entries_number ; j++){
                if (buckets[i].items[j] <= element){
                    rank += buckets[i].frequencys[j];
                }
            }
        }

        rank += kll.rank(element);

        return rank;
    }

    int_t quantile(float delta){
        if (delta < 0 or 1 < delta){
            throw std::invalid_argument("delta must belong to [0, 1]");
        }
        std::vector<std::pair<int_t, size_t>> data_kll = kll.data();
        std::vector<std::pair<int_t, size_t>> data_hot_filter;
        size_t quantile_pos = 0;

        // Collect all the elements and its frequencys from the hot filter
        size_t buckets_number = buckets.size();
        for (size_t i=0 ; i<buckets_number ; i++){
            size_t bucket_size = buckets[i].items.size();
            for (size_t j=0 ; j<bucket_size ; j++){
                data_hot_filter.push_back(std::make_pair(buckets[i].items[j], buckets[i].frequencys[j]));
                quantile_pos += buckets[i].frequencys[j];
            }
        }

        // Sort the created vector
        std::sort(data_hot_filter.begin(), data_hot_filter.end(), [](const std::pair<int_t, size_t>& a, const std::pair<int_t, size_t>& b){
            return a.first < b.first;
        });

        quantile_pos += kll.getSketch_size();
        quantile_pos = static_cast<size_t>(std::round(delta * quantile_pos));

        size_t hot_filter_size = data_hot_filter.size(), kll_size = data_kll.size();
        size_t all_data_size =  hot_filter_size + kll_size;
        size_t i=0, j=0, count=0;
        int last_picked = -1;
        for (size_t x=0; x<all_data_size ; x++){
            if (count <= quantile_pos){
                if (i < hot_filter_size and j < kll_size){
                    if (data_hot_filter[i].first <= data_kll[j].first){
                        count += data_hot_filter[i].second;
                        i++;
                        last_picked = HOT_FILTER;
                    } else {
                        count += data_kll[j].second;
                        j++;
                        last_picked = CLASIC_KLL;
                    }
                } else if (i >= hot_filter_size){
                    count += data_kll[j].second;
                    j++;
                    last_picked = CLASIC_KLL;
                } else {
                    count += data_hot_filter[i].second;
                    i++;
                    last_picked = HOT_FILTER;
                }
            } else {
                if (last_picked == HOT_FILTER){
                    return data_hot_filter[i-1].first;
                }
                if (last_picked == CLASIC_KLL){
                    return data_kll[j-1].first;
                }
                if (data_hot_filter[0].first <= data_kll[0].first){
                    return data_hot_filter[0].first;
                }
                return data_kll[0].first;
            }
        }
        if (last_picked == HOT_FILTER){
            return data_hot_filter[hot_filter_size-1].first;
        }
        return data_kll[kll_size-1].first;
        
    }

    /**
     * @brief Determines used memory of the object.
     * 
     * @return size_t: used memory in bytes.
     */
    size_t memory(){
        size_t memory_used = kll.memory(), hot_filter_size = buckets.size();

        // used memory by the buckets and its elements.
        for (uint32_t i=0 ; i<hot_filter_size ; i++){
            memory_used += buckets[i].memory() + sizeof(buckets[i]);
        }

        // used static memory
        memory_used += sizeof(*this);

        return memory_used;
    }

private:
    size_t hash(size_t x) {
        x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
        x = x ^ (x >> 31);
        return (x % n_buckets);
    }
    
};

// int main(){
//     int num = 3;
//     float quantile = 0.2;
//     CooledKLL kll(4, 2, 4, 10, 0.6);
//     for (int i=1 ; i<= 100 ; i++){
//         for (int j=0 ; j<10 ; j++){
//             kll.insert(i);
//         }
//     }

//     int_t quantile_answer = kll.quantile(quantile);

//     std::cout << "Rank " << num << ": " << kll.rank(num) << std::endl;
//     std::cout << "Quantile " << quantile << ": " << quantile_answer << std::endl;
// }