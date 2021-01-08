//
// Created by jarek on 30.12.20.
//

#ifndef DBSCAN_CACHE_HPP
#define DBSCAN_CACHE_HPP

#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <functional>
#include <chrono>

using namespace std;
using namespace std::chrono;

/**
 * @brief Special cache object with very specific pruning strategy. Every value has 2 arithmetic identifiers.
 * For example these can be identifiers of 2 different objects, and the value can be a result of comparing these 2
 * objects. We assume that the user of this cache knows that at some point a subset of values stored in the cache can be
 * removed.
 * @tparam T Type of the stored value.
 */
template<class T>
class Cache {
protected:
    /**
     * @brief A hash table of hash tables. We store a value identified by 2 identifiers. The first one is smaller, than
     * the second one.
     */
    unordered_map<int, unordered_map<int, T>> data;

    /**
     * @brief Heap like structure with the smallest element at the top. It is used to allow fast removal of elements
     * from <b>data</b> with indexes less than some given number. It is a pointer, because priority_queue doeas not
     * have a clear() method, instead we delete and create the new structure as opposed to pop()-ing all elements from
     * the structure. Deleting and creating should be faster.
     */
    priority_queue<int, vector<int>, greater<int>> *indexes;

    /**
     * @brief Internal counter of values with unique pairs of identifiers stored in the cache.
     */
    size_t elements_counter;

    /**
     * @brief The maximum number of values that were stored in the cache from the cache creation or call to
     * <b>clear</b> method.
     */
    size_t max_elements_counter;

    /**
     * @brief A counter of cache misses - when the value is not in the cache.
     */
    mutable int misses_counter;

    /**
     * @brief A counter of cache hits - when the value was found in the cache.
     */
    mutable int hits_counter;

    /**
     * @brief Internal helper parameter used to measure time spent in cache by the caller.
     */
    mutable time_point<high_resolution_clock> start_time;

    /**
     * @brief Counter of the total time spent in cache by the caller.
     */
    mutable duration<double> time_of_operation;

    /**
     * @brief Internal: start measuring time spent in cache.
     */
    void start() const;

    /**
     * @brief Internal: stop measuring time spent in cache and update total time counter.
     */
    void stop() const;

    /**
     * @brief Increments emelents counter by 1 and possibly updates maximum element counter.
     */
    void incrementElementsCounter();

public:
    /**
     * @brief Constructor. Resets all counters and allocates <b>indexes</b> structure.
     */
    Cache();

    /**
     * @brief Destructor. Frees <b>indexs</b> structure.
     */
    ~Cache();

    /**
     * @brief Retrieves a value from the cache identified by 2 indexes.
     * @param i First index of the value to be retrieved.
     * @param j Second index of the value to be retrieved.
     * @param value Returned, retrieved value.
     * @return True if value was found, false otherwise.
     */
    bool get(int i, int j, T &value) const;

    /**
     * @brief Adds value identified by two indexes to the cache. If value already exists in the cache, then
     * the operation i NO OP.
     * @param i First index of the value to be added to the cache.
     * @param j Second index of the value to be added to the cache.
     * @param value Value to be added.
     */
    void add(int i, int j, T value);

    /**
     * @brief Returns number of currently stored elements in the cache.
     * @return Count of elements in the cache.
     */
    size_t size() const;

    /**
     * @brief Returns the maximum number of elements that were stored in the cache from the creation of the cache or
     * from last call to <b>clear</b> method.
     * @return Count of elements in the cache.
     */
    size_t maxSize() const;

    /**
     * @brief Getter for number of misses counted from the creation of the cache or from last call to <b>clear</b>
     * method.
     * @return Value of internal <b>misses_counter</b>.
     */
    int misses() const;

    /**
     * @brief Getter for number of hits counted from the creation of the cache or from last call to <b>clear</b>
     * method.
     * @return Value of internal <b>hits_counter</b>.
     */
    int hits() const;

    /**
     * @brief Removes from the cache all elements which smaller index is smaller than a given value.
     * @param i A threshold to compare with the smaller indices of the values to decide if value should be removed from
     * the cache or not.
     */
    void purgeSmaller(int i);

    /**
     * @brief Removes all data and resets all counters.
     */
    void clear();

    /**
     * @brief Getter for time duration of the operations (get, add purgeSmaller) on the cache.
     * @return Value of the internal counter <b>time_of_operation</b>.
     */
    duration<double> time() const;
};

//------------------ IMPLEMENTATION --------------------

template<class T>
Cache<T>::Cache() : elements_counter(0), max_elements_counter(0), misses_counter(0), hits_counter(0),
                    time_of_operation(duration<double>::zero())
{
    indexes = new priority_queue<int, vector<int>, greater<int>>();
}

template<class T>
Cache<T>::~Cache()
{
    delete indexes;
}

template<class T>
bool Cache<T>::get(int i, int j, T &value) const
{
    start();
    if(i == j) {
        stop();
        misses_counter++;
        return false;
    } else if(i > j) {
        swap(i, j);
    }
    auto iter_1 = data.find(i);
    if(iter_1 == data.end()) {
        misses_counter++;
        stop();
        return false;
    }
    auto &map = iter_1->second;
    auto iter_2 = map.find(j);
    if(iter_2 == map.end()) {
        misses_counter++;
        stop();
        return false;
    }
    value = iter_2->second;
    hits_counter++;
    stop();
    return true;
}

template<class T>
void Cache<T>::add(int i, int j, T value)
{
    start();
    if(i == j) {
        stop();
        return;
    } else if(i > j) {
        swap(i, j);
    }
    auto iter_1 = data.find(i);
    if(iter_1 == data.end()) {
        unordered_map<int, T> map;
        map.emplace(j, value);
        data.emplace(i, map);
        indexes->push(i);
        incrementElementsCounter();
    } else {
        auto ret = iter_1->second.emplace(j, value);
        if(ret.second) {
            incrementElementsCounter();
        }
    }
    stop();
}

template<class T>
size_t Cache<T>::size() const
{
    return elements_counter;
}

template<class T>
size_t Cache<T>::maxSize() const
{
    return max_elements_counter;
}

template<class T>
int Cache<T>::hits() const
{
    return hits_counter;
}

template<class T>
int Cache<T>::misses() const
{
    return misses_counter;
}

template<class T>
duration<double> Cache<T>::time() const
{
    return time_of_operation;
}

template<class T>
void Cache<T>::purgeSmaller(int i)
{
    start();
    while(!indexes->empty() && indexes->top() < i) {
        auto iter = data.find(indexes->top());
        elements_counter -= iter->second.size();
        data.erase(iter);
        indexes->pop();
    }
    stop();
}

template<class T>
void Cache<T>::clear()
{
    data.clear();
    delete indexes;
    indexes = new priority_queue<int, vector<int>, greater<int>>();
    elements_counter = 0;
    max_elements_counter = 0;
    hits_counter = 0;
    misses_counter = 0;
}

template<class T>
void Cache<T>::incrementElementsCounter()
{
    elements_counter++;
    if(elements_counter > max_elements_counter) {
        max_elements_counter = elements_counter;
    }
}

template<class T>
void Cache<T>::start() const
{
    start_time = high_resolution_clock::now();
}

template<class T>
void Cache<T>::stop() const
{
    time_of_operation += (high_resolution_clock::now() - start_time);
}

#endif //DBSCAN_CACHE_HPP
