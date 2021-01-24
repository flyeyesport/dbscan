//
// Created by jarek on 18.12.20.
//

#ifndef DBSCAN_DBSCAN_HPP
#define DBSCAN_DBSCAN_HPP


#include <vector>
#include <chrono>
#include <boost/function_output_iterator.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <type_traits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <boost/filesystem.hpp>
#include <algorithm>
#include "DataPoint.hpp"
#include "ClusterPoint.hpp"
#include "Cache.hpp"

using namespace std;
using namespace std::chrono;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;


/**
 * @brief DBSCAN algorithm implementation with extensions:
 * - border points can be assigned to only one or to many clusters
 * - different index options to speed up calculations
 * @tparam T Type of the dimensions of the vectors of data to be clustered. It must be an arithmetic type.
 * @tparam Dim Number of dimensions of data vectors.
 * @tparam RStarMaxNodeElements Maximum number of elements stored in single node in R*-tree index (when it is used).
 */
template<class T, int Dim, size_t RStarMaxNodeElements = 16>
class DBSCAN {
protected:
    /**
     * @brief Internal type of the data point used in R*-tree index.
     */
    typedef bg::model::point<T, Dim, bg::cs::cartesian> RPoint;

    /**
     * @brief Internal type with a mapping between points in R*-tree index and our internal collection of points.
     */
    typedef std::pair<RPoint, int> RRecord;

    /**
     * @brief Name of the file from which data was read. It is written to the "stat" file when results are saved.
     */
    string input_file_name;

    /**
     * @brief Name of the dataset. It can be set when reading data. It is used in file names when saving results.
     */
    string dataset_name;

    /**
     * @brief One of two main parameters of the algorithm. It is the minimum number of points that are close enough
     * (@see <b>eps</b>) to a given point to classify this point as a core point.
     */
    int min_neighbours;

    /**
     * @brief One of two main parameters of the algorithm. It is the radius in which we search for neighbouring
     * points for a given point. If a number of neighbours is less or equal some threshold (@see <b> min_neighbours</b>)
     * then we treat such point as a core point.
     */
    T eps;

    /**
     * @brief Boolean flag used to select one of two versions of the algorithm: DBSCAN (value set to false), and DBSCAN+
     * (value set to true). DBSCAN+ is an extension of DBSCAN in which border points can be assigned to more than one
     * cluster.
     */
    bool border_points_in_many_clusters;

    /**
     * @brief Distance used to calculate distances between data points.
     */
    function<T(const DataPoint<T, Dim> &, const DataPoint<T, Dim> &)> distance_func;

    /**
     * @brief Some distance functions (like for example cosine distance) can't operate with zero data point (vectors with all
     * dimensions set to 0). This boolean flag allows us to remove such points from the set before clustering. All
     * removed points are stored in <b>zeros</b> list.
     */
    bool treat_zero_as_noise;

    /**
     * @brief List of zero data points (vectors with all dimensions set to 0). @see <b>treat_zero_as_noise</b>.
     */
    vector<pair<int, const DataPoint<T, Dim> *>> zeros;

    /**
     * @brief Boolean flag indicating if all data point should be normalized (by changing lengths of vectors to 1)
     * before clustering.
     */
    bool normalize_values;

    /**
     * @brief Boolean flag indicating if all data point should be normalized (by calculating z-scores of vectors)
     * before clustering.
     */
    bool calculate_z_score;

    /**
     * @brief Boolean flag indicating if R*-tree index should be used. If it is set to true, then triangular inequality
     * index will not be used. It is possible to not use any of these indexes.
     */
    bool use_rstar_tree;

    /**
     * @brief Boolean flag indicating if triangular inequality index should be used. If it is set to true, then R*-tree
     * index will not be used. It is possible to not use any of these indexes.
     */
    bool use_triangular_inequality_index;

    /**
     * @brief R*-tree index. Used when <b>use_rstar_tree</b> flag is set to true.
     */
    bgi::rtree<RRecord, bgi::rstar<RStarMaxNodeElements>> index;

    /**
     * @brief Center point used by triangular inequality index. It is used when <b>use_triangular_inequality_index</b>
     * flag is set to true.
     */
    DataPoint<T, Dim> ti_center;

    /**
     * @brief Timer: number of seconds used for reading data from a file.
     */
    duration<double> reading_input_file_time;

    /**
     * @brief Timer: number of seconds used for calculating distances between all points in dataset and the center
     * point (@see <b>ti_center</b>) used only if triangular inequality index is used
     * (@see <b>use_triangular_inequality_index</b>).
     */
    duration<double> calc_dist_to_ti_center_time;

    /**
     * @brief Timer: number of seconds used to normalize data vectors. Only used when
     * <b>normalize_values</b> or <b>calculate_z_score</b> flag is set to true.
     */
    duration<double> calc_norm_time;

    /**
     * @brief Timer: number of seconds used to sort data points based on distances to center point
     * (@see <b>ti_center</b>). It is used only if <b>use_triangular_inequality_index</b> flag is set to true.
     */
    duration<double> sort_ti_time;

    /**
     * @brief Timer: number of seconds used to calculate distances between data points.
     */
    duration<double> calc_dist_time;

    /**
     * @brief Simple counter which incremented value is added to the name of the saved files with results. Each time
     * a save operation is performed the counter is incremented. It allows to save results of the algorithm on the same
     * data with slightly different settings in files with names that differ only by this counter.
     */
    int save_counter;

    /**
     * @brief Cache object used optionally for triangle inequality index. To use it set <b>use_cache</b> flag to true.
     */
    Cache<T> cache;

    /**
     * @brief Boolean flag indicating if cache for triangular inequality index will be used.
     */
    bool use_cache;

    /**
     * @brief Path to directory where files with results will be written.
     */
    string output_dir;

    /**
     * @brief Boolean flag indicating if epsilon algorithm parameter should be guessed/calculated or should be used a
     * passed <b>eps</b> parameter.
     */
    bool guess_epsilon;

    /**
     * @brief Boolean flag indicating if distance function returns values that when smaller indicate that data points
     * are closer to each other or the opposite. For example Euclidean distance is smaller for closer points but cosine
     * distance is bigger for closer points.
     */
    bool is_less_closer;

    /**
     * @brief Counter of the number of calculations of the distances between points and the ti_center point when
     * building a TI index.
     */
    int dist_calc_for_index_counter;

    /**
     * @brief Internal method returning a set of points that are in the radius of <b>eps</b> from the given point.
     * The method may use different indexes to speed up searches. @see <b>use_rstar_tree</b> and
     * <b>use_triangular_inequality_index</b>. A distance function <b>distance_func</b> is used to calculate distances
     * between points.
     * @param points A set of points to be searched for neighbours.
     * @param center_index An index of the point in a set <b>points</b> that is given as the center.
     * @return A list of indexes of the points that are no farther away (in a sense of a distance function
     * <b>distance_func</b>) than <b>eps</b>. The list does not contain <p>center_index</b> value.
     */
    vector<int> regionQuery(vector<ClusterPoint<DataPoint<T, Dim>, T>> &points, int center_index);

    /**
     * @brief Internal method used to expand a cluster when some first core point of a cluster is found. The method
     * searches for new points and adds them to the cluster if they meet cluster conditions.
     * @param points A set of all data points in which the method searches for points to be included in a new cluster.
     * @param center_index Index of the first point of the cluster. It is an index of the point in the <b>points</b>
     * set.
     * @param neighbours A first set of indexes of the points that are neighbours of the center point. Each of these
     * points should be checked if it can be included into a cluster. For each of these points a new set of neighbours
     * is constructed to repeat eventual inclusion in the cluster. The process is repeated as long as there are
     * neighbours that were not checked.
     * @param cluster_number A number of the claster to be created. It is a label of the cluster.
     */
    void expandCluster(vector<ClusterPoint<DataPoint<T, Dim>, T>> &points, int center_index, vector<int> &neighbours,
                       int cluster_number);

    /**
     * @brief Internal method that verifies if two points are in the distance equal or less than <b>eps</b> from each
     * other. A distance function <b>distance_func</b> is used to calculate distances between points. For each points an
     * internal counter of comparisons is incremented.
     * @param p1 First point to be checked.
     * @param p2 Second point to be checked.
     * @return True if points are close enough, false otherwise.
     */
    bool closeEnough(ClusterPoint<DataPoint<T, Dim>, T> &p1, ClusterPoint<DataPoint<T, Dim>, T> &p2);

    /**
     * @brief Internal mehod which runs the DBSCAN algorithm in one of the selected variants. It labels all the data
     * points with cluster numbers and a category: noise, border point, core point.
     * @param data A set of points to be clustered.
     */
    void findClusters(vector<ClusterPoint<DataPoint<T, Dim>, T>> &data);

    /**
     * @brief Internal method that optionally normalizes data, filters out zero data points and builds one of the
     * indexes to speed up further processing. @see boolean flags of this class to learn possible configuration options
     * for this method.
     * @param data A set of data points to be optionally normalized, filtered and from which one of the indexes can be
     * built.
     * @return A set of extended data points ready to be used by clustering algorithm. This extension allows labeling
     * of points.
     */
    vector<ClusterPoint<DataPoint<T, Dim>, T>> buildIndex(const vector<pair<DataPoint<T, Dim>, vector<int>>> &data);

    /**
     * @brief Internal method that saves results of the clustering algorithm into 2 files:
     * - a csv file with the list of
     *   all data points with cluster numbers and some other extra data (original index, values of all dimensions,
     *   number of comparisons with other data points, class: noise, core point, border point), and cluster number
     *   (or -1 for noise points).
     * - a txt file with some of the satistics of and configuration parameters used by the clustering algorithm.
     * @param points Processed and clusterd points to be saved.
     * @param rand_index Rand index to be saved in stat file.
     * @param rand_tp Count of pairs that were in the same cluster in result and ground truth sets to be saved in stat
     * file.
     * @param rand_tn Count of pairs that were in other clusters in result and ground truth sets to be saved in stat
     * file.
     * @param rand_all Count of all pairs in ground truth set to be saved in stat file.
     * @param reading_input_file_time Duration of the reading process of data from file. It is saved in statistics text
     * file.
     * @param clustering_time Duration of the clustering process. It is saved in statistics text file.
     */
    void save(const vector<ClusterPoint<DataPoint<T, Dim>, T>> &points, double rand_index,
              int rand_tp, int rand_tn, int rand_all, duration<double> reading_input_file_time,
              duration<double> clustering_time);

    /**
     * @brief Calculates Rand index. It also sorts results and adds to the results points that were removed before
     * clustering because they were zeros.
     * @param data Set with ground truth data.
     * @param points Result set. It will be sorted by this method.
     * @param rand_tp Returned count of pairs that were in the same cluster in result and ground truth sets.
     * @param rand_tn Returned count of pairs that were in other clusters in result and ground truth sets.
     * @param rand_all Returned count of all pairs in ground truth set.
     * @return Rand index.
     */
    double sortResultsAndCalculateRandIndex(const vector<pair<DataPoint<T, Dim>, vector<int>>> &data,
                                            vector<ClusterPoint<DataPoint<T, Dim>, T>> &points, int &rand_tp,
                                            int &rand_tn, int &rand_all);

    /**
     * @brief Helper method which creates a structre of a point (RPoint) used in R*-tree index. It is compiled only for
     * points with 2 dimensions. For now R*-tree index can operate on point of 2 or 3 dimensions.
     * @tparam C Type of all the dimensions of the data point.
     * @param point Data point to be converted to the RPoint.
     * @return A new RPoint to be inserted into R*-tree index.
     */
    template<typename C = T, enable_if_t<is_equal_to_dim<C, 2, Dim>::value, int> = 0>
    RPoint createRPoint(const DataPoint<T, Dim> &point) const
    {
        return RPoint(point[0], point[1]);
    }

    /**
     * @brief Helper method which creates a structre of a point (RPoint) used in R*-tree index. It is compiled only for
     * points with 3 dimensions. For now R*-tree index can operate on point of 2 or 3 dimensions.
     * @tparam C Type of all the dimensions of the data point.
     * @param point Data point to be converted to the RPoint.
     * @return A new RPoint to be inserted into R*-tree index.
     */
    template<typename C = T, enable_if_t<is_equal_to_dim<C, 3, Dim>::value, int> = 0>
    RPoint createRPoint(const DataPoint<T, Dim> &point) const
    {
        return RPoint(point[0], point[1], point[2]);
    }

    /**
     * @brief Helper method which for data points with 16 should not be used. Use only for points with 2 or 3
     * dimensions.
     * @tparam C Type of all the dimensions of the data point.
     * @param point Data point to be converted to the RPoint.
     * @return A new RPoint to be inserted into R*-tree index.
     */
    template<typename C = T, enable_if_t<is_equal_to_dim<C, 16, Dim>::value, int> = 0>
    RPoint createRPoint(const DataPoint<T, Dim> &) const
    {
        throw runtime_error("Do not use R*-tree for other dimensions than 2 or 3!");
    }

    /**
     * @brief Internal method which sets <b>min_neighbours</b> to 5 and guesses the value of <b>eps</b>.
     * @param data Set of points to be to calculate best <b>eps</b>.
     * @param data Set of points to be to calculate best <b>eps</b>.
     * @return True if operation was successful, false otherwise.
     */
    bool guessEpsilon(vector<ClusterPoint<DataPoint<T, Dim>, T>> &data);

    /**
     * @brief
     * @param data
     * @param index
     * @param k
     * @return
     */
    T calcDistanceToKthClosest(vector<ClusterPoint<DataPoint<T, Dim>, T>> &data, int index, int k);

public:
    /**
     * @brief Constructor.
     * @param min_neighbours A threshold: minimum number of neighbours in a radius of <b>eps</b> to treat a point in
     * the center as a core point. The center point is counted to the number of neighbours.
     * @param eps A radius of the circle in which the algorithm searches for neighbours for each point candidate.
     * If the number of found neighbours is equal or more than the <b>min_neighbours</b> threshold, then the candidate
     * is marked as a core point.
     * @param distance_func A distance function used to calculate distances between data points.
     */
    DBSCAN(int min_neighbours, double eps,
           function<T(const DataPoint<T, Dim> &, const DataPoint<T, Dim> &)> distance_func);

    /**
     * @brief default destructor.
     */
    ~DBSCAN() = default;

    /**
     * @brief Runs a clustering algorithm and saves results into 2 files (data points with assigned cluster numbers
     * in csv format, and statistics and parameters used in txt format). It optionally normalizes data, removes
     * zero vectors, builds an index, clusters data and saves results.
     * @param data A set of data points to be clustered.
     * @return A set of clustered points. Data points are extended with cluster numbers and some additional info
     * (original index, values of all dimensions, number of comparisons with other data points, class: noise, core
     * point, border point), and cluster number (or -1 for noise points).
     */
    vector<ClusterPoint<DataPoint<T, Dim>, T>> findClusters(const vector<pair<DataPoint<T, Dim>, vector<int>>> &data);

    /**
     * @brief Reads data points from file. Each point is represented by a line in the file. Each line must have the
     * same number of values delimited by commas.
     * @param file_name Path to the file with data.
     * @param dataset_name A name of the dataset to be used to construct distinguishable file names with results.
     * @param data Returned list of pairs of data point and ground truth cluster index.
     * @return True if operation was successful, false otherwise.
     */
    bool
    readData(const string &file_name, const string &dataset_name, vector<pair<DataPoint<T, Dim>, vector<int>>> &data);

    /**
     * @brief Sets an internal flag to use a R*-tree index. All other indexes are switched off. Notice that this index
     * can be used only for data with 2 or 3 dimensions.
     */
    void useRStarTreeIndex();

    /**
     * @brief Sets an internal flag to use a triangular inequality index. All other indexes are switched off.
     */
    void useTriangularInequalityIndex();

    /**
     * @brief Switches off all of the indexes. No index is used to speed up searching for point's neighbours.
     */
    void useNoIndex();

    /**
     * @brief Sets an internal flag to normalize all data vectors to the lenght of 1. All other methods of
     * normalization are switched off.
     */
    void useDataNormalization();

    /**
     * @brief Sets an internal flag to NOT normalize all data vectors to the lenght of 1.
     */
    void useNoDataNormalization();

    /**
     * @brief Sets an internal flag to normalize all data vectors using z-score. All other methods of
     * normalization are switched off.
     */
    void useZScoreNormalization();

    /**
     * @brief Sets an internal flag to NOT normalize all data vectors using z-score.
     */
    void useNoZScoreNormalization();

    /**
     * @brief Sets an internal flag to filter out all zero vectors (vectors with all dimensions set to 0) from data to
     * be clustered. Such data points are treated as noise. It is useful for some distance functions that give undefined
     * results for such input data.
     */
    void treatZeroVectorsAsNoise();

    /**
     * @brief Sets an internal flag to NOT filter out all zero vectors (vectors with all dimensions set to 0) from data
     * to be clustered.
     */
    void treatZeroVectorsAsNotNoise();

    /**
     * @brief Sets an internal flag to run the clustering algorithm in a way that allows border points to be assigned
     * to many clusters. This version of algorithm is called DBSCAN+.
     */
    void assignBorderPointsToPossiblyManyClusters();

    /**
     * @brief Sets an internal flag to run the clustering algorithm in a way that allows border points to be assigned
     * to only one cluster. The order of processing the points determines to which cluster the point will be assigned.
     * This is the basic version of the algorithm, and is called DBSCAN.
     */
    void assignBorderPointsToOnlyOneCluster();

    /**
     * @brief Sets coordinates of the point used to calculate triangle inequality index. It can be any point.
     * By default it is set to the center of coordinate system.
     * @param point New coordinates to be set a a center point for triangle inequality index.
     */
    void setTriangularInequalityCenter(const DataPoint<T, Dim> &point);

    /**
     * @brief Sets new value of the <b>eps</b> algorithm parameter. @see <b>eps</b>.
     * @param new_epsilon New value to be set.
     */
    void setEpsilon(T new_epsilon);

    /**
     * @brief Sets new value of the <b>min_neighbours</b> algorithm parameter. @see <b>min_neighbours</b>.
     * @param new_min_neighbours New value to be set.
     */
    void setMinNeighbours(int new_min_neighbours);

    /**
     * @brief Sets an internal flag to use additional cache for triangle inequality index.
     */
    void useCache();

    /**
     * @brief Sets an internal flag to NOT use additional cache for triangle inequality index.
     */
    void useNoCache();

    /**
     * @brief Sets output directory into which files with results and statistics will be written after running the
     * algorithm.
     * @param folder A new path to the directory in which results will be saved.
     */
    void setOutputFolder(const string &folder);

    /**
     * @brief Sets internal flag to guess the value of the epsilon algorithm parameter instead of using some given value.
     */
    void setGuessEpsilon();

    /**
     * @brief Changes distance function used to calculate distances between data points to a new one.
     * @param distance A new function to be used.
     */
    void useDistanceFunction(function<T(const DataPoint<T, Dim> &, const DataPoint<T, Dim> &)> distance);

    /**
     * @brief Sets internal flag to define how to treat values returned by <b>distance_func</b>. Some distance
     * functions return bigger values for closer data points (for example cosine distance) and some return smaller
     * values for closer points (this is more intuitive, for example Euclidean distance).
     */
    void setSmallerIsCloser();

    /**
     * @brief Sets internal flag to define how to treat values returned by <b>distance_func</b>. Some distance
     * functions return bigger values for closer data points (for example cosine distance) and some return smaller
     * values for closer points (this is more intuitive, for example Euclidean distance).
     */
    void setBiggerIsCloser();
};

//--------------------------------------------- IMPLEMENTATION -------------------------------------------

template<class T, int Dim, size_t RStarMaxNodeElements>
DBSCAN<T, Dim, RStarMaxNodeElements>::DBSCAN(int min_neighbours, double eps,
                                             function<T(const DataPoint<T, Dim> &,
                                                        const DataPoint<T, Dim> &)> distance_func)
        : min_neighbours(min_neighbours), eps(eps), border_points_in_many_clusters(false),
          distance_func(distance_func), treat_zero_as_noise(false), normalize_values(false), calculate_z_score(false),
          use_rstar_tree(false), use_triangular_inequality_index(false), ti_center(1), // ti_center(vector<T>(Dim, 1)),
          reading_input_file_time(duration<double>::zero()),
          calc_dist_to_ti_center_time(duration<double>::zero()), calc_norm_time(duration<double>::zero()),
          sort_ti_time(duration<double>::zero()), calc_dist_time(duration<double>::zero()), save_counter(0),
          use_cache(false), guess_epsilon(false), is_less_closer(true), dist_calc_for_index_counter(0)
{
    if(min_neighbours < 2) {
        throw out_of_range("Minimum number of neighbours must be larger than 1.");
    }
}

template<class T, int Dim, size_t RStarMaxNodeElements>
vector<ClusterPoint<DataPoint<T, Dim>, T>>
DBSCAN<T, Dim, RStarMaxNodeElements>::buildIndex(const vector<pair<DataPoint<T, Dim>, vector<int>>> &data)
{
    vector<ClusterPoint<DataPoint<T, Dim>, T>> result;
    if(data.empty()) {
        return result;
    }
    index.clear();
    zeros.clear();
    cache.clear();
    calc_dist_time = duration<double>::zero();
    calc_norm_time = duration<double>::zero();
    calc_dist_to_ti_center_time = duration<double>::zero();

    int j = 0;
    for(int i = 0; i < (int) data.size(); i++) {
        const auto &point = data[i].first;
        if(!treat_zero_as_noise || !point.isZero()) {
            result.emplace_back(j++, i, &(point));
        } else {
            zeros.emplace_back(i, &(point));
        }
    }

    auto norm_start = high_resolution_clock::now();
    if(!result.empty()) {
        if(normalize_values) {
            for(auto &elem : result) {
                elem.normalize();
            }
        } else if(calculate_z_score) {
            DataPoint<T, Dim> avg;
            for(int i = 0; i < (int) result.size(); i++) {
                avg = avg + *(result[i].getData());
            }
            avg = avg / (T) result.size();

            DataPoint<T, Dim> std;
            for(int i = 0; i < (int) result.size(); i++) {
                auto val = *(result[i].getData()) - avg;
                val.abs();
                std = std + val;
            }
            std = std / (T) result.size();

            for(auto &elem : result) {
                elem.calculateZScore(avg, std);
            }
        } else {
            for(auto &elem : result) {
                elem.removeNormalization();
            }
        }
    }
    auto norm_stop = high_resolution_clock::now();
    calc_norm_time = norm_stop - norm_start;

    if(guess_epsilon) {
        if(!guessEpsilon(result)) {
            throw runtime_error("Guessing parameter values failed. Probably the data set has not enough data points.");
        } else {
            cout << "Guessed parameters: epsilon: " << eps << ", min-neigbours: " << min_neighbours << "." << endl;
        }
    }

    if(use_rstar_tree) {
        for(int i = 0; i < (int) result.size(); i++) {
            index.insert(make_pair(createRPoint(*(result[i].getData())), i));
        }
    } else if(use_triangular_inequality_index) {
        dist_calc_for_index_counter = 0;
        for(auto &elem : result) {
            auto ti_center_dist_start = high_resolution_clock::now();
            T dist = distance_func(ti_center, *(elem.getData()));
            auto ti_center_dist_stop = high_resolution_clock::now();
            calc_dist_to_ti_center_time += ti_center_dist_stop - ti_center_dist_start;
            elem.setTriangularInequalityDistance(dist);
            dist_calc_for_index_counter++;
        }

        auto sort_ti_start = high_resolution_clock::now();
        if(is_less_closer) {
            sort(result.begin(), result.end(), [](const ClusterPoint<DataPoint<T, Dim>, T> &lhs,
                                                  const ClusterPoint<DataPoint<T, Dim>, T> &rhs) -> bool {
                return lhs.getTriangularInequalityDistance() < rhs.getTriangularInequalityDistance();
            });
        } else {
            sort(result.begin(), result.end(), [](const ClusterPoint<DataPoint<T, Dim>, T> &lhs,
                                                  const ClusterPoint<DataPoint<T, Dim>, T> &rhs) -> bool {
                return lhs.getTriangularInequalityDistance() > rhs.getTriangularInequalityDistance();
            });
        }
        //reindex
        int i = 0;
        for(auto &elem : result) {
            elem.setDataIndex(i++);
        }
        auto sort_ti_stop = high_resolution_clock::now();
        sort_ti_time = sort_ti_stop - sort_ti_start;
    }
    return result;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
vector<int>
DBSCAN<T, Dim, RStarMaxNodeElements>::regionQuery(vector<ClusterPoint<DataPoint<T, Dim>, T>> &points,
                                                  int center_index)
{
    vector<int> result;
    auto &center = points[center_index];
    if(use_rstar_tree) {
        auto selector = boost::make_function_output_iterator([this, &result, &points, &center, &center_index](
                const bgi::rtree<RRecord, bgi::rstar<RStarMaxNodeElements>>::value_type &elem) {
            if(elem.second != center_index && closeEnough(points[elem.second], center)) {
                result.push_back(elem.second);
            }
        });

        bg::model::box<RPoint> query_box(createRPoint(*(center.getData()) - eps),
                                         createRPoint(*(center.getData()) + eps));
        index.query(bgi::intersects(query_box), selector);
    } else if(use_triangular_inequality_index) {
        T backward_threshold = center.getTriangularInequalityDistance() - eps;
        int index = center_index - 1;
        while(index >= 0) {
            if(points[index].getTriangularInequalityDistance() < backward_threshold) {
                break;
            }
            if(closeEnough(points[index], center)) {
                result.push_back(index);
            }
            index--;
        }
        if(use_cache) {
            cache.purgeSmaller(index + 1);
        }
        T forward_threshold = center.getTriangularInequalityDistance() + eps;
        for(int i = center_index + 1; i < (int) points.size(); i++) {
            if(points[i].getTriangularInequalityDistance() > forward_threshold) {
                break;
            }
            if(closeEnough(points[i], center)) {
                result.push_back(i);
            }
        }
    } else {
        for(int i = 0; i < (int) points.size(); i++) {
            if(i != center_index && closeEnough(points[i], center)) {
                result.push_back(i);
            }
        }
    }
    return result;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
void
DBSCAN<T, Dim, RStarMaxNodeElements>::expandCluster(vector<ClusterPoint<DataPoint<T, Dim>, T>> &points,
                                                    int center_index,
                                                    vector<int> &neighbours, int cluster_number)
{
    points[center_index].addAsCorePointToCluster(cluster_number);
    for(int i = 0; i < (int) neighbours.size(); i++) {
        auto n_i = neighbours[i];
        if(points[n_i].isNoise()) {
            points[n_i].addAsBorderPointToCluster(cluster_number);;
        }
        if((!border_points_in_many_clusters && points[n_i].isUndefined()) ||
           (border_points_in_many_clusters && !points[n_i].belongsToCluster(cluster_number))) {
            points[n_i].addAsBorderPointToCluster(cluster_number);
            auto neighbour_neighbours = regionQuery(points, n_i);
            if((int) neighbour_neighbours.size() >= min_neighbours - 1) {
                neighbours.reserve(neighbours.size() + neighbour_neighbours.size());
                neighbours.insert(neighbours.end(), neighbour_neighbours.begin(), neighbour_neighbours.end());
                points[n_i].setAsCorePoint();
            }
        }
    }
}


template<class T, int Dim, size_t RStarMaxNodeElements>
vector<ClusterPoint<DataPoint<T, Dim>, T>>
DBSCAN<T, Dim, RStarMaxNodeElements>::findClusters(const vector<pair<DataPoint<T, Dim>, vector<int>>> &data)
{
    auto result = buildIndex(data);
    auto start = high_resolution_clock::now();
    findClusters(result);
    auto stop = high_resolution_clock::now();
    auto clustering_duration = stop - start;
    int rand_tp, rand_tn, rand_all;
    auto rand_index = sortResultsAndCalculateRandIndex(data, result, rand_tp, rand_tn, rand_all);
    save(result, rand_index, rand_tp, rand_tn, rand_all, reading_input_file_time, clustering_duration);
    return result;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
double
DBSCAN<T, Dim, RStarMaxNodeElements>::sortResultsAndCalculateRandIndex(const vector<pair<DataPoint<T, Dim>, vector<int>>> &data,
                                                                       vector<ClusterPoint<DataPoint<T, Dim>, T>> &points,
                                                                       int &rand_tp, int &rand_tn, int &rand_all)
{
    for(const auto &zero : zeros) {
        points.emplace_back(-1, zero.first, zero.second);
    }
    sort(points.begin(), points.end(),
         [](const ClusterPoint<DataPoint<T, Dim>, T> &lhs, const ClusterPoint<DataPoint<T, Dim>, T> &rhs) -> bool {
             return lhs.getOriginalDataIndex() < rhs.getOriginalDataIndex();
         });

    int same_in_both = 0;
    int same_in_data_other_in_points = 0;
    int other_in_data_same_in_points = 0;
    int other_in_both = 0;

    for(unsigned int i = 0; i < points.size(); i++) {
        for(unsigned int j = i + 1; j < points.size(); j++) {
            bool same_in_data;
            if(border_points_in_many_clusters) {
                same_in_data = false;
                unsigned int k = 0;
                while(!same_in_data && k < data[i].second.size()) {
                    unsigned int l = 0;
                    const auto &d1 = data[i].second[k];
                    while(!same_in_data && l < data[j].second.size()) {
                        const auto &d2 = data[j].second[l];
                        if(d1 == d2) {
                            same_in_data = true;
                        }
                        l++;
                    }
                    k++;
                }
            } else {
                same_in_data = data[i].second[0] == data[j].second[0];
            }
            if(same_in_data) {
                if(points[i].assignedToTheSameCluster(points[j])) {
                    same_in_both++;
                } else {
                    same_in_data_other_in_points++;
                }
            } else {
                if(points[i].assignedToTheSameCluster(points[j])) {
                    other_in_data_same_in_points++;
                } else {
                    other_in_both++;
                }
            }
        }
    }

    rand_tp = same_in_both;
    rand_tn = other_in_both;
    rand_all = same_in_both + other_in_both + same_in_data_other_in_points + other_in_data_same_in_points;
    return (double(rand_tp + rand_tn)) / rand_all;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
void DBSCAN<T, Dim, RStarMaxNodeElements>::findClusters(vector<ClusterPoint<DataPoint<T, Dim>, T>> &data)
{
    if((int) data.size() < min_neighbours) {
        for(auto &point : data) {
            point.setNoise();
        }
        return;
    }

    int cluster_number = 0;
    for(int i = 0; i < (int) data.size(); i++) {
        if(data[i].isUndefined()) {
            auto neighbours = regionQuery(data, i);
            if((int) neighbours.size() <
               min_neighbours - 1) { //we don't have i-th point in the neighbours set, so subtract 1
                data[i].setNoise();
            } else {
                expandCluster(data, i, neighbours, cluster_number++);
            }
        }
    }
}

template<class T, int Dim, size_t RStarMaxNodeElements>
void DBSCAN<T, Dim, RStarMaxNodeElements>::useRStarTreeIndex()
{
    if(Dim != 2 && Dim != 3) {
        throw runtime_error("R*-tree index can be used only for data points with 2 or 3 dimensions.");
    }
    use_rstar_tree = true;
    use_triangular_inequality_index = false;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
void DBSCAN<T, Dim, RStarMaxNodeElements>::useTriangularInequalityIndex()
{
    use_rstar_tree = false;
    use_triangular_inequality_index = true;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
void DBSCAN<T, Dim, RStarMaxNodeElements>::useNoIndex()
{
    use_rstar_tree = false;
    use_triangular_inequality_index = false;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
bool DBSCAN<T, Dim, RStarMaxNodeElements>::closeEnough(ClusterPoint<DataPoint<T, Dim>, T> &p1,
                                                       ClusterPoint<DataPoint<T, Dim>, T> &p2)
{
    T distance;
    if(use_triangular_inequality_index && use_cache) {
        if(cache.get(p1.getDataIndex(), p2.getDataIndex(), distance)) {
            if(is_less_closer) {
                return (distance <= eps);
            } else {
                return (distance >= eps);
            }
        } else {
            p1.distanceCalculated();
            p2.distanceCalculated();
            auto start = high_resolution_clock::now();
            distance = distance_func(*(p1.getData()), *(p2.getData()));
            auto stop = high_resolution_clock::now();
            calc_dist_time += stop - start;
            cache.add(p1.getDataIndex(), p2.getDataIndex(), distance);
            if(is_less_closer) {
                return (distance <= eps);
            } else {
                return (distance >= eps);
            }
        }
    } else {
        p1.distanceCalculated();
        p2.distanceCalculated();
        auto start = high_resolution_clock::now();
        distance = distance_func(*(p1.getData()), *(p2.getData()));
        auto stop = high_resolution_clock::now();
        calc_dist_time += stop - start;
        if(is_less_closer) {
            return (distance <= eps);
        } else {
            return (distance >= eps);
        }
    }
}

template<class T, int Dim, size_t RStarMaxNodeElements>
void
DBSCAN<T, Dim, RStarMaxNodeElements>::save(const vector<ClusterPoint<DataPoint<T, Dim>, T>> &points, double rand_index,
                                           int rand_tp, int rand_tn, int rand_all,
                                           duration<double> reading_input_file_time, duration<double> clustering_time)
{
    auto start = high_resolution_clock::now();
    string method_name;
    if(use_triangular_inequality_index) {
        method_name += "TI-";
    }
    method_name += "DBSCAN";
    if(border_points_in_many_clusters) {
        method_name += "+";
    }

    ostringstream counter_oss;
    counter_oss << setfill('0') << setw(3) << save_counter++;

    stringstream eps_oss;
    eps_oss.setf(ios::fixed);
    eps_oss.precision(8);
    eps_oss << eps;

    string common_name_part = counter_oss.str() + "_" + method_name + "_" + dataset_name + "_D" +
                              to_string(Dim) + "_R" + to_string(points.size()) + "_m" +
                              to_string(min_neighbours) + "_e" + eps_oss.str();
    auto out_name = boost::filesystem::path(output_dir) / (string("out_") + common_name_part + ".csv");

    ofstream out(out_name);
    out << "\"point id\", ";
    for(int i = 1; i <= Dim; i++) {
        out << "\"d" << i << "\", ";
    }
    if(normalize_values || calculate_z_score) {
        for(int i = 1; i <= Dim; i++) {
            out << "\"d" << i << " normalized\", ";
        }
    }
    out << "\"distance calculations count\", \"point type\", \"cluster label\"" << endl;

    for(auto &point : points) {
        out << point;
    }
    auto stop = high_resolution_clock::now();
    duration<double> save_duration = stop - start;

    auto total_duration = save_duration + reading_input_file_time + clustering_time;
    auto stat_name = boost::filesystem::path(output_dir) / (string("stat_") + common_name_part + ".txt");
    ofstream stat(stat_name);

    int max_cluster = -1;
    int noise_points_count = 0;
    int core_points_count = 0;
    int border_points_count = 0;
    int distance_calc_count = 0;
    for(auto &point : points) {
        if(point.isNoise()) {
            noise_points_count++;
        } else {
            auto cluster_id = point.getLastClusterId();
            if(cluster_id > max_cluster) {
                max_cluster = cluster_id;
            }
            if(point.isCorePoint()) {
                core_points_count++;
            } else {
                border_points_count++;
            }
        }
        distance_calc_count += point.getDistanceCalculationCount();
    }

    distance_calc_count /= 2;

    stat << "algorithm used:\t\t\t\t\t\t\t\t\t" << method_name << endl
         << "with:" << endl
         << (use_rstar_tree ? " - R*-tree index\n" : "")
         << (use_triangular_inequality_index ? " - triangular inequality index\n" : "")
         << (!use_rstar_tree && !use_triangular_inequality_index ? " - no index\n" : "")
         << (calculate_z_score ? " - Z-score normalization\n" : "")
         << (normalize_values ? " - vector normalization\n" : "")
         << (treat_zero_as_noise ? " - zero vectors treated as noise\n" : "")
         << "name of the input file:\t\t\t\t\t\t\t\t" << input_file_name << endl
         << "# of dimensions of a point:\t\t\t\t\t\t\t" << Dim << endl
         << "# of points in the input file:\t\t\t\t\t\t\t" << points.size() << endl
         << "Eps:\t\t\t\t\t\t\t\t\t\t" << eps << endl
         << "minPts:\t\t\t\t\t\t\t\t\t\t" << min_neighbours << endl
         << "values of dimensions of a reference point:\t\t\t\t\t" << ti_center << endl
         << "time spent on:\t\t\t\t\t\t\t\t\t" << endl
         << " - reading the input file:\t\t\t\t\t\t\t" << reading_input_file_time.count() << " seconds" << endl
         << " - calculation of distances to a reference point:\t\t\t\t" << calc_dist_to_ti_center_time.count()
         << " seconds" << endl
         << " - normalization of vectors:\t\t\t\t\t\t\t" << calc_norm_time.count() << " seconds" << endl
         << " - sorting of points w.r.t. their distances to the reference point:\t\t" << sort_ti_time.count()
         << " seconds" << endl
         << " - calculation of distances between points:\t\t\t\t\t" << calc_dist_time.count() << " seconds" << endl
         << " - clustering:\t\t\t\t\t\t\t\t\t" << clustering_time.count() << " seconds" << endl
         << " - saving results:\t\t\t\t\t\t\t\t" << save_duration.count() << " seconds" << endl
         << " - total:\t\t\t\t\t\t\t\t\t" << total_duration.count() << " seconds" << endl
         << "# of discovered clusters:\t\t\t\t\t\t\t" << (max_cluster + 1) << endl
         << "# of discovered noise points:\t\t\t\t\t\t\t" << noise_points_count << endl
         << "# of discovered core points:\t\t\t\t\t\t\t" << core_points_count << endl
         << "# of discovered border points:\t\t\t\t\t\t\t" << border_points_count << endl
         << "Rand # pairs in same clusters:\t\t\t\t\t\t\t" << rand_tp << endl
         << "Rand # pairs in other clusters:\t\t\t\t\t\t\t" << rand_tn << endl
         << "Rand # all pairs:\t\t\t\t\t\t\t\t" << rand_all << endl
         << "Rand index:\t\t\t\t\t\t\t\t\t" << rand_index << endl
         << "avg # of calculations of distance/similarity of a point to other points:\t"
         << (points.empty() ? 0 : ((double) distance_calc_count + dist_calc_for_index_counter) / points.size()) << endl;
    if(use_triangular_inequality_index && use_cache) {
        stat << "cache:" << endl
             << " - time used by cache:\t\t\t\t\t\t\t\t" << cache.time().count() << " seconds" << endl
             << " - misses:\t\t\t\t\t\t\t\t\t" << cache.misses() << endl
             << " - hits:\t\t\t\t\t\t\t\t\t" << cache.hits() << endl
             << " - maximum number of elements stored in cache:\t\t\t\t\t" << cache.maxSize() << endl;
    }

}

template<class T, int Dim, size_t RStarMaxNodeElements>
bool DBSCAN<T, Dim, RStarMaxNodeElements>::readData(const string &file_name, const string &dataset_name,
                                                    vector<pair<DataPoint<T, Dim>, vector<int>>> &data)
{
    auto start = high_resolution_clock::now();
    this->dataset_name = dataset_name;
    ifstream str(file_name);
    if(!str.is_open()) {
        return false;
    }
    input_file_name = file_name;
    string line = "";

    if(!getline(str, line)) {
        cerr << "Malformed input file - no header row." << endl;
        return false;
    }
    //we expect Dim comas because in the input file we have always extra ground truth column
    if(count(line.begin(), line.end(), ',') != Dim) {
        cerr << "Malformed input file. Header contains unexpected number of columns." << endl;
        return false;
    }

    while(getline(str, line)) {
        stringstream line_stream(line);
        string elem = "";
        vector<T> vals;
        vals.reserve(Dim - 1);
        int val_index = 0;
        vector<int> ground_truth;
        while(getline(line_stream, elem, ',')) {
            if(val_index < Dim) {
                T val;
                try {
                    if(is_same<T, float>::value) {
                        val = stof(elem);
                    } else if(is_same<T, double>::value) {
                        val = stod(elem);
                    } else if(is_same<T, long double>::value) {
                        val = stold(elem);
                    } else if(is_same<T, int>::value) {
                        val = stoi(elem);
                    } else if(is_same<T, long>::value) {
                        val = stol(elem);
                    } else if(is_same<T, unsigned long>::value) {
                        val = stoul(elem);
                    } else if(is_same<T, long long>::value) {
                        val = stoll(elem);
                    } else if(is_same<T, unsigned long long>::value) {
                        val = stoull(elem);
                    } else {
                        return false;
                    }
                } catch(const exception &e) {
                    cerr << "Malformed input file." << endl;
                    return false;
                }
                vals.push_back(val);
            } else if(val_index >= Dim) {
                ground_truth.push_back(stoi(elem));
            }
            val_index++;
        }
        if(vals.empty()) {
            //just skip the line
        }
        if(ground_truth.empty()) {
            cerr << "Malformed input file. No ground truth values." << endl;
            return false;
        }
        data.emplace_back(DataPoint<T, Dim>(vals), ground_truth);
    }
    auto stop = high_resolution_clock::now();
    reading_input_file_time = stop - start;
    return true;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
void DBSCAN<T, Dim, RStarMaxNodeElements>::setTriangularInequalityCenter(const DataPoint<T, Dim> &point)
{
    ti_center = point;
}


template<class T, int Dim, size_t RStarMaxNodeElements>
void DBSCAN<T, Dim, RStarMaxNodeElements>::useDataNormalization()
{
    normalize_values = true;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
void DBSCAN<T, Dim, RStarMaxNodeElements>::useNoDataNormalization()
{
    normalize_values = false;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
void DBSCAN<T, Dim, RStarMaxNodeElements>::useZScoreNormalization()
{
    calculate_z_score = true;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
void DBSCAN<T, Dim, RStarMaxNodeElements>::useNoZScoreNormalization()
{
    calculate_z_score = false;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
void DBSCAN<T, Dim, RStarMaxNodeElements>::treatZeroVectorsAsNoise()
{
    treat_zero_as_noise = true;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
void DBSCAN<T, Dim, RStarMaxNodeElements>::treatZeroVectorsAsNotNoise()
{
    treat_zero_as_noise = false;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
void DBSCAN<T, Dim, RStarMaxNodeElements>::assignBorderPointsToOnlyOneCluster()
{
    border_points_in_many_clusters = false;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
void DBSCAN<T, Dim, RStarMaxNodeElements>::assignBorderPointsToPossiblyManyClusters()
{
    border_points_in_many_clusters = true;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
void DBSCAN<T, Dim, RStarMaxNodeElements>::setEpsilon(T new_epsilon)
{
    guess_epsilon = false;
    eps = new_epsilon;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
void DBSCAN<T, Dim, RStarMaxNodeElements>::setMinNeighbours(int new_min_neighbours)
{
    min_neighbours = new_min_neighbours;
}


template<class T, int Dim, size_t RStarMaxNodeElements>
void DBSCAN<T, Dim, RStarMaxNodeElements>::useCache()
{
    use_cache = true;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
void DBSCAN<T, Dim, RStarMaxNodeElements>::useNoCache()
{
    use_cache = false;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
void DBSCAN<T, Dim, RStarMaxNodeElements>::setOutputFolder(const string &folder)
{
    output_dir = folder;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
void DBSCAN<T, Dim, RStarMaxNodeElements>::setGuessEpsilon()
{
    guess_epsilon = true;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
void DBSCAN<T, Dim, RStarMaxNodeElements>::useDistanceFunction(
        function<T(const DataPoint<T, Dim> &, const DataPoint<T, Dim> &)> distance)
{
    distance_func = distance;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
void DBSCAN<T, Dim, RStarMaxNodeElements>::setSmallerIsCloser()
{
    is_less_closer = true;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
void DBSCAN<T, Dim, RStarMaxNodeElements>::setBiggerIsCloser()
{
    is_less_closer = false;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
bool DBSCAN<T, Dim, RStarMaxNodeElements>::guessEpsilon(vector<ClusterPoint<DataPoint<T, Dim>, T>> &data)
{
    const int section_count = 5;
    const int section_width = 10;
    if(data.size() < section_count * section_width) {
        return false;
    }
    //for now just the simplest implementation... it probably can be done faster.
    vector<T> distances(data.size());
    const int number_of_neighbours = Dim * 2;
    for(int i = 0; i < (int) data.size(); i++) {
        distances[i] = calcDistanceToKthClosest(data, i, number_of_neighbours);
    }
    if(is_less_closer) {
        sort(distances.begin(), distances.end(), [](const T &lhs, const T &rhs) -> bool {
            return lhs < rhs;
        });
    } else {
        sort(distances.begin(), distances.end(), [](const T &lhs, const T &rhs) -> bool {
            return lhs > rhs;
        });
    }
    T scale = (T) (distances.size()) / (distances[distances.size() - 1] - distances[0]);
    T first = distances[0];
    for(int i = 0; i < (int) distances.size(); i++) {
        distances[i] = (distances[i] - first) * scale;
    }
    const T threshold = (T) distances.size() / 10;
    int i = 0;
    while(i + section_width < (int) distances.size() && distances[i + section_width] - distances[i] < threshold) {
        i++;
    }
    eps = distances[i] / scale + first;
    min_neighbours = number_of_neighbours;
    ofstream out("distances.csv");
    for(auto &distance : distances) {
        out << distance << endl;
    }
    return true;
}

template<class T, int Dim, size_t RStarMaxNodeElements>
T DBSCAN<T, Dim, RStarMaxNodeElements>::calcDistanceToKthClosest(vector<ClusterPoint<DataPoint<T, Dim>, T>> &data,
                                                                 int index, int k)
{
    auto smallest_comparer = [this](const T &lhs, const T &rhs) -> bool {
        if(this->is_less_closer) {
            return lhs < rhs;
        } else {
            return lhs > rhs;
        }
    };
    priority_queue<T, vector<T>, decltype(smallest_comparer)> smallest(smallest_comparer);
    const auto &selected = data[index];
    for(int i = 0; i < (int) data.size(); i++) {
        if(i != index) {
            auto distance = distance_func(*(selected.getData()), *(data[i].getData()));
            smallest.push(distance);
            if((int) smallest.size() > k) {
                smallest.pop();
            }
        }
    }
    if(smallest.empty()) {
        return 0;
    } else {
        return smallest.top();
    }
}

#endif //DBSCAN_DBSCAN_HPP
