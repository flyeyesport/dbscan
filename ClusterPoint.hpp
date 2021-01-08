//
// Created by jarek on 19.12.20.
//

#ifndef DBSCAN_CLUSTERPOINT_HPP
#define DBSCAN_CLUSTERPOINT_HPP

#include <vector>

using namespace std;

/**
 * @brief Type of the point recognised by the clustering algorithm. It can be:
 * - a CORE point, a point with enough neighbours in neighbour area of radius <b>eps</b> - it belongs to some cluster,
 * - a BORDER point, note a CORE point but its neighbour which still belongs to som cluster (or possibly many clusters
 *   for DBSCAN+ version of the algorithm),
 * - a NOISE point which does not belong to any cluster.
 */
enum class PointType {
    CORE = 1,
    BORDER = 0,
    NOISE = -1
};

/**
 * @brief Special identifier of the cluster (rather group) of points that were identified by the clustering algorithm
 * as noise. Do not confuse with PointType::NOISE. It must be negative, because it is compared with cluster identifiers
 * which are greater or equal to 0.
 */
const int NOISE = -1;

/**
 * @brief Special identifier of the cluster (rather group) of points that were not yet processed by the clustering
 * algorithm. It must be negative, because it is compared with cluster identifiers which are greater or equal to 0.
 */
const int UNDEFINED = -2;

/**
 * @brief A class representing an extension to the data point allowing the clustering algorithm to operate on data
 * points without modifying original data points.
 * @tparam T Type of the original data point.
 * @tparam U Arithmetic type used by the original data point in all of its dimensions.
 */
template<class T, class U>
class ClusterPoint {
protected:

    /**
     * @brief An index of the point in a set of selected, extended data points. Points of this class are stored in a
     * structure with strict ordering. For efficiency reasons each such a point knows its position it this structure.
     */
    int index;

    /**
     * @brief An index of the original data point on which this point is constructed as an extension. This information
     * is used only to print data or may be used to sort extended points in an order of original data points before
     * returning results.
     */
    int original_index;

    /**
     * @brief A constant pointer to original data point. If data is not normalized it allows to save space and avoid
     * unnecessary copying of data from original points to extended ones.
     */
    const T *data_point;

    /**
     * @brief Internal counter counting all comparisons of distances with other points. Note that every time a distance
     * is calculated between 2 points, both of them increase their counters.
     */
    int dist_calc_counter;

    /**
     * @brief Type of the point calculated by the clustering algorithm. @see <b>PointType</b> type for more info.
     */
    PointType point_type;

    /**
     * @brief Internal list of identifiers/labels of clusters to which a point was assigned by the clustering algorithm.
     * Note that it is a list instead of single value, because border points in DBSCAN+ version of the algorithm can be
     * assigned to more than one cluster.
     */
    vector<int> cluster_ids;

    /**
     * @brief Boolean flag indicating if data vector was normalized (length of vector set to 1). If this flag is set
     * to true, then other normalizing methods are switched off.
     */
    bool is_normalized;

    /**
     * @brief Boolean flag indicating if data vector was normalized using z-score. If this flag is set
     * to true, then other normalizing methods are switched off.
     */
    bool is_z_score_calculated;

    /**
     * @brief Pointer to a copy of original data point that was normalized with one of the available methods.
     * If the normalization is switched off, then this pointer is set to nullptr, and original data is used.
     */
    T *modified_data_point;

    /**
     * @brief Calculated distance to some center point used by triangular inequality index. It is only used when
     * triangular inequality index is used by the clustering algorithm.
     */
    U ti_distance;

    /**
     * @brief Internal method assigns this point to a cluster with a given number/label.
     * @param cluster_id A label of the cluster to which a point should be assigned.
     */
    void addPointToCluster(int cluster_id);

public:
    /**
     * @brief Constructor.
     * @param index Index of this point in a strictly sorted structure of all such points.
     * @param original_index Index of the original data point which this point is extending.
     * @param data_point Pointer to the original data point. It is used to access data when it is not modified by
     * normalization.
     */
    ClusterPoint(int index, int original_index, const T *data_point);

    /**
     * @brief Copy constructor. Note that <b>modified_data_point</b> must be treated specially.
     * @param other Data to be copied.
     */
    ClusterPoint(const ClusterPoint<T, U> &other);

    /**
     * @brief Move constructor. Note that <b>modified_data_point</b> must be treated specially.
     * @param other Data to be moved.
     */
    ClusterPoint(ClusterPoint<T, U> &&other);

    /**
     * @brief Copy operator. Note that <b>modified_data_point</b> must be treated specially.
     * @param other Data to be copied.
     * @return This object with data copied from <b>other</b>.
     */
    ClusterPoint<T, U> &operator=(const ClusterPoint<T, U> &other);

    /**
     * @brief Move operator. Note that <b>modified_data_point</b> must be treated specially.
     * @param other Data to be moved.
     * @return This object with data copied from <b>other</b>.
     */
    ClusterPoint<T, U> &operator=(ClusterPoint<T, U> &&other);

    /**
     * @brief Destructor. Frees <b>modified_data_point</b>, which could be allocated.
     */
    ~ClusterPoint();

    /**
     * @brief Marks the point as noise.
     */
    void setNoise();

    /**
     * @brief Tests if a noise was added to the group/cluster of noise points.
     * @return True if the point was added to the group of noise points, false otherwise.
     */
    bool isNoise() const;

    /**
     * @brief Tests if the point was already processed by the algorithm or not.
     * @return True if the point was already processed by the algorithm, false otherwise.
     */
    bool isUndefined() const;

    /**
     * @brief Tests if the point was identified by the clustering algorithm as a CORE point.
     * @return True if the point is a CORE point, false otherwise.
     */
    bool isCorePoint() const;

    /**
     * @brief Marks the point as CORE point and adds to the cluster with a given label.
     * @param cluster_id A label of the cluster to which a point should be added.
     */
    void addAsCorePointToCluster(int cluster_id);

    /**
     * @brief Marks the point as BORDER point and adds to the cluster with a given label.
     * @param cluster_id A label of the cluster to which a point should be added.
     */
    void addAsBorderPointToCluster(int cluster_id);

    /**
     * @brief Marks the point as CORE point.
     */
    void setAsCorePoint();

    /**
     * @brief Getter for <b>index</b> value.
     * @return Index of this object in a strictly ordered structure where all extended points are stored.
     */
    int getDataIndex() const;

    /**
     * @brief Getter for <b>original_index</b> value.
     * @return Index of the original point in the set of data points which is a strongly ordered structure and keeps
     * all original data points.
     */
    int getOriginalDataIndex() const;

    /**
     * @brief Setter for <b>index</b>. It is used only if the order of the points in the strictly ordered set of
     * points has changed.
     * @param data_index New value of the index to be set.
     */
    void setDataIndex(int data_index);

    /**
     * @brief Getter for original data point.
     * @return Constant pointer to the original data point.
     */
    const T *getData() const;

    /**
     * @brief Checks if the point was assigned to the cluster with a given label.
     * @param cluster_id A label of the cluster.
     * @return True if point was assigned to the given cluster, false otherwise.
     */
    bool belongsToCluster(int cluster_id) const;

    /**
     * @brief Increases internal counter of times that the distance between this point and some other point was
     * calculated. It should be called every time a distance is calculated.
     */
    void distanceCalculated();

    /**
     * @brief Getter for counter of distance calculations.
     * @return Number of distance calculations.
     */
    int getDistanceCalculationCount() const;

    /**
     * @brief Getter for the cluster label that was assigned to this point. In case of BORDER points and DBSCAN+
     * version of the clustering algorithm the last assigned cluster label is returned.
     * @return The cluster label assigned to the point or NOISE or UNDEFINED.
     */
    int getLastClusterId() const;

    /**
     * @brief Normalizes data vector. Internally it creates a copy of the original data and normalizes this copy.
     */
    void normalize();

    /**
     * @brief Normalizes data vector using z-score. Internally it creates a copy of the original data and normalizes
     * this copy.
     */
    void calculateZScore(const T &avg, const T &std);

    /**
     * @brief Restores data to its original values.
     */
    void removeNormalization();

    /**
     * @brief Setter for <b>ti_distance</b>.
     * @param distance A new value to be set.
     */
    void setTriangularInequalityDistance(U distance);

    /**
     * @brief Setter for <b>ti_distance</b>.
     * @param distance A new value to be set.
     */
    U getTriangularInequalityDistance() const;

    /**
     * @brief An operator used to print an extended point and its original data point.
     * @param str Stream to print to.
     * @param obj An object/point from which to take data to print.
     * @return a reference to a modified stream to which data was printed.
     */
    friend ostream &operator<<(ostream &str, const ClusterPoint<T, U> &obj)
    {
        str << obj.original_index << ", " << *(obj.data_point) << ", ";
        if(obj.modified_data_point != nullptr) {
            str << *(obj.modified_data_point) << ", ";
        }

        str << obj.dist_calc_counter << ", ";
        switch(obj.point_type) {
            case PointType::NOISE:
                str << "-1";
                break;
            case PointType::BORDER:
                str << "0";
                break;
            case PointType::CORE:
                str << "1";
                break;
        }

        if(obj.isNoise()) {
            str << ", -1" << endl;
        } else {
            for(auto cluster_id : obj.cluster_ids) {
                str << ", " << cluster_id;
            }
            str << endl;
        }
        return str;
    }

};

//------------------------------------- IMPLEMENTATION ----------------------------------------

template<class T, class U>
void ClusterPoint<T, U>::addPointToCluster(int cluster_id)
{
    if(cluster_ids.back() == NOISE || cluster_ids.back() == UNDEFINED) {
        cluster_ids.back() = cluster_id;
    } else {
        cluster_ids.push_back(cluster_id);
    }
}

template<class T, class U>
ClusterPoint<T, U>::ClusterPoint(int index, int original_index, const T *data_point)
        : index(index), original_index(original_index), data_point(data_point), dist_calc_counter(0),
          point_type(PointType::NOISE), is_normalized(false), is_z_score_calculated(false),
          modified_data_point(nullptr), ti_distance(0)
{
    cluster_ids.push_back(UNDEFINED);
}

template<class T, class U>
ClusterPoint<T, U>::ClusterPoint(const ClusterPoint<T, U> &other)
        : index(other.index), original_index(other.original_index), data_point(other.data_point), dist_calc_counter(other.dist_calc_counter),
          point_type(other.point_type), cluster_ids(other.cluster_ids), is_normalized(other.is_normalized),
          is_z_score_calculated(other.is_z_score_calculated),
          modified_data_point(other.modified_data_point == nullptr ? nullptr : new T(*(other.modified_data_point))),
          ti_distance(other.ti_distance)
{}

template<class T, class U>
ClusterPoint<T, U>::ClusterPoint(ClusterPoint<T, U> &&other)
        : index(other.index), original_index(other.original_index), data_point(other.data_point),
          dist_calc_counter(other.dist_calc_counter),
          point_type(other.point_type), cluster_ids(move(other.cluster_ids)),
          is_normalized(other.is_normalized),
          is_z_score_calculated(other.is_z_score_calculated),
          modified_data_point(exchange(other.modified_data_point, nullptr)),
          ti_distance(other.ti_distance)
{}


template<class T, class U>
ClusterPoint<T, U> &ClusterPoint<T, U>::operator=(const ClusterPoint<T, U> &other)
{
    if(this != &other) {
        delete modified_data_point;
        index = other.index;
        original_index = other.original_index;
        data_point = other.data_point;
        dist_calc_counter = other.dist_calc_counter;
        point_type = other.point_type;
        cluster_ids = other.cluster_ids;
        is_normalized = other.is_normalized;
        is_z_score_calculated = other.is_z_score_calculated;
        if(other.modified_data_point == nullptr) {
            modified_data_point = nullptr;
        } else {
            modified_data_point = new T(*(other.modified_data_point));
        }
        ti_distance = other.ti_distance;
    }
    return *this;
}

template<class T, class U>
ClusterPoint<T, U> &ClusterPoint<T, U>::operator=(ClusterPoint<T, U> &&other)
{
    if(this != &other) {
        delete modified_data_point;
        index = other.index;
        original_index = other.original_index;
        data_point = other.data_point;
        dist_calc_counter = other.dist_calc_counter;
        point_type = other.point_type;
        cluster_ids = move(other.cluster_ids);
        is_normalized = other.is_normalized;
        is_z_score_calculated = other.is_z_score_calculated;
        modified_data_point = exchange(other.modified_data_point, nullptr);
        ti_distance = other.ti_distance;
    }
    return *this;
}

template<class T, class U>
ClusterPoint<T, U>::~ClusterPoint()
{
    delete modified_data_point;
}

template<class T, class U>
void ClusterPoint<T, U>::setNoise()
{
    cluster_ids.clear();
    cluster_ids.push_back(NOISE);
    point_type = PointType::NOISE;
}

template<class T, class U>
bool ClusterPoint<T, U>::isNoise() const
{
    return (cluster_ids.back() == NOISE);
}

template<class T, class U>
bool ClusterPoint<T, U>::isUndefined() const
{
    return (cluster_ids.back() == UNDEFINED);
}

template<class T, class U>
void ClusterPoint<T, U>::addAsCorePointToCluster(int cluster_id)
{
    point_type = PointType::CORE;
    addPointToCluster(cluster_id);
}

template<class T, class U>
void ClusterPoint<T, U>::addAsBorderPointToCluster(int cluster_id)
{
    point_type = PointType::BORDER;
    addPointToCluster(cluster_id);
}

template<class T, class U>
void ClusterPoint<T, U>::setAsCorePoint()
{
    point_type = PointType::CORE;
}

template<class T, class U>
int ClusterPoint<T, U>::getDataIndex() const
{
    return index;
}

template<class T, class U>
int ClusterPoint<T, U>::getOriginalDataIndex() const
{
    return original_index;
}

template<class T, class U>
void ClusterPoint<T, U>::setDataIndex(int data_index)
{
    index = data_index;
}

template<class T, class U>
const T *ClusterPoint<T, U>::getData() const
{
    if(is_normalized || is_z_score_calculated) {
        return modified_data_point;
    } else {
        return data_point;
    }
}

template<class T, class U>
bool ClusterPoint<T, U>::belongsToCluster(int cluster_id) const
{
    for(int i = ((int) cluster_ids.size()) - 1; i >= 0; i--) {
        if(cluster_ids[i] == cluster_id) {
            return true;
        }
    }
    return false;
}

template<class T, class U>
void ClusterPoint<T, U>::distanceCalculated()
{
    dist_calc_counter++;
}

template<class T, class U>
int ClusterPoint<T, U>::getLastClusterId() const
{
    return cluster_ids.back();
}

template<class T, class U>
bool ClusterPoint<T, U>::isCorePoint() const
{
    return (point_type == PointType::CORE);
}

template<class T, class U>
int ClusterPoint<T, U>::getDistanceCalculationCount() const
{
    return dist_calc_counter;
}

template<class T, class U>
void ClusterPoint<T, U>::calculateZScore(const T &avg, const T &std)
{
    if(is_z_score_calculated) {
        return;
    }
    is_z_score_calculated = true;
    delete modified_data_point;
    modified_data_point = new T(*data_point);
    modified_data_point->calculateZScore(avg, std);
}

template<class T, class U>
void ClusterPoint<T, U>::normalize()
{
    if(is_normalized) {
        return;
    }
    is_normalized = true;
    delete modified_data_point;
    modified_data_point = new T(*data_point);
    modified_data_point->normalize();
}

template<class T, class U>
void ClusterPoint<T, U>::removeNormalization()
{
    if(!is_normalized && !is_z_score_calculated) {
        return;
    }
    is_normalized = false;
    is_z_score_calculated = false;
    delete modified_data_point;
    modified_data_point = nullptr;
}

template<class T, class U>
void ClusterPoint<T, U>::setTriangularInequalityDistance(U distance)
{
    ti_distance = distance;
}

template<class T, class U>
U ClusterPoint<T, U>::getTriangularInequalityDistance() const
{
    return ti_distance;
}

#endif //DBSCAN_CLUSTERPOINT_HPP
