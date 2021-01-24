//
// Created by jarek on 18.12.20.
//

#ifndef DBSCAN_DATAPOINT_HPP
#define DBSCAN_DATAPOINT_HPP

#include <vector>
#include <stdlib.h>
#include <cmath>

using namespace std;

/**
 * @brief A simple structure used in templates to verify if <b>Dim</b> template param has some given value.
 * @tparam Dummy Not used, but needed for templates in our case to compile.
 * @tparam N A given value to which we compare <b>Dim</b> template param.
 * @tparam Dim A template parameter which value we verify and compare to N.
 */
template<typename Dummy, int N, int Dim>
struct is_equal_to_dim
{
    static const bool value = (N == Dim);
};

/**
 * @brief Class representing a data point with variable, but known at compile time number of dimensions. Each dimension
 * has the same type.
 * @tparam T Type of each dimension.
 * @tparam Dim Number of dimensions.
 */
template <class T, int Dim, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
class DataPoint {
protected:

    /**
     * @brief list of values of all dimensions.
     */
    vector<T> data;

public:

    /**
     * @brief Empty, dummy constructor. Used for example by stl containers when initializing with sizes greater than 0
     * with empty elements to be replaced by some meaningful data in the future.
     */
    DataPoint() : data(vector<T>(Dim, 0)) {}

    /**
     * @brief Constructor. Copies data from vector to our internal structure.
     * @param data A list of values for every dimension.
     */
    explicit DataPoint(const vector<T> &data) : data(data) {
        if(data.size() != Dim) {
            throw runtime_error("Number of dimensions passed to constructor of DataPoint is different than expected.");
        }
    }

    DataPoint(T val) : data(vector<T>(Dim, 0)) {
        data[0] = val;
    }

    /**
     * @brief A version of the constructor for data point with only one dimension.
     * @tparam C Type of the data.
     * @param x Value of the data.
     */
    template<typename C = T, enable_if_t<is_equal_to_dim<C, 1, Dim>::value, int> = 0>
    explicit DataPoint(T x) : data({x}) {}

    /**
     * @brief A version of the constructor for data point with exactly two dimensions.
     * @tparam C Type of the data in each dimension.
     * @param x Value of the data in the first dimension.
     * @param y Value of the data in the second dimension.
     */
    template<typename C = T, enable_if_t<is_equal_to_dim<C, 2, Dim>::value, int> = 0>
    DataPoint(T x, T y) : data({x, y}) {}

    /**
     * @brief A version of the constructor for data point with exactly three dimensions.
     * @tparam C Type of the data in each dimension.
     * @param x Value of the data in the first dimension.
     * @param y Value of the data in the second dimension.
     * @param z Value of the data in the third dimension.
     */
    template<typename C = T, enable_if_t<is_equal_to_dim<C, 3, Dim>::value, int> = 0>
    DataPoint(T x, T y, T z) : data({x, y, z}) {}

    /**
     * @brief An operator used to print data point.
     * @param str Stream to print to.
     * @param obj An object/point from which to take data to print.
     * @return a reference to a modified stream to which data was printed.
     */
    friend ostream& operator<<(ostream &str, const DataPoint<T, Dim> &obj)
    {
        bool first = true;
        for(auto &val : obj.data) {
            if(first) {
                first = false;
            } else {
                str << ", ";
            }
            str << val;
        }
        return str;
    }

    /**
     * @brief Access operator for value of selected dimension.
     * @param i The number of dimension for which data should be returned.
     * @return Reference to the value of the given dimension.
     */
    const T &operator[](int i) const
    {
        return data[i];
    }

    /**
     * @brief Tests if the data point has all dimensions set to value 0.
     * @return True if all dimensions are equal to zero, false otherwise.
     */
    bool isZero() const
    {
        for(auto &val : data) {
            if(val != 0) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Sets values of all dimensions to 0.
     */
    void setZero()
    {
        for(auto &val : data) {
            val = 0;
        }
    }

    /**
     * @brief Sets values of all dimensions it their absolute values.
     */
    void abs()
    {
        for(auto &val : data) {
            val = std::abs(val);
        }
    }

    /**
     * @brief Subtract operator. It does not modify this point.
     * @param val A value to be subtracted from all dimensions.
     * @return A new point with all values reduced by <b>val</b>.
     */
    DataPoint<T, Dim> operator-(const T &val) const
    {
        DataPoint<T, Dim> result(*this);
        for(auto &elem : result.data) {
            elem -= val;
        }
        return result;
    }

    /**
     * @brief Add operator. It does not modify this point.
     * @param val A value to be added to all dimensions.
     * @return A new point with all values increased by <b>val</b>.
     */
    DataPoint<T, Dim> operator+(const T &val) const
    {
        DataPoint<T, Dim> result(*this);
        for(auto &elem : result.data) {
            elem += val;
        }
        return result;
    }

    /**
     * @brief Division operator. It does not modify this point.
     * @param val A value by which all dimensions should be divided.
     * @return A new point with all values divided by <b>val</b>.
     */
    DataPoint<T, Dim> operator/(const T &val) const
    {
        DataPoint<T, Dim> result(*this);
        for(auto &elem : result.data) {
            elem /= val;
        }
        return result;
    }

    /**
     * @brief Subtract operator. It does not modify this point.
     * @param other A point to be subtracted from this point.
     * @return A new point with all values set to results of subtraction for each dimension.
     */
    DataPoint<T, Dim> operator-(const DataPoint<T, Dim> &other) const
    {
        DataPoint<T, Dim> result(*this);
        int i = 0;
        for(auto &elem : result.data) {
            elem -= other.data[i];
            i++;
        }
        return result;
    }

    /**
     * @brief Add operator. It does not modify this point.
     * @param other A point to be added to this point.
     * @return A new point with all values set to results of addition for each dimension.
     */
    DataPoint<T, Dim> operator+(const DataPoint<T, Dim> &other) const
    {
        DataPoint<T, Dim> result(*this);
        int i = 0;
        for(auto &elem : result.data) {
            elem += other.data[i];
            i++;
        }
        return result;
    }

    /**
     * @brief Normalizes each dimension by calculating z-score. Modifies this point.
     * @param avg Average used in z-score formula.
     * @param std Standard deviation used in z-score formula.
     */
    void calculateZScore(const DataPoint<T, Dim> &avg, const DataPoint<T, Dim> &std)
    {
        int i = 0;
        for(auto &elem : data) {
            elem = (elem - avg.data[i]) / std.data[i];
            i++;
        }
    }

    /**
     * @brief Normalizes each dimension. Modifies this point.
     */
    void normalize()
    {
        T sum = 0;

        for(auto elem : data) {
            sum += elem * elem;
        }

        T denominator = sqrt(sum);

        if(denominator == 0) {
            throw logic_error("The point to be normalized is a zero vector (or very close to it).");
        }

        for(auto &elem : data) {
            elem /= denominator;
        }
    }


};


#endif //DBSCAN_DATAPOINT_HPP
