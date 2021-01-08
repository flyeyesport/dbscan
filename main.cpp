#include <iostream>
#include <vector>
#include <boost/program_options.hpp>
#include "DBSCAN.hpp"

using namespace std;
using namespace boost::program_options;

bool detectDataDimension(const string &input_file, int &dim)
{
    ifstream str(input_file);
    if(!str.is_open()) {
        return false;
    }
    string line = "";
    if(!getline(str, line)) {
        return false;
    }
    dim = count(line.begin(), line.end(), ',');
    return true;
}

int main(int argc, char **argv)
{
    auto euclidean_distance = [] < typename T,
    int Dim>(const DataPoint<T, Dim> &point_1,
    const DataPoint<T, Dim> &point_2) -> T{
            T sum = 0;
            for(int i = 0; i < Dim; i++) {
                auto d = point_1[i] - point_2[i];
                sum += d * d;
            }
            return sqrt(sum);
    };

    auto cosine_distance = [] < typename T,
    int Dim>(const DataPoint<T, Dim> &point_1,
    const DataPoint<T, Dim> &point_2) -> T{
            T enumerator = 0;
            for(int j = 0; j < Dim; j++) {
                enumerator += point_1[j] * point_2[j];
            }

            T sum_1 = 0;
            for(int j = 0; j < Dim; j++) {
                sum_1 += point_1[j] * point_1[j];
            }
            T sum_2 = 0;
            for(int j = 0; j < Dim; j++) {
                sum_2 += point_2[j] * point_2[j];
            }

            T denominator = sqrt(sum_1) * sqrt(sum_2);
            return enumerator / denominator;
    };


    options_description desc("This is a simple program to test DBSCAN algorithm. Usage");
    desc.add_options()
            ("help,h", "Show help and exit.")
            ("input-file,i", value < string > ()->required(), "Input file with data to be processed.")
            ("output-dir,o", value < string > ()->required(),
             "Directory to which 2 files with results and statistics will be written.")
            ("name,n", value < string > ()->required(),
             "Name of the dataset. It is used as part of the names of output files.")
            ("eps,e", value < double > (),
             "Radius of the neighbourhoods of the core points in which the algorithm searches for cluster member"
             "candidates.")
            ("guess-eps", bool_switch(),
             "Algorithm tries to guess the radius of the neighbourhoods of the core points in which the algorithm "
             "searches for cluster member candidates.")
            ("min-neighbours,m", value < int > (),
             "Minimum number of neighbours in an area with radius <eps> to classify a point as core point.")
            ("use-euclidean-distance", bool_switch(), "Use Euclidean distance to calculate distance between data "
                                                      "points.")
            ("use-cosine-distance", bool_switch(), "Use cosine distance to calculate distance between data points.")
            ("bigger-is-closer", bool_switch(), "Use this switch if you want to invert the meaning of the used "
                                                "distance function. For example cosine distance returns bigger values "
                                                "for closer points.")
             ("use-norm", bool_switch(), "Normalize data points before clustering.")
            ("use-z-score", bool_switch(), "Normalize data points using z-score before clustering.")
            ("use-plus", bool_switch(), "Use DBSCAN+ version of the algorithm instead of classical DBSCAN.")
            ("use-ti", bool_switch(), "Use triangular inequality to speed up processing.")
            ("use-rstar", bool_switch(),
             "Use R*-tree index to speed up processing of two- and three-dimensional data.")
            ("use-cache", bool_switch(), "Use additional cache for triangular inequality index.");


    variables_map vm;

    bool opt_guess_eps, opt_use_euclidean_distance, opt_use_cosine_distance, opt_use_norm, opt_use_z_score, opt_use_plus,
            opt_use_ti, opt_use_rstar, opt_use_cache, opt_bigger_is_closer;

    try {
        store(parse_command_line(argc, argv, desc), vm);
        if(vm.count("help")) {
            cout << desc << endl;
            return 1;
        }
        opt_guess_eps = vm["guess-eps"].as<bool>();
        opt_use_euclidean_distance = vm["use-euclidean-distance"].as<bool>();
        opt_use_cosine_distance = vm["use-cosine-distance"].as<bool>();
        opt_use_norm = vm["use-norm"].as<bool>();
        opt_use_z_score = vm["use-z-score"].as<bool>();
        opt_use_plus = vm["use-plus"].as<bool>();
        opt_use_ti = vm["use-ti"].as<bool>();
        opt_use_rstar = vm["use-rstar"].as<bool>();
        opt_use_cache = vm["use-cache"].as<bool>();
        opt_bigger_is_closer = vm["bigger-is-closer"].as<bool>();

        notify(vm);

        if(vm.count("eps") == 0 && !opt_guess_eps) {
            throw runtime_error("exactly one option '--eps' or '--guess-eps' is required but none was used");
        } else if(vm.count("min-neighbours") == 0 && !opt_guess_eps) {
            throw runtime_error("exactly one option '--min-neighbours' or '--guess-eps' is required but none was used");
        } else if(vm.count("eps") == 1 && opt_guess_eps) {
            throw runtime_error("exactly one option '--eps' or '--guess-eps' is required but both were used");
        } else if(vm.count("min-neighbours") == 1 && opt_guess_eps) {
            throw runtime_error("exactly one option '--min-neighbours' or '--guess-eps' is required but both were used");
        }
        if(!opt_use_euclidean_distance && !opt_use_cosine_distance) {
            throw runtime_error("exactly one option '--use-euclidean-distance' or '--use-cosine-distance' is required "
                                "but none was used");
        } else if(opt_use_euclidean_distance && opt_use_cosine_distance) {
            throw runtime_error("exactly one option '--use-euclidean-distance' or '--use-cosine-distance' is required "
                                "but both were used");
        }
        if(opt_use_norm && opt_use_z_score) {
            throw runtime_error("at most one option '--use-norm' or '--use-z-score' can be used at a time but both "
                                "were used");
        }
        if(opt_use_ti && opt_use_rstar) {
            throw runtime_error("at most one option '--use-ti' or '--use-rstar' can be used at a time but both "
                                "were used");
        }
        if(opt_use_cache && !opt_use_ti) {
            throw runtime_error("option '--use-cache' can be used only with option '--use-ti'");
        }
    } catch(const exception &ex) {
        cerr << ex.what() << endl << endl << desc << endl;
        return 1;
    }

    int dim;
    if(!detectDataDimension(vm["input-file"].as<string>(), dim)) {
        cout << "Error while reading data header." << endl;
        return 1;
    }

    if(dim == 2) {
        DBSCAN<double, 2> dbscan(2, 1, euclidean_distance);
        vector<DataPoint<double, 2>> data;
        if(!dbscan.readData(vm["input-file"].as<string>(), vm["name"].as<string>(), data)) {
            cout << "Error while reading data." << endl;
            return 1;
        }
        dbscan.setOutputFolder(vm["output-dir"].as<string>());
        if(opt_guess_eps) {
            dbscan.setGuessEpsilon();
        } else {
            dbscan.setEpsilon(vm["eps"].as<double>());
            dbscan.setMinNeighbours(vm["min-neighbours"].as<int>());
        }
        if(opt_use_euclidean_distance) {
            dbscan.useDistanceFunction(euclidean_distance);
            dbscan.setSmallerIsCloser();
        } else {
            dbscan.useDistanceFunction(cosine_distance);
            dbscan.setBiggerIsCloser();
        }
        if(opt_use_norm) {
            dbscan.useDataNormalization();
        } else if(opt_use_z_score) {
            dbscan.useZScoreNormalization();
        }
        if(opt_use_plus) {
            dbscan.assignBorderPointsToPossiblyManyClusters();
        }
        if(opt_use_ti) {
            dbscan.useTriangularInequalityIndex();
        } else if(opt_use_rstar) {
            dbscan.useRStarTreeIndex();
        }
        if(opt_use_cache) {
            dbscan.useCache();
        }
        if(opt_bigger_is_closer) {
            dbscan.setBiggerIsCloser();
        } else {
            dbscan.setSmallerIsCloser();
        }

        dbscan.findClusters(data);
    } else if(dim == 16) {
        DBSCAN<double, 16> dbscan(2, 1, euclidean_distance);
        vector<DataPoint<double, 16>> data;
        if(!dbscan.readData(vm["input-file"].as<string>(), vm["name"].as<string>(), data)) {
            cout << "Error while reading data." << endl;
            return 1;
        }
        dbscan.setOutputFolder(vm["output-dir"].as<string>());
        if(opt_guess_eps) {
            dbscan.setGuessEpsilon();
        } else {
            dbscan.setEpsilon(vm["eps"].as<double>());
            dbscan.setMinNeighbours(vm["min-neighbours"].as<int>());
        }
        if(opt_use_euclidean_distance) {
            dbscan.useDistanceFunction(euclidean_distance);
            dbscan.setSmallerIsCloser();
        } else {
            dbscan.useDistanceFunction(cosine_distance);
            dbscan.setBiggerIsCloser();
        }
        if(opt_use_norm) {
            dbscan.useDataNormalization();
        } else if(opt_use_z_score) {
            dbscan.useZScoreNormalization();
        }
        if(opt_use_plus) {
            dbscan.assignBorderPointsToPossiblyManyClusters();
        }
        if(opt_use_ti) {
            dbscan.useTriangularInequalityIndex();
        } else if(opt_use_rstar) {
            dbscan.useRStarTreeIndex();
        }
        if(opt_use_cache) {
            dbscan.useCache();
        }
        if(opt_bigger_is_closer) {
            dbscan.setBiggerIsCloser();
        } else {
            dbscan.setSmallerIsCloser();
        }

        dbscan.findClusters(data);
    } else {
        cout << "This simple version of the demo program accepts only data of dimensions 2 and 16. It is easy and "
                "possible to construct a program that accepts other dimensions but this is only a demo." << endl;
        return 1;
    }


/*
    // Construct an object implementing DBSCAN algorithm with many
    // variants
    // - the first template type (double in this example) can be any
    //   arithmetic type used to represent data points
    // - the second template value is the number of dimensions
    //   of the data points(2 in this example)
    // - the first parameter is an eps - a radius of the area around
    //   core points to search for more cluster point candidates
    //   (4 in our example)
    // - the second parameter is  the minimum number of points that
    //   are in a distance no greater than eps described above to
    //   declare a given point as a core point
    //   (10 in our example).
    // - the last parameter is a pointer to the distance function
    //   function to be used to calculate distance between data points
    //   (euclidean_distance in our example).
    DBSCAN<double, 2> dbscan(4, 10, euclidean_distance);

    // A list of points to be read from the file and then processed
    // by our algorithm. Notice that in our example we expect data
    // points with two dimensions of type double, matching
    // the template parameters of DBSCAN object.
    vector<DataPoint<double, 2>> data;

    // Read data point from file into our list. Reading and
    // processing is deliberately split to allow processing
    // of the same dataset multiple times but with different
    // parameters.
    // Parameters:
    // - path to the file with data ("input.csv" in our example)
    // - name of the dataset used to construct names of output files
    //   to distinguish different processed datasets.
    // - reference to the list of data points - it is filled with
    // data from input file.
    if(!dbscan.readData("input.csv", "complex9", data)) {
        cout << "Error while reading data." << endl;
        return 1;
    }

    // Define many different, optional parameters of the algorithm,
    // for example here we switch from DBSCAN to DBSCAN+
    dbscan.assignBorderPointsToPossiblyManyClusters();

    // Run the algorithm. The results are returned by this method
    // (in this example we ignore them) and are also saved into 2
    // files: with data and with statistics.
    dbscan.findClusters(data);

    // We can change the parameters of the algorithm..
    dbscan.setEpsilon(15);
    dbscan.setMinNeighbours(3);
    dbscan.useTriangularInequalityIndex();

    // .. and run it again to produce 2 next output files with results
    // but reusing the same data.
    dbscan.findClusters(data);
*/


    return 0;
}
