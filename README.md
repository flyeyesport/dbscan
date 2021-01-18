# dbscan
DBSCAN implementation with many variants: DBSCAN+, TI-DBSCAN, R*-tree, euclidean and cosine distances, norm and z-score.

You need [boost](https://www.boost.org) library to compile. I tested it with version 1.73 of the library.

## build
```
mkdir build; cd build; cmake ..; make
```

## run, play, test
Use **dbscan** sample program to test the algorithm. It expects an input CSV file with 3 or 17 columns and a header. The algorithm works on data points of any dimension >= 2, but the sample program only accepts dimensions of 2 and 16. For historical reasons the program expects 1 last, extra column in input file which it then ignores.

Run **dbscan** without parameters to see all the options.

It was tested only on linux, but should compile and work on other platforms. 
