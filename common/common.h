#ifndef COMMON_H
#define COMMON_H

#include <mpi.h>
#include <caliper/cali.h>
#include <caliper/cali-manager.h>
#include <adiak.hpp>
#include <vector>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <limits.h>
#include <memory>
#include <cmath> 
#include <unordered_map>
#include <queue>
#include <functional>
#include <utility>

enum DataType {
    sorted,
    reverseSorted,
    rands,
    perturbed
};
/* Define Caliper region names using macros */
#define CALI_DATA_INIT_RUNTIME "data_init_runtime"
#define CALI_CORRECTNESS_CHECK "correctness_check"
#define CALI_COMM "comm"
#define CALI_COMM_SMALL "comm_small"
#define CALI_COMM_LARGE "comm_large"
#define CALI_COMP "comp"
#define CALI_COMP_SMALL "comp_small"
#define CALI_COMP_LARGE "comp_large"


// input data generation
std::vector<int> generate_array(int size, const DataType dataType);
std::vector<int> generate_data_chunk(int start, int end, const DataType dataType);

// sort algorithms
std::vector<int> sampleSort(std::vector<int>& inputArray);
std::vector<int> radixSort(std::vector<int>& inputArray);
std::vector<int> mpiMergeSort(std::vector<int>& inputArray);
std::vector<int> bitonicSort(std::vector<int>& inputArray);

// sort verification
bool sortVerify(const std::vector<int>& inputArray);

#endif // COMMON_H
