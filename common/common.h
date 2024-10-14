#ifndef COMMON_H
#define COMMON_H

#include <mpi.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <limits.h>

enum DataType {
    sorted,
    reverseSorted,
    rands,
    perturbed
};


// input data generation
std::vector<int> generate_array(int size, const DataType dataType);
std::vector<int> generate_data_chunk(int start, int end, const DataType dataType);

// sort algorithms
std::vector<int> sampleSort(const std::vector<int>& inputArray);
void lsd_radix_sort(std::vector<int>& data, int rank, int size);
int get_digit(int number, int digit_pos);

// sort verification
bool sortVerify(const std::vector<int>& inputArray);

#endif // COMMON_H