#ifndef COMMON_H
#define COMMON_H

#include <mpi.h>
#include <vector>
#include <iostream>

enum DataType {
    sorted,
    reverseSorted,
    rands,
    perturbed
};

bool sortVerify(const std::vector<int>& array);
std::vector<int> generate_random_array(int size, const DataType dataType);
std::vector<int> generate_data_chunk(int start, int end, const DataType dataType);

#endif // COMMON_H