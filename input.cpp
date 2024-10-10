#include <mpi.h>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <numeric>  // For std::iota
#include <random>  // For random number generation
#include <stdexcept>  // For std::invalid_argument


enum DataType {
    sorted,
    reverseSorted,
    random,
    perturbed
};

bool verify_sorted(const std::vector<int>& array);
std::vector<int> generate_random_array(int size, const DataType dataType);
std::vector<int> generate_data_chunk(int start, int end, const DataType dataType);

std::vector<int> generate_random_array(int size, const DataType dataType) {
    int chunkSize = size / MPI::COMM_WORLD.Get_size();
    int start = rank * chunkSize;
    int end = (rank == MPI::COMM_WORLD.Get_size() - 1) ? size : start + chunkSize;

    std::vector<int> localDataChunk = generate_data_chunk(start, end, dataType);
    std::vector<int> gatheredDataChunk;

    if (rank == 0) {
        gatheredDataChunk.resize(size);
    }

    MPI_Gather(localDataChunk.data(), localDataChunk.size(), MPI_INT,
                    gatheredDataChunk.data(), localDataChunk.size(), MPI_INT,
                    0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::vector<int> finalData;
        for (const auto& chunk : gatheredDataChunk) {
            finalData.insert(finalData.end(), chunk.begin(), chunk.end());
        }
        std::cout << "Generated " << dataType << " data of size " << size << " across " << size << " nodes." << std::endl;
    }
}

std::vector<int> generate_data_chunk(int start, int end, const DataType dataType) {
    std::vector<int> data(end - start);
    if (dataType == sorted) {
        std::iota(data.begin(), data.end(), start);
    } else if (dataType == reverseSorted) {
        std::iota(data.rbegin(), data.rend(), start);
    } else if (dataType == random) {
        std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<int> dist(0, end);
        std::generate(data.begin(), data.end(), [&]() { return dist(rng); });
    } else if (dataType == perturbed) {
        std::iota(data.begin(), data.end(), start);
        std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<int> dist(0, end);
        int perturbation = (end - start) * 0.01;
        for (int i = 0; i < perturbation; ++i) {
            data[dist(rng) % data.size()] = dist(rng);
        }
    } else {
        throw std::invalid_argument("Unknown data type");
    }
    return data;
}

#endif // VERIFY_SORTED_H