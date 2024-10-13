#include <mpi.h>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>


enum DataType {
    sorted,
    reverseSorted,
    rands,
    perturbed
};

std::vector<int> generate_data_chunk(int start, int end, const DataType dataType);

std::vector<int> generate_random_array(int arraySize, const DataType dataType) {
    int rank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    int chunkSize = arraySize / numProcs;
    int remainder = arraySize % numProcs;
    int start = rank * chunkSize + std::min(rank, remainder);
    int end = start + chunkSize + (rank < remainder ? 1 : 0);

    std::vector<int> localDataChunk = generate_data_chunk(start, end, dataType);
    std::vector<int> gatheredDataChunk;

    if (rank == 0) {
        gatheredDataChunk.resize(arraySize);
    }

    std::vector<int> counts(numProcs), displacements(numProcs);
    for (int i = 0; i < numProcs; ++i) {
        counts[i] = (i < remainder) ? chunkSize + 1 : chunkSize;
        displacements[i] = (i > 0) ? displacements[i - 1] + counts[i - 1] : 0;
    }

    MPI_Gatherv(localDataChunk.data(), localDataChunk.size(), MPI_INT,
                gatheredDataChunk.data(), counts.data(), displacements.data(), MPI_INT,
                0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::vector<int> finalData;
	    std::cout << "Generated " << dataType << " data of size " << arraySize << " across " << numProcs << " nodes." << std::endl;
        return gatheredDataChunk;
    }
    return {};
}

std::vector<int> generate_data_chunk(int start, int end, const DataType dataType) {
    std::vector<int> data(end - start);
    if (dataType == sorted) {
        std::iota(data.begin(), data.end(), start);
    } else if (dataType == reverseSorted) {
        std::iota(data.rbegin(), data.rend(), start);
    } else if (dataType == rands) {
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

// // Main function for testing
// int main(int argc, char *argv[]) {
//     MPI_Init(&argc, &argv);

//     int rank, size;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &size);

//     // Example usage
//     int array_size = 1000; // Example array size
//     DataType data_type = random; // Example data type

//     std::vector<int> data = generate_random_array(array_size, data_type);

//     if (rank == 0) {
//         std::cout << "Generated array: ";
//         for (int i = 0; i < data.size(); i++) {
//             std::cout << data[i] << " ";
//         }
//         std::cout << std::endl;
//     }

//     MPI_Finalize();
//     return 0;
// }