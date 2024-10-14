#include "common.h"


bool sortVerify(const std::vector<int>& inputArray) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int N = 0;
    if (rank == 0) {
        N = inputArray.size();
    }
    // Broadcast N to all processes
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int localSize = N / size;
    int remainder = N % size;
    if (rank < remainder) {
        localSize += 1;
    }

    std::vector<int> local_array(localSize > 0 ? localSize : 1);
    std::vector<int> sendCounts(size, 0);
    std::vector<int> sendDispls(size, 0);
    if (rank == 0) {
        int sum = 0;
        for (int i = 0; i < size; ++i) {
            sendCounts[i] = N / size + (i < remainder ? 1 : 0);
            sendDispls[i] = sum;
            sum += sendCounts[i];
        }
    }

    const int* sendbuf = nullptr;
    int* sendCountsPtr = nullptr;
    int* sendDisplsPtr = nullptr;
    if (rank == 0) {
        sendbuf = inputArray.data();
        sendCountsPtr = sendCounts.data();
        sendDisplsPtr = sendDispls.data();
    }

    int dummy;
    int* recvbuf = localSize > 0 ? local_array.data() : &dummy;

    // Scatter the inputArray to all processes
    MPI_Scatterv(
        sendbuf,              // Send buffer (root process)
        sendCountsPtr,        // Send counts
        sendDisplsPtr,        // Displacements
        MPI_INT,              // Data type
        recvbuf,              // Receive buffer
        localSize,            // Receive count
        MPI_INT,              // Data type
        0,                    // Root process
        MPI_COMM_WORLD        // Communicator
    );

    // Local sorting verification
    for (int i = 0; i < localSize - 1; ++i) {
        if (local_array[i] > local_array[i + 1]) {
            return false;
        }
    }

    // Boundary verification
    int local_min = (localSize > 0) ? local_array.front() : INT_MAX;
    int local_max = (localSize > 0) ? local_array.back() : INT_MIN;

    int prev_max = INT_MIN;
    if (rank < size - 1) {
        MPI_Send(&local_max, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    }
    if (rank > 0) {
        MPI_Recv(&prev_max, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    bool boundary_sorted = true;
    if (rank > 0 && localSize > 0) {
        if (prev_max > local_min) {
            boundary_sorted = false;
        }
    }

    // Global reduction to determine if the array is sorted
    bool global_result = false;
    MPI_Allreduce(
        &boundary_sorted, // sendbuf
        &global_result,   // recvbuf
        1,                // count
        MPI_C_BOOL,       // datatype
        MPI_LAND,         // op: reduction operation (logical AND)
        MPI_COMM_WORLD    // comm
    );
    return global_result;
}