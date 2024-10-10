#include "common.h"

bool sortVerify(const std::vector<int>& inputArray) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int N = inputArray.size();

    int localSize = N / size;
    int remainder = N % size;
    if (rank < remainder) {
        localSize += 1;
    }

    std::vector<int> local_array(localSize);
    std::vector<int> sendCounts(size);
    std::vector<int> sendDispls(size);

    if (rank == 0) {
        int sum = 0;
        for (int i = 0; i < size; ++i) {
            sendCounts[i] = N / size;
            if (i < remainder) {
                sendCounts[i] += 1;
            }
            sendDispls[i] = sum;
            sum += sendCounts[i];
        }
    }

    // Scatter the inputArray to all processes
    MPI_Scatterv(
        inputArray.data(),       // Send buffer (root process)
        sendCounts.data(),       // Send counts
        sendDispls.data(),       // Displacements
        MPI_INT,                 // Data type
        local_array.data(),      // Receive buffer
        localSize,               // Receive count
        MPI_INT,                 // Data type
        0,                       // Root process
        MPI_COMM_WORLD           // Communicator
    );

    // Local sorting verification
    for (int i = 0; i < localSize - 1; ++i) {
        if (local_array[i] > local_array[i + 1]) {
            return false; // local array is not sorted
        }
    }
    //if sorted then front is min and back is max
    int local_min = local_array.front(); 
    int local_max = local_array.back();

    // Boundary verification
    bool boundary_sorted = true;
    int prev_max = 0;

    if (rank < size - 1) { // send the last element from the local array to the next process
        MPI_Send(&local_max, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    }

    if (rank > 0) { // receive the last element from the previous process
        MPI_Recv(&prev_max, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (prev_max > local_min) {
            boundary_sorted = false;
        }
    }

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