#include "../common/common.h"

std::vector<int> sampleSort(const std::vector<int>& inputArray) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    MPI_Comm_size(MPI_COMM_WORLD, &size);  

    int N;
    if (rank == 0) {
        N = inputArray.size();
    }
    // Step 1: Broadcast the total number of elements to all processes
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Step 2: Compute sendCounts and sendDispls on all processes for MPI_Scatterv
    int base_count = N / size;
    int remainder = N % size;

    std::vector<int> scatterSendCounts(size);
    std::vector<int> scatterSendDispls(size);
    int sum = 0;
    for (int i = 0; i < size; ++i) {
        scatterSendCounts[i] = base_count + (i < remainder ? 1 : 0);
        scatterSendDispls[i] = sum;
        sum += scatterSendCounts[i];
    }

    // Step 3: Distribute data among processes using MPI_Scatterv
    int localSize = scatterSendCounts[rank];
    std::vector<int> localArray(localSize);
    MPI_Scatterv(
        rank == 0 ? inputArray.data() : nullptr, // Send buffer
        scatterSendCounts.data(),                // Send counts
        scatterSendDispls.data(),                // Displacements
        MPI_INT,                                 // Data type
        localArray.data(),                       // Receive buffer
        localSize,                               // Receive count
        MPI_INT,                                 // Data type
        0,                                       // Root process
        MPI_COMM_WORLD                           // Communicator
    );

    // Step 4: Each process sorts its local data
    std::sort(localArray.begin(), localArray.end());

    // Step 5: Each process selects local samples
    int s = size - 1; // Number of samples per process
    std::vector<int> localSamples;
    if (!localArray.empty()) {
        double interval = static_cast<double>(localArray.size()) / (s + 1);
        for (int i = 1; i <= s; ++i) {
            int idx = static_cast<int>(i * interval);
            if (idx >= localArray.size()) {
                idx = localArray.size() - 1;
            }
            localSamples.push_back(localArray[idx]);
        }
    }
    int localSampleSize = localSamples.size();

    // Step 6: Gather sample sizes at the root process
    std::vector<int> recvSampleSizes(size);
    MPI_Gather(
        &localSampleSize,       // Send buffer
        1,                      // Send count
        MPI_INT,                // Data type
        recvSampleSizes.data(), // Receive buffer
        1,                      // Receive count
        MPI_INT,                // Data type
        0,                      // Root process
        MPI_COMM_WORLD          // Communicator
    );

    // Step 7: Compute displacements and total number of samples at the root
    std::vector<int> sampleRecvDispls(size, 0);
    int totalSamples = 0;
    if (rank == 0) {
        for (int i = 0; i < size; ++i) {
            sampleRecvDispls[i] = totalSamples;
            totalSamples += recvSampleSizes[i];
        }
    }

    // Step 8: Gather all samples at the root process
    std::vector<int> gatheredSamples;
    if (rank == 0) {
        gatheredSamples.resize(totalSamples);
    }
    MPI_Gatherv(
        localSamples.data(),     // Send buffer
        localSampleSize,         // Send count
        MPI_INT,                 // Data type
        gatheredSamples.data(),  // Receive buffer
        recvSampleSizes.data(),  // Receive counts
        sampleRecvDispls.data(), // Displacements
        MPI_INT,                 // Data type
        0,                       // Root process
        MPI_COMM_WORLD           // Communicator
    );

    // Step 9: Root process selects splitters from gathered samples
    std::vector<int> splitters(size - 1);
    if (rank == 0) {
        std::sort(gatheredSamples.begin(), gatheredSamples.end());
        double interval = static_cast<double>(gatheredSamples.size()) / size;
        for (int i = 1; i < size; ++i) {
            int idx = static_cast<int>(i * interval);
            if (idx >= gatheredSamples.size()) {
                idx = gatheredSamples.size() - 1;
            }
            splitters[i - 1] = gatheredSamples[idx];
        }
    }

    // Step 10: Broadcast splitters to all processes
    MPI_Bcast(
        splitters.data(), // Data buffer
        size - 1,         // Number of elements
        MPI_INT,          // Data type
        0,                // Root process
        MPI_COMM_WORLD    // Communicator
    );

    // Step 11: Each process partitions its local data based on splitters
    std::vector<std::vector<int>> buckets(size);
    for (const auto& val : localArray) {
        int idx = std::upper_bound(splitters.begin(), splitters.end(), val) - splitters.begin();
        buckets[idx].push_back(val);
    }

    // Step 12: Prepare sendCountsAlltoall and sdisplsAlltoall for MPI_Alltoallv
    std::vector<int> sendCountsAlltoall(size);
    std::vector<int> sdisplsAlltoall(size);
    int sendTotal = 0;
    for (int i = 0; i < size; ++i) {
        sendCountsAlltoall[i] = buckets[i].size();
    }
    sdisplsAlltoall[0] = 0;
    for (int i = 1; i < size; ++i) {
        sdisplsAlltoall[i] = sdisplsAlltoall[i - 1] + sendCountsAlltoall[i - 1];
    }
    sendTotal = sdisplsAlltoall[size - 1] + sendCountsAlltoall[size - 1];

    // Step 13: Concatenate all buckets into sendBuf
    std::vector<int> sendBuf(sendTotal);
    for (int i = 0; i < size; ++i) {
        std::copy(buckets[i].begin(), buckets[i].end(), sendBuf.begin() + sdisplsAlltoall[i]);
    }

    // Step 14: Exchange bucket sizes using MPI_Alltoall
    std::vector<int> recvCountsAlltoall(size);
    MPI_Alltoall(
        sendCountsAlltoall.data(), // Send buffer
        1,                         // Send count
        MPI_INT,                   // Data type
        recvCountsAlltoall.data(), // Receive buffer
        1,                         // Receive count
        MPI_INT,                   // Data type
        MPI_COMM_WORLD             // Communicator
    );

    // Step 15: Prepare recv displacements for MPI_Alltoallv
    std::vector<int> rdisplsAlltoall(size);
    int recvTotal = 0;
    rdisplsAlltoall[0] = 0;
    for (int i = 0; i < size; ++i) {
        recvTotal += recvCountsAlltoall[i];
        if (i > 0) {
            rdisplsAlltoall[i] = rdisplsAlltoall[i - 1] + recvCountsAlltoall[i - 1];
        }
    }

    // Step 16: Allocate receive buffer
    std::vector<int> recvBuf(recvTotal);

    // Step 17: Exchange data using MPI_Alltoallv
    MPI_Alltoallv(
        sendBuf.data(),               // Send buffer
        sendCountsAlltoall.data(),    // Send counts
        sdisplsAlltoall.data(),       // Send displacements
        MPI_INT,                      // Data type
        recvBuf.data(),               // Receive buffer
        recvCountsAlltoall.data(),    // Receive counts
        rdisplsAlltoall.data(),       // Receive displacements
        MPI_INT,                      // Data type
        MPI_COMM_WORLD                // Communicator
    );

    // Step 18: Each process sorts its received data
    std::sort(recvBuf.begin(), recvBuf.end());

    // Step 19: Gather sorted data sizes at the root process
    int localSortedSize = recvBuf.size();
    std::vector<int> sortedRecvSizes(size);
    MPI_Gather(
        &localSortedSize,   // Send buffer
        1,                  // Send count
        MPI_INT,            // Data type
        sortedRecvSizes.data(), // Receive buffer
        1,                  // Receive count
        MPI_INT,            // Data type
        0,                  // Root process
        MPI_COMM_WORLD      // Communicator
    );

    // Step 20: Compute displacements for final gathering at the root
    std::vector<int> sortedRecvDispls(size);
    if (rank == 0) {
        sortedRecvDispls[0] = 0;
        for (int i = 1; i < size; ++i) {
            sortedRecvDispls[i] = sortedRecvDispls[i - 1] + sortedRecvSizes[i - 1];
        }
    }

    // Step 21: Gather sorted data back to the root process
    std::vector<int> sortedData;
    if (rank == 0) {
        sortedData.resize(N);
    }
    MPI_Gatherv(
        recvBuf.data(),        // Send buffer
        localSortedSize,       // Send count
        MPI_INT,               // Data type
        sortedData.data(),     // Receive buffer
        sortedRecvSizes.data(),// Receive counts
        sortedRecvDispls.data(),// Displacements
        MPI_INT,               // Data type
        0,                     // Root process
        MPI_COMM_WORLD         // Communicator
    );

    // Step 22: Return the sorted data on the root process
    if (rank == 0) {
        return sortedData;
    } else {
        return std::vector<int>(); // Empty vector for non-root processes
    }
}
