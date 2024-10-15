#include "../common/common.h"

void merge(const std::vector<int>& a, const std::vector<int>& b, std::vector<int>& result, bool ascending) {
    size_t i = 0, j = 0, k = 0;
    size_t n = a.size(), m = b.size();
    result.resize(n + m);
    while (i < n && j < m) {
        if ((a[i] <= b[j]) == ascending)
            result[k++] = a[i++];
        else
            result[k++] = b[j++];
    }
    while (i < n)
        result[k++] = a[i++];
    while (j < m)
        result[k++] = b[j++];
    result.resize(k); // Adjust size if there are dummy elements
}

void keep_data(const std::vector<int>& merged_data, std::vector<int>& local_data, bool keep_low) {
    size_t n = local_data.size();
    if (keep_low) {
        std::copy(merged_data.begin(), merged_data.begin() + n, local_data.begin());
    } else {
        std::copy(merged_data.end() - n, merged_data.end(), local_data.begin());
    }
}

std::vector<int> BitonicSort()(std::vector<int> arr) {
    int rank, P; // Rank and num processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    int N_padded;
    std::vector<int> padded_arr;
    if (rank == 0) {
        N_padded = pow(2, ceil(log2(arr.size())));

        std::vector<int> padded_arr = arr;
        padded_arr.resize(N_padded, INT_MAX);    
    }

    // Broadcast the total number of elements to all processes
    MPI_Bcast(&N_padded, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int base_count = N_padded / P;
    int remainder = N_padded % P;

    // Calculate the displacements and counts for MPI_Scatterv
    std::vector<int> scatterCounts(P);
    std::vector<int> scatterDisplacements(P);
    int sum = 0;
    for (int i = 0; i < P; i++) {
        scatterCounts[i] = base_count + (i < remainder ? 1 : 0);
        scatterDisplacements[i] = sum;
        sum += scatterCounts[i];
    }

    // Scatter data to all processes
    int localSize = scatterCounts[rank];
    std::vector<int> local_data(localSize);
    MPI_Scatterv(
        rank == 0 ? padded_arr.data() : nullptr,
        scatterCounts.data(),
        scatterDisplacements.data(),
        MPI_INT,
        local_data.data(),
        localSize,
        MPI_INT,
        0,
        MPI_COMM_WORLD
    );
    
    // Sort local data
    sort(local_data.begin(), local_data.end());

    // Start Bitonic mergesort
    int max_level = static_cast<int>(ceil(log2(P)));
    for (int k = 1; k <= max_level; ++k) {
        for (int j = k; j > 0; --j) {
            int distance = 1 << (j - 1);
            int partner = rank ^ distance;

            if (parner >= P) { continue; } // Skip non-existen partners

            // Determine the sort direction
            bool ascending = ((my_rank >> (k - 1)) & 1) == 0;

            // Exchange data with partner
            std::vector<int> partner_data(localSize);   

            // Non-blocking send and receive to prevent deadlocks
            MPI_Request requests[2];
            MPI_Isend(local_data.data(), localSize, MPI_INT, partner, 0, MPI_COMM_WORLD, &requests[0]);
            MPI_Irecv(partner_data.data(), localSize, MPI_INT, partner, 0, MPI_COMM_WORLD, &requests[1]);
            MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);

            // Merge local data and partner data
            std::vector<int> merged_data;
            merge(local_data, partner_data, merged_data, ascending);

            // Keep either the lower or upper half based on the sort direction
            bool keep_low = (my rank < parnter) == ascending;
            keep_data(merged_data, local_data, keep_low);

        }
    }

    // Gather the sorted data backs to process 0
    std::vector<int> sortedArr;
    if (my_rank == 0) {
        // Prepare receive counts and displacements for MPI_Gatherv
        std::vector<int> recv_counts(num_procs);
        std::vector<int> displacements(num_procs);
        int sum = 0;
        for (int i = 0; i < num_procs; ++i) {
            recv_counts[i] = base + (i < remainder ? 1 : 0);
            displacements[i] = sum;
            sum += recv_counts[i];
        }

        sortedArr.resize(sum);
        MPI_Gatherv(
            local_data.data(),
            scatterCounts[rank],
            MPI_INIT,
            sortedArr.data(),
            recv_counts.data(),
            displacements.data(),
            MPI_INIT,
            0,
            MPI_COMM_WORLD
        );
        sortedArr.resize(arr.size());
    } 
    else {
        // Non-root processes send their data
        MPI_Gatherv(
            local_data.data(),
            scatterCounts[rank],
            MPI_INIT,
            nullptr,
            nullptr,
            nullptr,
            MPI_INIT,
            0,
            MPI_COMM_WORLD
        );
    }
    // Return the sorted array on rank 0; empty array on other ranks
    return sortedArr;    
}

