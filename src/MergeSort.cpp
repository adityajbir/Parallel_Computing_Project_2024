#include "../common/common.h"

// Merge function for merge sort
void merge(std::vector<int>& arr, int left, int mid, int right) {
    int n1 = mid - left + 1;
    int n2 = right - mid;

    // Temp arrays
    std::vector<int> lhs(n1);
    std::vector<int> rhs(n2);

    // Copy data to temp arrays
    std::copy(arr.begin() + left, arr.begin() + mid + 1, lhs.begin());
    std::copy(arr.begin() + mid + 1, arr.begin() + right + 1, rhs.begin());

    // Merge temp arrays back into arr[left..right]
    int i = 0, j = 0, k = left;

    while (i < n1 && j < n2) {
        if (lhs[i] <= rhs[j]) {
            arr[k++] = lhs[i++];
        } else {
            arr[k++] = rhs[j++];
        }
    }

    // Copy remaining elements
    while (i < n1)
        arr[k++] = lhs[i++];
    while (j < n2)
        arr[k++] = rhs[j++];
}

// Merge sort function
void mergeSort(std::vector<int>& arr, int left, int right) {
    if (left < right) {
        // Start computation small region
        CALI_MARK_BEGIN(CALI_COMP_SMALL);

        int mid = left + (right - left) / 2;
        // Sort first and second halves
        mergeSort(arr, left, mid);
        mergeSort(arr, mid + 1, right);
        // Merge the sorted halves
        merge(arr, left, mid, right);

        // End computation small region
        CALI_MARK_END(CALI_COMP_SMALL);
    }
}

std::vector<int> multiwayMerge(const std::vector<std::vector<int>>& sorted_subarrays) {
    // Use a min-heap to perform multi-way merge
    typedef std::pair<int, std::pair<int, int>> HeapNode; // (value, (array_index, element_index))
    std::priority_queue<HeapNode, std::vector<HeapNode>, std::greater<HeapNode>> minHeap;

    // Initialize heap with the first element of each subarray
    for (size_t i = 0; i < sorted_subarrays.size(); ++i) {
        if (!sorted_subarrays[i].empty()) {
            minHeap.push({sorted_subarrays[i][0], {i, 0}});
        }
    }

    std::vector<int> merged_array;
    while (!minHeap.empty()) {
        HeapNode current = minHeap.top();
        minHeap.pop();

        int value = current.first;
        int array_index = current.second.first;
        int element_index = current.second.second;

        merged_array.push_back(value);

        // If there are more elements in the same subarray, add the next element to the heap
        if (element_index + 1 < sorted_subarrays[array_index].size()) {
            minHeap.push({sorted_subarrays[array_index][element_index + 1], {array_index, element_index + 1}});
        }
    }

    return merged_array;
}

// MPI merge sort function
std::vector<int> mpiMergeSort(std::vector<int>& inputArray) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n; // Total number of elements

    if (rank == 0) {
        n = inputArray.size();
    }

    // Broadcast the total number of elements
    CALI_MARK_BEGIN(CALI_COMM);
    CALI_MARK_BEGIN(CALI_COMM_SMALL);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    CALI_MARK_END(CALI_COMM_SMALL);
    CALI_MARK_END(CALI_COMM);

    // Determine the size of data each process will handle
    int local_n = n / size;
    int remainder = n % size;

    std::vector<int> counts(size);
    std::vector<int> displs(size);

    for (int i = 0; i < size; ++i) {
        counts[i] = local_n + (i < remainder ? 1 : 0);
        displs[i] = (i == 0) ? 0 : displs[i - 1] + counts[i - 1];
    }

    // Allocate local data
    std::vector<int> local_data(counts[rank]);

    // Scatter the data using MPI_Scatterv
    CALI_MARK_BEGIN(CALI_COMM);
    CALI_MARK_BEGIN(CALI_COMM_LARGE);
    MPI_Scatterv(rank == 0 ? inputArray.data() : nullptr, counts.data(), displs.data(), MPI_INT,
                 local_data.data(), counts[rank], MPI_INT, 0, MPI_COMM_WORLD);
    CALI_MARK_END(CALI_COMM_LARGE);
    CALI_MARK_END(CALI_COMM);

    // Local merge sort
    CALI_MARK_BEGIN(CALI_COMP);
    CALI_MARK_BEGIN(CALI_COMP_LARGE);
    mergeSort(local_data, 0, counts[rank] - 1);
    CALI_MARK_END(CALI_COMP_LARGE);
    CALI_MARK_END(CALI_COMP);

    // Gather the sorted subarrays at the root process
    std::vector<int> gathered_data;
    if (rank == 0) {
        gathered_data.resize(n);
    }

    CALI_MARK_BEGIN(CALI_COMM);
    CALI_MARK_BEGIN(CALI_COMM_LARGE);
    MPI_Gatherv(local_data.data(), 
                counts[rank], 
                MPI_INT,
                gathered_data.data(), 
                counts.data(), 
                displs.data(), 
                MPI_INT, 0, 
                MPI_COMM_WORLD
                );
    CALI_MARK_END(CALI_COMM_LARGE);
    CALI_MARK_END(CALI_COMM);

    // Root process performs final merge
    if (rank == 0) {
        // Split gathered data into sorted subarrays based on counts and displs
        std::vector<std::vector<int>> sorted_subarrays(size);
        for (int i = 0; i < size; ++i) {
            sorted_subarrays[i].assign(gathered_data.begin() + displs[i],
                                       gathered_data.begin() + displs[i] + counts[i]);
        }

        // Merge all sorted subarrays
        CALI_MARK_BEGIN(CALI_COMP);
        CALI_MARK_BEGIN(CALI_COMP_LARGE);
        std::vector<int> sortedArray = multiwayMerge(sorted_subarrays);
        CALI_MARK_END(CALI_COMP_LARGE);
        CALI_MARK_END(CALI_COMP);

        return sortedArray;
    } else {
        // Other processes return an empty vector
        return std::vector<int>();
    }
}