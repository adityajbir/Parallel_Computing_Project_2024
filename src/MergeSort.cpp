#include "../common/common.h"

// Merge function (basically a splitter)
void merge(std::vector<int>& arr, int left, int mid, int right) {
    int n1 = mid - left + 1;
    int n2 = right - mid;

    // Temp arrays needed for copying data over 
    std::vector<int> lhs(n1);
    std::vector<int> rhs(n2);

    std::copy(arr.begin() + left, arr.begin() + mid + 1, lhs.begin());
    std::copy(arr.begin() + mid + 1, arr.begin() + right + 1, rhs.begin());

    
    int i = 0, j = 0, k = left;
    //merge temps into original
    while (i < n1 && j < n2) {
        if (lhs[i] <= rhs[j]) {
            arr[k++] = lhs[i++];
        } else {
            arr[k++] = rhs[j++];
        }
    }

    // Copy remaining elements if left over 

    while (i < n1)
        arr[k++] = lhs[i++];
    while (j < n2)
        arr[k++] = rhs[j++];
}

// Merge sort function(gatherer)
void mergeSort(std::vector<int>& arr, int left, int right) {
    if (left < right) {
        
        int mid = left + (right - left) / 2;
        // Sort first and second halves
        mergeSort(arr, left, mid);
        mergeSort(arr, mid + 1, right);
        // Merge the sorted halves
        merge(arr, left, mid, right);

        
    }
}

std::vector<int> multiwayMerge(const std::vector<std::vector<int>>& sorted_subarrays) {
    // Use a mininmum-heap to perform multi-way merge

    //define heap 
    typedef std::pair<int, std::pair<int, int>> HeapNode; 
    std::priority_queue<HeapNode, std::vector<HeapNode>, std::greater<HeapNode>> minHeap;

    // Initialize heap with the first element of each subarray
    for (size_t i = 0; i < sorted_subarrays.size(); ++i) {
        if (!sorted_subarrays[i].empty()) {
            minHeap.push({sorted_subarrays[i][0], {i, 0}});
        }
    }

    std::vector<int> merged_array;

    // add elements to final array as long as it is not emoty
    while (!minHeap.empty()) {
        HeapNode current = minHeap.top();
        minHeap.pop();

        int value = current.first;
        int array_index = current.second.first;
        int element_index = current.second.second;

        merged_array.push_back(value);

        // If there are more elements from same subarray, add the  element to heap
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


    // Total elements
    int n; 

    //fetches size of the input array we need from argument
    if (rank == 0) {
        n = inputArray.size();
    }

    //comm call    
    CALI_MARK_BEGIN(CALI_COMM);
    CALI_MARK_BEGIN(CALI_COMM_SMALL);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    CALI_MARK_END(CALI_COMM_SMALL);
    CALI_MARK_END(CALI_COMM);

    // Determine the size of data each process will handle(should be close to equal for each)
    int local_n = n / size;
    int remainder = n % size;

    //count and displacement vecs for distrubition
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

    // Local merge sort call
    CALI_MARK_BEGIN(CALI_COMP);
    CALI_MARK_BEGIN(CALI_COMP_LARGE);
    mergeSort(local_data, 0, counts[rank] - 1);
    CALI_MARK_END(CALI_COMP_LARGE);
    CALI_MARK_END(CALI_COMP);

    // Gather the sorted subarrays from workers at the root process
    std::vector<int> gathered_data;
    if (rank == 0) {
        gathered_data.resize(n);
    }

    //comm call
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

        // default return call 
        return std::vector<int>();
    }
}