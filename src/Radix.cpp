#include "../common/common.h"

// Function to extract the digit at a specific position
int get_digit(int number, int position) {
    return (number / static_cast<int>(std::pow(10, position))) % 10;
}

// Function to build histogram for current position
std::vector<int> build_histogram(const std::vector<int>& data, int position) {
    std::vector<int> histogram(10, 0);
    for (int number : data) {
        int digit = get_digit(number, position);
        histogram[digit]++;
    }
    return histogram;
}

// Function to perform in-place MSD Radix Sort on a local data segment
void in_place_msd_radix_sort(std::vector<int>& data, int l, int max_digits) {
    if (data.size() <= 1 || l >= max_digits) return;

    // Step 1: Build histogram for current digit position
    std::vector<int> histogram = build_histogram(data, l);

    // Step 2: Calculate heads and tails for buckets
    std::vector<int> heads(10, 0), tails(10, 0);
    tails[0] = histogram[0];
    for (int i = 1; i < 10; i++) {
        heads[i] = heads[i - 1] + histogram[i - 1];
        tails[i] = tails[i - 1] + histogram[i];
    }

    // Step 3: Sort elements in place based on the current digit
    for (int i = 0; i < 10; i++) {
        while (heads[i] < tails[i]) {
            int elem = data[heads[i]];
            int bucket = get_digit(elem, l);
            if (bucket == i) {
                heads[i]++;
            } else {
                std::swap(data[heads[i]], data[heads[bucket]]);
                heads[bucket]++;
            }
        }
    }

    // Step 4: Recursively sort each bucket if there are further digits to process
    for (int i = 0; i < 10; i++) {
        int start = (i == 0) ? 0 : tails[i - 1];
        int end = tails[i];
        if (end - start > 1) {
            std::vector<int> sub_data(data.begin() + start, data.begin() + end);
            in_place_msd_radix_sort(sub_data, l + 1, max_digits);
            std::copy(sub_data.begin(), sub_data.end(), data.begin() + start);
        }
    }
}

int calculateMaxDigits(const std::vector<int>& numbers) {
    if (numbers.empty()) return 0;

    int max_value = *max_element(numbers.begin(), numbers.end());
    return max_value > 0 ? static_cast<int>(log10(max_value)) + 1 : 1;
}

std::vector<int> radixSort(std::vector<int> &arr) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Broadcast the total number of elements to all processes
    int array_size = arr.size();

    CALI_MARK_BEGIN(CALI_COMM);
    CALI_MARK_BEGIN(CALI_COMM_SMALL);

    MPI_Bcast(&array_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    CALI_MARK_END(CALI_COMM_SMALL);
    CALI_MARK_END(CALI_COMM);

    // Scatter data to all processes
    int local_n = array_size / size;
    std::vector<int> local_data(local_n);

    CALI_MARK_BEGIN(CALI_COMM);
    CALI_MARK_BEGIN(CALI_COMM_LARGE);

    MPI_Scatter(arr.data(), local_n, MPI_INT, local_data.data(), local_n, MPI_INT, 0, MPI_COMM_WORLD);

    CALI_MARK_END(CALI_COMM_LARGE);
    CALI_MARK_END(CALI_COMM);

    // Each process performs in-place MSD radix sort on its segment
    CALI_MARK_BEGIN(CALI_COMP);
    CALI_MARK_BEGIN(CALI_COMP_LARGE);
    int max_digits = calculateMaxDigits(arr); // Get the maximum number of digits in the array.
    CALI_MARK_END(CALI_COMP_LARGE);
    CALI_MARK_END(CALI_COMP);

    CALI_MARK_BEGIN(CALI_COMP);
    CALI_MARK_BEGIN(CALI_COMP_LARGE);
    in_place_msd_radix_sort(local_data, 0, max_digits);
    CALI_MARK_END(CALI_COMP_LARGE);
    CALI_MARK_END(CALI_COMP);

    // Vector to gather data from all processes
    CALI_MARK_BEGIN(CALI_COMP);
    CALI_MARK_BEGIN(CALI_COMP_SMALL);

    std::vector<int> sortedData;
    if(rank == 0) {
        sortedData.resize(array_size);
    }

    // Gather sorted segments back to the root process
    MPI_Gather(local_data.data(), local_n, MPI_INT, sortedData.data(), local_n, MPI_INT, 0, MPI_COMM_WORLD);

    CALI_MARK_END(CALI_COMP_SMALL);
    CALI_MARK_END(CALI_COMP);

    CALI_MARK_BEGIN(CALI_COMP);
    CALI_MARK_BEGIN(CALI_COMP_SMALL);
    // Only the root process performs final merging
    if (rank == 0) {
        std::sort(sortedData.begin(), sortedData.end());  // Merge step: final sort
    }
    CALI_MARK_END(CALI_COMP_SMALL);
    CALI_MARK_END(CALI_COMP);

    return sortedData;
}
