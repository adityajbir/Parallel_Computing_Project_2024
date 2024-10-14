#include "../common/common.h"

#define NUM_DIGITS 10 // Base-10 for radix sort

// Function prototypes
int get_digit(int number, int digit_pos);

// Main LSD Radix Sort function
void lsd_radix_sort(std::vector<int>& data, int rank, int size) {
    int n = data.size();
    int max_num = *std::max_element(data.begin(), data.end());
    int max_digit_pos = static_cast<int>(std::log10(max_num)) + 1;

    for (int digit_pos = 0; digit_pos < max_digit_pos; digit_pos++) {
        std::vector<int> histogram(NUM_DIGITS, 0);
        std::vector<int> local_histogram(NUM_DIGITS, 0);

        // Step 1: Build local histogram for each process
        int start = rank * (n / size);
        int end = (rank == size - 1) ? n : start + (n / size);

        // Debugging output
        std::cout << "Rank " << rank << " processing range [" << start << ", " << end << ")" << std::endl;

        for (int i = start; i < end; i++) {
            int digit = get_digit(data[i], digit_pos);
            local_histogram[digit]++;
        }

        // Step 2: Combine local histograms into a global histogram
        MPI_Allreduce(local_histogram.data(), histogram.data(), NUM_DIGITS, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        // Step 3: Compute prefix sum of the global histogram
        std::vector<int> prefix_sum(NUM_DIGITS, 0);
        for (int i = 1; i < NUM_DIGITS; i++) {
            prefix_sum[i] = prefix_sum[i - 1] + histogram[i - 1];
        }

        // Step 4: Distribute data based on prefix sum
        std::unique_ptr<int[]> sorted_data(new int[n]);
        std::unique_ptr<int[]> local_sorted_data(new int[end - start]);
        std::vector<int> local_count(NUM_DIGITS, 0);

        for (int i = start; i < end; i++) {
            int digit = get_digit(data[i], digit_pos);
            int pos = prefix_sum[digit] + local_count[digit]++;
            local_sorted_data[pos - start] = data[i];
        }

        MPI_Allgather(local_sorted_data.get(), end - start, MPI_INT, sorted_data.get(), end - start, MPI_INT, MPI_COMM_WORLD);

        // Copy sorted data back to original array
        for (int i = 0; i < n; i++) {
            data[i] = sorted_data[i];
        }
    }
}

int get_digit(int number, int digit_pos) {
    for (int i = 0; i < digit_pos; i++) {
        number /= 10;
    }
    return number % 10;
}