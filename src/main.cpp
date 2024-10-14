#include "../common/common.h"

int main(int argc, char** argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);

    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2) {
        if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " <array_size>" << std::endl;
        }
        MPI_Finalize(); // Finalize MPI before exiting
        return 1;
    }

    int array_size = std::stoi(argv[1]);

    std::vector<int> arr = generate_array(array_size, sorted);
    if(rank == 0) {
        for(int i = 0; i < array_size; i++){
              std::cout << arr[i] << " ";
        }
    }
    std::cout << std::endl;

    // Call the sortVerify function
    bool isSorted = sortVerify(arr);

    // Only the master process (rank 0) prints the result
    if (rank == 0) {
       if (isSorted) {
            std::cout << "The array is sorted." << std::endl;
        } else {
            std::cout << "The array is not sorted." << std::endl;
        }
    }

    // Finalize MPI
    MPI_Finalize();

    return 0;
}