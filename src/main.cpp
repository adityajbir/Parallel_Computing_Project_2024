#include "../common/common.h"

int main(int argc, char** argv) {
    // Start Caliper main region
    CALI_MARK_BEGIN(CALI_MAIN);

    // Initialize MPI
    CALI_MARK_BEGIN("MPI_Init");
    MPI_Init(&argc, &argv);
    CALI_MARK_END("MPI_Init");

    int rank;
    int size;

    CALI_MARK_BEGIN("MPI_Comm_rank");
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    CALI_MARK_END("MPI_Comm_rank");

    CALI_MARK_BEGIN("MPI_Comm_size");
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    CALI_MARK_END("MPI_Comm_size");

    // Initialize Adiak
    adiak::init(MPI_COMM_WORLD);

    // Collect metadata
    adiak::launchdate();    // Launch date of the job
    adiak::libraries();     // Libraries used
    adiak::cmdline(argc, argv); // Command line used to launch the job
    adiak::clustername();   // Name of the cluster

    if (argc < 2) {
        if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " <arraySize>" << std::endl;
        }
        CALI_MARK_BEGIN("MPI_Finalize");
        MPI_Finalize(); // Finalize MPI before exiting
        CALI_MARK_END("MPI_Finalize");
        return 1;
    }

    int arraySize = std::stoi(argv[1]);
    if (argc > 2 && std::stoi(argv[2]) != -1) {
        arraySize = std::pow(2, std::stoi(argv[2]));
    }

    // Collect remaining metadata
    std::string algorithm = "Sample Sort"; // Replace with your actual algorithm
    adiak::value("algorithm", algorithm);
    adiak::value("programming_model", "MPI");
    std::string data_type = "int";
    adiak::value("data_type", data_type);
    int size_of_data_type = sizeof(int);
    adiak::value("size_of_data_type", size_of_data_type);
    adiak::value("input_size", arraySize);
    std::string input_type = "Random"; // Change based on input data
    adiak::value("input_type", input_type);
    adiak::value("num_procs", size);
    std::string scalability = "strong"; // or "weak" depending on your experiment
    adiak::value("scalability", scalability);
    adiak::value("group_num", 2);
    adiak::value("implementation_source", "handwritten");

    // Data initialization
    CALI_MARK_BEGIN(CALI_DATA_INIT_RUNTIME);
    std::vector<int> arr = generate_array(arraySize, rands);
    CALI_MARK_END(CALI_DATA_INIT_RUNTIME);

    if (rank == 0) {
        std::cout << "Input Array: " << std::endl;
        for (int i = 0; i < arraySize; i++) {
            std::cout << arr[i] << " ";
        }
        std::cout << std::endl << "#######################################" << std::endl;
    }

    std::vector<int> sortedArr = sampleSort(arr); // Testing sampleSort function

    if (rank == 0) {
        std::cout << "Sorted Array: " << std::endl;
        for (int i = 0; i < arraySize; i++) {
            std::cout << sortedArr[i] << " ";
        }
        std::cout << std::endl << "#######################################" << std::endl;
    }

    // Correctness check
    CALI_MARK_BEGIN(CALI_CORRECTNESS_CHECK);
    bool isSorted = sortVerify(sortedArr);
    CALI_MARK_END(CALI_CORRECTNESS_CHECK);

    // Only the master process (rank 0) prints the result
    if (rank == 0) {
        if (isSorted) {
            std::cout << "The array is sorted." << std::endl;
        } else {
            std::cout << "The array is not sorted." << std::endl;
        }
    }

    // End Caliper main region
    CALI_MARK_END(CALI_MAIN);

    // Finalize MPI
    CALI_MARK_BEGIN("MPI_Finalize");
    MPI_Finalize();
    CALI_MARK_END("MPI_Finalize");

    return 0;
}
