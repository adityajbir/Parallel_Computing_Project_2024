#include "../common/common.h"

int main(int argc, char** argv) {
    // Start Caliper main region
    CALI_CXX_MARK_FUNCTION

    // Initialize MPI
    MPI_Init(&argc, &argv);

    int rank;
    int size;

    CALI_MARK_BEGIN("MPI_Comm_rank");
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    CALI_MARK_END("MPI_Comm_rank");

    CALI_MARK_BEGIN("MPI_Comm_size");
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    CALI_MARK_END("MPI_Comm_size");

    // Initialize Adiak
    adiak::init(NULL);

    // Collect metadata
    adiak::launchdate();    // Launch date of the job
    adiak::libraries();     // Libraries used
    adiak::cmdline();       // Command line used to launch the job
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

    int arraySize = std::pow(2, std::stoi(argv[1]));
    std::string inputTypeStr = argv[2];
    std::string algorithmStr = argv[3];
    // Map inputTypeStr to DataType enum
    std::unordered_map<std::string, DataType> inputTypeMap = {
        {"sorted", sorted},
        {"reverseSorted", reverseSorted},
        {"random", rands},
        {"perturbed", perturbed}
    };

    if (inputTypeMap.find(inputTypeStr) == inputTypeMap.end()) {
        if (rank == 0) {
            std::cerr << "Invalid input type: " << inputTypeStr << std::endl;
        }
        CALI_MARK_BEGIN("MPI_Finalize");
        MPI_Finalize(); // Finalize MPI before exiting
        CALI_MARK_END("MPI_Finalize");
        return 1;
    }

    DataType inputType = inputTypeMap[inputTypeStr];

    // Map algorithmStr to function pointers
    std::unordered_map<std::string, std::vector<int>(*)(std::vector<int>&)> algorithmMap = {
        {"sampleSort", sampleSort},
        {"radixSort", radixSort},
        {"mpiMergeSort", mpiMergeSort}
    };

    if (algorithmMap.find(algorithmStr) == algorithmMap.end()) {
        if (rank == 0) {
            std::cerr << "Invalid algorithm: " << algorithmStr << std::endl;
        }
        CALI_MARK_BEGIN("MPI_Finalize");
        MPI_Finalize(); // Finalize MPI before exiting
        CALI_MARK_END("MPI_Finalize");
        return 1;
    }

    auto algorithmFunc = algorithmMap[algorithmStr];

    // Collect remaining metadata
    adiak::value("algorithm", algorithmStr);
    adiak::value("programming_model", "MPI");
    adiak::value("data_type", "int");
    int size_of_data_type = sizeof(int);
    adiak::value("size_of_data_type", size_of_data_type);
    adiak::value("input_size", arraySize);
    adiak::value("input_type", inputTypeStr);
    adiak::value("num_procs", size);
    std::string scalability = "strong"; // or "weak" depending on your experiment
    adiak::value("scalability", scalability);
    adiak::value("group_num", 2);
    adiak::value("implementation_source", "handwritten");

    // Data initialization
    CALI_MARK_BEGIN(CALI_DATA_INIT_RUNTIME);
    std::vector<int> arr = generate_array(arraySize, inputType);
    CALI_MARK_END(CALI_DATA_INIT_RUNTIME);

    int printSize = std::min(40, arraySize);
    if (rank == 0) {
        std::cout << "Input Array: " << std::endl;
        for (int i = 0; i < printSize; i++) {
            std::cout << arr[i] << " ";
        }
        std::cout << std::endl << "#######################################" << std::endl;
    }

    std::vector<int>sortedArr= algorithmFunc(arr);

    if (rank == 0) {
        std::cout << "Sorted Array: " << std::endl;
        for (int i = 0; i < printSize; i++) {
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

    // Finalize MPI
    MPI_Finalize();

    return 0;
}
