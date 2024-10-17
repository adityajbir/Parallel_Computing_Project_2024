#include "../common/common.h"

void compare_low(int j, std::vector<int>& local, int elements_per_proc, int rank, int size) {
    CALI_MARK_BEGIN(CALI_COMM);
    int min;

    int partner = rank ^ (1 << j);

    // Send local maximum to partner
    MPI_Send(&local[elements_per_proc - 1], 1, MPI_INT, partner, 0, MPI_COMM_WORLD);

    // Receive min from partner
    MPI_Recv(&min, 1, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    CALI_MARK_END(CALI_COMM_SMALL);

    // Prepare buffers
    std::vector<int> buffer_send(elements_per_proc + 1);
    std::vector<int> buffer_receive(elements_per_proc + 1);

    int send_counter = 0;
    for (int i = elements_per_proc - 1; i >= 0; i--) {
        if (local[i] > min) {
            send_counter++;
            buffer_send[send_counter] = local[i];
        } else {
            break;
        }
    }

    CALI_MARK_BEGIN(CALI_COMM_LARGE);

    buffer_send[0] = send_counter;  // First element is the count
    MPI_Send(buffer_send.data(), send_counter + 1, MPI_INT, partner, 0, MPI_COMM_WORLD);

    MPI_Recv(buffer_receive.data(), elements_per_proc + 1, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    CALI_MARK_END(CALI_COMM_LARGE);
    CALI_MARK_END(CALI_COMM);

    CALI_MARK_BEGIN(CALI_COMP);
    CALI_MARK_BEGIN(CALI_COMP_SMALL);

    std::vector<int> temp_array = local;  // Copy local array

    int buffer_size = buffer_receive[0];
    int k = 1;  // Index in buffer_receive
    int m = 0;  // Index in temp_array

    for (int i = 0; i < elements_per_proc; i++) {
        if (m < temp_array.size() && (k > buffer_size || temp_array[m] <= buffer_receive[k])) {
            local[i] = temp_array[m];
            m++;
        } else if (k <= buffer_size) {
            local[i] = buffer_receive[k];
            k++;
        }
    }

    CALI_MARK_END(CALI_COMP_SMALL);

    CALI_MARK_BEGIN(CALI_COMP_LARGE);
    std::sort(local.begin(), local.end());
    CALI_MARK_END(CALI_COMP_LARGE);
    CALI_MARK_END(CALI_COMP);
}

void compare_high(int j, std::vector<int>& local, int elements_per_proc, int rank, int size) {
    CALI_MARK_BEGIN(CALI_COMM);
    int max;

    int partner = rank ^ (1 << j);

    // Receive max from partner
    MPI_Recv(&max, 1, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Send local minimum to partner
    MPI_Send(&local[0], 1, MPI_INT, partner, 0, MPI_COMM_WORLD);

    CALI_MARK_END(CALI_COMM_SMALL);

    std::vector<int> buffer_send(elements_per_proc + 1);
    std::vector<int> buffer_receive(elements_per_proc + 1);

    int send_counter = 0;
    for (int i = 0; i < elements_per_proc; i++) {
        if (local[i] < max) {
            send_counter++;
            buffer_send[send_counter] = local[i];
        } else {
            break;
        }
    }

    CALI_MARK_BEGIN(CALI_COMM_LARGE);

    MPI_Recv(buffer_receive.data(), elements_per_proc + 1, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    buffer_send[0] = send_counter;
    MPI_Send(buffer_send.data(), send_counter + 1, MPI_INT, partner, 0, MPI_COMM_WORLD);

    CALI_MARK_END(CALI_COMM_LARGE);
    CALI_MARK_END(CALI_COMM);

    CALI_MARK_BEGIN(CALI_COMP);
    CALI_MARK_BEGIN(CALI_COMP_SMALL);

    std::vector<int> temp_array = local;

    int buffer_size = buffer_receive[0];
    int k = 1;
    int m = elements_per_proc - 1;

    for (int i = elements_per_proc - 1; i >= 0; i--) {
        if (m >= 0 && (k > buffer_size || temp_array[m] >= buffer_receive[k])) {
            local[i] = temp_array[m];
            m--;
        } else if (k <= buffer_size) {
            local[i] = buffer_receive[k];
            k++;
        }
    }

    CALI_MARK_END(CALI_COMP_SMALL);

    CALI_MARK_BEGIN(CALI_COMP_LARGE);
    std::sort(local.begin(), local.end());
    CALI_MARK_END(CALI_COMP_LARGE);
    CALI_MARK_END(CALI_COMP);
}

std::vector<int> bitonicSort(std::vector<int>& inputArray) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Ensure the number of processes is a power of 2
    if ((size & (size - 1)) != 0) {
        if (rank == 0) {
            std::cerr << "Number of processes must be a power of 2 for Bitonic sort." << std::endl;
        }
        MPI_Finalize();
        exit(1);
    }

    int N = 0;
    if (rank == 0) {
        N = inputArray.size();
    }

    // Broadcast N to all processes
    CALI_MARK_BEGIN(CALI_COMM);
    CALI_MARK_BEGIN(CALI_COMM_SMALL);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    CALI_MARK_END(CALI_COMM_SMALL);
    CALI_MARK_END(CALI_COMM);

    // Compute elements_per_proc and prepare for scattering
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

    int elements_per_proc = scatterSendCounts[rank];
    std::vector<int> local(elements_per_proc);

    // Distribute data among processes
    CALI_MARK_BEGIN(CALI_COMM);
    CALI_MARK_BEGIN(CALI_COMM_LARGE);
    MPI_Scatterv(
        rank == 0 ? inputArray.data() : nullptr,
        scatterSendCounts.data(),
        scatterSendDispls.data(),
        MPI_INT,
        local.data(),
        elements_per_proc,
        MPI_INT,
        0,
        MPI_COMM_WORLD
    );
    CALI_MARK_END(CALI_COMM_LARGE);
    CALI_MARK_END(CALI_COMM);

    // Local sort
    CALI_MARK_BEGIN(CALI_COMP);
    CALI_MARK_BEGIN(CALI_COMP_LARGE);
    std::sort(local.begin(), local.end());
    CALI_MARK_END(CALI_COMP_LARGE);
    CALI_MARK_END(CALI_COMP);

    // Bitonic merge
    int dimension = static_cast<int>(std::log2(size));

    for (int i = 0; i < dimension; i++) {
        for (int j = i; j >= 0; j--) {
            if (((rank >> (i + 1)) % 2 == 0 && (rank >> j) % 2 == 0) ||
                ((rank >> (i + 1)) % 2 != 0 && (rank >> j) % 2 != 0)) {
                compare_low(j, local, elements_per_proc, rank, size);
            } else {
                compare_high(j, local, elements_per_proc, rank, size);
            }
        }
    }

    // Gather sorted data
    std::vector<int> gatherRecvCounts = scatterSendCounts;
    std::vector<int> gatherRecvDispls = scatterSendDispls;

    std::vector<int> sortedData;
    if (rank == 0) {
        sortedData.resize(N);
    }

    CALI_MARK_BEGIN(CALI_COMM);
    CALI_MARK_BEGIN(CALI_COMM_LARGE);

    MPI_Gatherv(
        local.data(),
        elements_per_proc,
        MPI_INT,
        sortedData.data(),
        gatherRecvCounts.data(),
        gatherRecvDispls.data(),
        MPI_INT,
        0,
        MPI_COMM_WORLD
    );

    CALI_MARK_END(CALI_COMM_LARGE);
    CALI_MARK_END(CALI_COMM);

    // Return sorted data on rank 0
    if (rank == 0) {
        return sortedData;
    } else {
        return std::vector<int>();
    }
}