#include "../common/common.h"

void compare_low(int j, std::vector<int>& local, int local_size, int rank, int size) {
    int min;

    int partner = rank ^ (1 << j);

    // Send local maximum to partner
    MPI_Send(&local[local_size - 1], 1, MPI_INT, partner, 0, MPI_COMM_WORLD);

    // Receive min from partner
    MPI_Recv(&min, 1, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Init send and rec buffers
    std::vector<int> buffer_send(local_size + 1);
    std::vector<int> buffer_receive(local_size + 1);

    int send_counter = 0;
    for (int i = local_size - 1; i >= 0; i--) {
        if (local[i] > min) {
            send_counter++;
            buffer_send[send_counter] = local[i];
        } else {
            break;
        }
    }

    buffer_send[0] = send_counter;  // First element is the count
    MPI_Send(buffer_send.data(), send_counter + 1, MPI_INT, partner, 0, MPI_COMM_WORLD);

    MPI_Recv(buffer_receive.data(), local_size + 1, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    std::vector<int> temp_array = local;  // Copy local array

    int buffer_size = buffer_receive[0];
    int k = 1;  // Index in buffer_receive
    int m = 0;  // Index in temp_array

    for (int i = 0; i < local_size; i++) {
        if (m < temp_array.size() && (k > buffer_size || temp_array[m] <= buffer_receive[k])) {
            local[i] = temp_array[m];
            m++;
        } else if (k <= buffer_size) {
            local[i] = buffer_receive[k];
            k++;
        }
    }
    std::sort(local.begin(), local.end());
}

void compare_high(int j, std::vector<int>& local, int local_size, int rank, int size) {
    int max;

    int partner = rank ^ (1 << j);

    // Receive max from partner
    MPI_Recv(&max, 1, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Send local minimum to partner
    MPI_Send(&local[0], 1, MPI_INT, partner, 0, MPI_COMM_WORLD);

    std::vector<int> buffer_send(local_size + 1);
    std::vector<int> buffer_receive(local_size + 1);

    int send_counter = 0;
    for (int i = 0; i < local_size; i++) {
        if (local[i] < max) {
            send_counter++;
            buffer_send[send_counter] = local[i];
        } else {
            break;
        }
    }

    MPI_Recv(buffer_receive.data(), local_size + 1, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    buffer_send[0] = send_counter;
    MPI_Send(buffer_send.data(), send_counter + 1, MPI_INT, partner, 0, MPI_COMM_WORLD);

    std::vector<int> temp_array = local;

    int buffer_size = buffer_receive[0];
    int k = 1;
    int m = local_size - 1;

    for (int i = local_size - 1; i >= 0; i--) {
        if (m >= 0 && (k > buffer_size || temp_array[m] >= buffer_receive[k])) {
            local[i] = temp_array[m];
            m--;
        } else if (k <= buffer_size) {
            local[i] = buffer_receive[k];
            k++;
        }
    }

    std::sort(local.begin(), local.end());
}

std::vector<int> bitonicSort(std::vector<int>& inputArray) {
    int rank, P;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    // Ensure the number of processes is a power of 2
    if ((P & (P - 1)) != 0) {
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

    // Compute local size and prepare for scattering
    int base_count = N / P;
    int remainder = N % P;

    std::vector<int> scatterSendCounts(P);
    std::vector<int> scatterSendDispls(P);
    int sum = 0;
    for (int i = 0; i < P; ++i) {
        scatterSendCounts[i] = base_count + (i < remainder ? 1 : 0);
        scatterSendDispls[i] = sum;
        sum += scatterSendCounts[i];
    }

    int local_size = scatterSendCounts[rank];
    std::vector<int> local(local_size);

    // Distribute data among processes
    CALI_MARK_BEGIN(CALI_COMM);
    CALI_MARK_BEGIN(CALI_COMM_LARGE);
    MPI_Scatterv(
        rank == 0 ? inputArray.data() : nullptr,
        scatterSendCounts.data(),
        scatterSendDispls.data(),
        MPI_INT,
        local.data(),
        local_size,
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

    // Bitonic merge network
    int dimension = static_cast<int>(std::log2(P));
    CALI_MARK_BEGIN(CALI_COMP);
    CALI_MARK_BEGIN(CALI_COMP_LARGE);
    for (int i = 0; i < dimension; i++) {
        for (int j = i; j >= 0; j--) {
            if (((rank >> (i + 1)) % 2 == 0 && (rank >> j) % 2 == 0) ||
                ((rank >> (i + 1)) % 2 != 0 && (rank >> j) % 2 != 0)) {
                compare_low(j, local, local_size, rank, P);
            } else {
                compare_high(j, local, local_size, rank, P);
            }
        }
    }
    CALI_MARK_END(CALI_COMP_LARGE);
    CALI_MARK_END(CALI_COMP);

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
        local_size,
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