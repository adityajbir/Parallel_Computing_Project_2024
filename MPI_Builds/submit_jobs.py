import subprocess
import time

input_sizes = [16, 18, 20, 22, 24, 26, 28]
input_types = ["sorted", "random", "reverseSorted", "perturbed"]
num_procs_list = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
alg_types = ["sampleSort", "radixSort", "mpiMergeSort", "bitonicSort"]
max_tasks_per_node = 32 

for input_size in input_sizes:
    for input_type in input_types:
        for num_procs in num_procs_list:
            for alg_type in alg_types:
                if num_procs <= max_tasks_per_node:
                    nodes = 1
                    ntasks_per_node = num_procs
                else:
                    nodes = num_procs // max_tasks_per_node
                    ntasks_per_node = max_tasks_per_node

                command = [
                    "sbatch",
                    f"--nodes={nodes}",
                    f"--ntasks-per-node={ntasks_per_node}",
                    "mpi.grace_job",
                    str(input_size),
                    input_type,
                    alg_type
                ]

                # For debugging purposes, you can print the command
                print("Submitting job:", ' '.join(command))

                # Submit the job
                subprocess.run(command)
                time.sleep(5)
