# CSCE 435 Group project

## 0. Group number: 2  

## 1. Group members:
1. First: Aditya Biradar (Merge Sort)
2. Second: Eyad Nazir 
3. Third: Eduardo Alvarez
4. Fourth: Juan Carrasco 

## 1a. Method of Communication:
- We will use normal phone messaging as our method of coomunication for this project

## 2. Project topic (e.g., parallel sorting algorithms)

- For this project, we will be implementing various sorting algorithms for parallel computing and do a comparitive analysis of the algorithms to identify pros and cons of each algorithm.


### 2a. Brief project description (what algorithms will you be comparing and on what architectures)
- Architecture: For all the algorithms below we will be implementing them with an MPI architecture
- Bitonic Sort:
- Sample Sort:
- Merge Sort:
- Radix Sort:

### 2b. Pseudocode for each parallel algorithm
- For MPI programs, include MPI calls you will use to coordinate between processes

- Bitonic Sort:
    code 

- Sample Sort:
    code 


- Merge Sort(Aditya - I will finish this later):
    MergeSort():
        Check length of array len(array):
            If length == 0 || length == 1:
                return base case of the array 
            else:
                continue
        Split array into two halves:
            split=len(array)/2 **round down for whole number**
        

        




- Radix Sort:
    code 


### 2c. Evaluation plan - what and how will you measure and compare
- Input sizes, Input types
- Strong scaling (same problem size, increase number of processors/nodes)
- Weak scaling (increase problem size, increase number of processors)
