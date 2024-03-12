#include <iostream>
#include<chrono>
#include "Solution.h"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    auto start_time = std::chrono::high_resolution_clock::now();
    SecondSolution usual(10000);
    usual.run(0.001);
    usual.print_result();
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "Time passed: " << duration.count() << " micsec" << std::endl;
    MPI_Finalize();
    return 0;
}