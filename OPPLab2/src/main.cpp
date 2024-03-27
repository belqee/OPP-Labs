#include <iostream>
#include<chrono>
#include "Solution.h"

int main(int argc, char **argv)
{
    omp_set_num_threads(2);

        auto start_time = std::chrono::high_resolution_clock::now();
        FirstSolution usual(9600);
        usual.run(0.0001);
        //usual.print_result();

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        std::cout << "Time passed: " << duration.count() << " micsec" << std::endl;
    return 0;
}