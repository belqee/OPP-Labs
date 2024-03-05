#include <iostream>

#include "Solution.h"

int main(int argc, char **argv)
{
    AllSolution usual(10000);
    usual.run(0.00001);
    usual.print_result();
    return 0;
}