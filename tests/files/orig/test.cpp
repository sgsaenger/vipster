#include <vector>
#include <list>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <string>
#include <iostream>
#include <iomanip>
#include <limits>
#include <array>
#include <stddef.h>
#include <bitset>
#include <sstream>
#include <map>
//#include <mpi.h>

struct test_A{
    int t;
    double a,b,c;
};

struct test_B{
    int t;
    double a,b,c;
    double q;
};

int main()
{
    std::map<size_t, double> m;
    using mm = std::map<size_t, double>;
    std::cout 
        << sizeof(mm::key_type) << std::endl
        << sizeof(mm::mapped_type) << std::endl
        << sizeof(mm::value_type) << std::endl
        << sizeof(mm::node_type) << std::endl
        ;
}
