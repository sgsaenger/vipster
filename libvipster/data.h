#ifndef DATA_H
#define DATA_H

#include <vector>
#include <numeric>

#include "vec.h"

namespace Vipster{

struct BaseData{
    std::string name{};
    virtual ~BaseData() = default;
};

template<size_t N, typename T>
struct DataGrid: public BaseData, private std::vector<T>
{
    using Extent = std::array<size_t, N>;
    using std::vector<T>::begin;
    using std::vector<T>::end;
    using std::vector<T>::at;
    using std::vector<T>::operator[];
    using iterator = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;
    Mat cell;
    Vec origin;
    const Extent extent;
    const size_t size;
    static constexpr size_t Dim = N;
    DataGrid(std::array<size_t, N> extent)
        : extent{extent},
          size{std::accumulate(extent.begin(), extent.end(),
                               1ul, std::multiplies<size_t>())}
    {
        std::vector<T>::resize(size);
    }
    // TODO: c++17?
//    template<typename ...Args,
//             typename = typename std::enable_if<sizeof...(Args) == N>::type>
//    DataGrid(Args... args)
//        : extent{args...}, size{(... * args)}
//    {
//        std::vector<T>::resize(size);
//    }
};

template<typename T>
struct DataGrid2D: public DataGrid<2, T>
{
    T& operator()(size_t x, size_t y)
    {
        return (*this)[this->extent[0]*y+x];
    }
    T operator()(size_t x, size_t y)const
    {
        return (*this)[this->extent[0]*y+x];
    }
    using DataGrid<2, T>::DataGrid;
};

template<typename T>
struct DataGrid3D: public DataGrid<3, T>
{
    T& operator()(size_t x, size_t y, size_t z)
    {
        return (*this)[this->extent[1]*this->extent[0]*z +
                       this->extent[0]*y +
                       x];
    }
    T operator()(size_t x, size_t y, size_t z)const
    {
        return (*this)[this->extent[1]*this->extent[0]*z +
                       this->extent[0]*y +
                       x];
    }
    using DataGrid<3, T>::DataGrid;
};

using DataGrid2D_f = DataGrid2D<float>;
using DataGrid3D_f = DataGrid3D<float>;
using DataGrid3D_v = DataGrid3D<Vec>;

}

#endif // DATA_H
