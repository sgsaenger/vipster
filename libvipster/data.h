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
// shape information
    using Extent = std::array<size_t, N>;
    using value_type = T;
    static constexpr size_t Dim = N;
// constructors
    DataGrid(Extent extent)
        : extent{extent},
          size{std::accumulate(extent.begin(), extent.end(),
                               1ul, std::multiplies<size_t>())}
    {
        std::vector<T>::resize(size);
    }
    template<typename ...Args,
             typename = std::enable_if_t<sizeof...(Args) == N>,
             typename = std::enable_if_t<std::conjunction_v<std::is_integral<Args>...>>>
    DataGrid(Args... args)
        : extent{args...}, size{(... * args)}
    {
        std::vector<T>::resize(size);
    }
// expose data with real dimension
private:
    template<typename ...Args,
             typename = std::enable_if_t<sizeof...(Args) == N>,
             typename = std::enable_if_t<std::conjunction_v<std::is_integral<Args>...>>>
    constexpr size_t index(Args... args) const
    {
        size_t idx{}, i{1};
        ((idx += args * std::accumulate(extent.begin()+i++, extent.end(), 1, std::multiplies<size_t>())), ...);
        return idx;
    }
    constexpr size_t index(const Extent& arg) const
    {
        size_t idx{arg[Dim-1]};
        for(size_t i=0; i<Dim-1; ++i){
            idx += arg[i] * std::accumulate(extent.begin()+i+1, extent.end(), 1, std::multiplies<size_t>());
        }
        return idx;
    }
public:
    template <typename ...Args>
    T& operator()(Args... args)
    {
        return (*this)[index(args...)];
    }
    T& operator()(const Extent& arg)
    {
        return (*this)[index(arg)];
    }
    template <typename ...Args>
    T operator()(Args... args) const
    {
        return (*this)[index(args...)];
    }
    T operator()(const Extent& arg) const
    {
        return (*this)[index(arg)];
    }
// expose flat data via vector-interface
public:
    using std::vector<T>::begin;
    using std::vector<T>::end;
    using std::vector<T>::at;
    using std::vector<T>::operator[];
    using iterator = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;
// data-members
public:
    Mat cell{};
    Vec origin{};
    const Extent extent;
    const size_t size;
};

using DataGrid2D_f = DataGrid<2, double>;
using DataGrid3D_f = DataGrid<3, double>;
using DataGrid3D_v = DataGrid<3, Vec>;

}

#endif // DATA_H
