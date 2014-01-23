#ifndef LATAN_GLOBAL_HPP_
#define	LATAN_GLOBAL_HPP_

#include <map>
#include <string>
#include <latan/Eigen/Dense>

#define LATAN_BEGIN_CPPDECL namespace Latan {
#define LATAN_END_CPPDECL }

// attribute to switch off unused warnings with gcc
#ifdef __GNUC__
#define __dumb __attribute__((unused))
#else
#define __dumb
#endif

LATAN_BEGIN_CPPDECL

// Environment
namespace Env
{
    extern const std::string fullName;
    extern const std::string name;
    extern const std::string version;
}

// string conversions
template <typename T>
T strTo(std::string str)
{
    T buf;
    std::istringstream stream(str);
    
    stream >> buf;
    
    return buf;
}

template <typename T>
std::string strFrom(T x)
{
    std::ostringstream stream;
    
    stream << x;
    
    return stream.str();
}

// key test in maps
template <typename K, typename T>
bool keyExists(const K &key, const std::map<K, T> &map)
{
    return (map.find(key) != map.end());
}

LATAN_END_CPPDECL

#include <latan/Exceptions.hpp>

#endif