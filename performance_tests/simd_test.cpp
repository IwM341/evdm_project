#include <array>
#include <cassert>
#include <cstddef>
#include <experimental/simd>
#include <functional>
#include <iostream>
#include <numeric>
namespace stdx = std::experimental;

float rnd(){
    return rand()/(RAND_MAX+0.0f);
}

template <typename T>
void println(std::string_view name, T const& a)
{
    std::cout << name << ": ";
    for (std::size_t i{}; i != std::size(a); ++i)
        std::cout << a[i] << ' ';
    std::cout << '\n';
}

int main()
{
    using V = stdx::native_simd<double>;
    using W = stdx::fixed_size_simd<float, 8>;
    
    W X([](size_t i){return rnd();});
    W Y([](size_t i){return rnd();});
    W Z = X + 0.5f*Y;
    
    println("X",X);
    println("Y",Y);
    println("Z",Z);
    
    return 0;//system("g++ main.cpp -O3 -S -fverbose-asm -march=native -mavx2 -o out.asm");
}