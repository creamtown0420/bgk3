#include <iostream>
#include <boost/math/quadrature/tanh_sinh.hpp>

namespace bmq = boost::math::quadrature;

template <typename F>
double dde3d(F func, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, double epsilone)
{
    double result = 0;
    bmq::tanh_sinh<double> integrator(epsilone);

    auto f1 = [&](double x)
    {
        auto f2 = [&](double y)
        {
            auto f3 = [&](double z)
            {
                return func(x, y, z);
            };
            return integrator.integrate(f3, z_min, z_max);
        };
        return integrator.integrate(f2, y_min, y_max);
    };

    result = integrator.integrate(f1, x_min, x_max);
    return result;
}

// int main()
// {
//     auto func = [](double x, double y, double z)
//     {
//         return x * x + y * y + z * z;
//     };
//     double epsilone;
//     double result = dde3d(func, -1, 1, -2, 2, -3, 3, epsilone);
//     std::cout << "Result: " << result << std::endl;

//     return 0;
// }