#include <iostream>
#include <functional>
#include <boost/math/quadrature/gauss_kronrod.hpp>
namespace bmq = boost::math::quadrature;

template <typename F>
double quad3(F func, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, double epsilon)
{
    double result = 0;
    bmq::gauss_kronrod<double, 15> integrator;

    auto f1 = [&](double x)
    {
        auto f2 = [&](double y)
        {
            auto f3 = [&](double z)
            {
                return func(x, y, z);
            };

            return integrator.integrate(f3, z_min, z_max, epsilon);
        };

        return integrator.integrate(f2, y_min, y_max, epsilon);
    };

    result = integrator.integrate(f1, x_min, x_max, epsilon);

    return result;
}

template <typename F>
double quad2(F func, double x_min, double x_max, double y_min, double y_max, double epsilon)
{
    double result = 0;
    bmq::gauss_kronrod<double, 15> integrator;

    auto f1 = [&](double x)
    {
        auto f2 = [&](double y)
        {
            return func(x, y);
        };

        return integrator.integrate(f2, y_min, y_max, epsilon);
    };

    result = integrator.integrate(f1, x_min, x_max, epsilon);

    return result;
}