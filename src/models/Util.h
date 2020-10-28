#pragma once

#include <vector>

extern "C" double r2f1_adsj(const double& a, const double& b, const double& c,
                            const double& z);
extern "C" double rf1_adsj(const double& a, const double& b1, const double& b2,
                           const double& c, const double& x, const double& y);

namespace miho {

template <typename T>
int sgn(const T& val);

bool is_approx(const double& x, const double& y);

bool is_approx_int(const double& x, const double& y);

int binomial(int n, int k);

int factorial(int n);

template <typename T>
T pochhammer(const T& val, int n);

double hypergeometric_2F1(const double& a, const double& b, const double& c,
                          const double& z);

double appell_F1(const double& a, const double& b1, const double& b2,
                 const double& c, const double& z1, const double& z2);

double nintegrate_1D(std::function<double(double)> func, const double& x_low,
                     const double& x_upp);

std::vector<double> new_a_list(const std::vector<double>& delta);

}  // namespace miho
