#pragma once

namespace miho {

template <typename T>
int sgn(const T& val);

bool is_approx(const double& x, const double& y);

int binomial(int n, int k);

int factorial(int n);

template <typename T>
T pochhammer(const T& val, int n);

double hypergeometric_2F1(const double& a, const double& b, const double& c,
                          const double& z);

double appell_F1(const double& a, const double& b1, const double& b2,
                 const double& c, const double& z1, const double& z2);

}  // namespace miho
