#include "Util.h"

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace miho {

template <typename T>
int sgn(const T& val) {
  return (T(0) < val) - (val < T(0));
}
// Explicit instantiations.
template int sgn<int>(const int&);
template int sgn<double>(const double&);

bool is_approx(const double& x, const double& y) {
  // 1.0 is here because epsilon is the smallest difference
  // that can be distinguished from 1.0
  double max_val = std::max({1.0, std::fabs(x), std::fabs(y)});
  return std::fabs(x - y) <= std::numeric_limits<double>::epsilon() * max_val;
}

int binomial(int n, int k) {
  if (k == 0 || k == n) return 1;
  return binomial(n - 1, k) + binomial(n - 1, k - 1);
}

int factorial(int n) {
  if (n == 0 || n == 1) return 1;
  return n * factorial(n - 1);
}

template <typename T>
T pochhammer(const T& val, int n) {
  T res = T(1);
  for (auto i = 0; i < n; ++i) {
    res *= val + T(i);
  }
  return res;
}
// Explicit instantiations.
template int pochhammer<int>(const int&, int);
template double pochhammer<double>(const double&, int);

// https://journal.r-project.org/archive/2015/RJ-2015-022/RJ-2015-022.pdf
double hypergeometric_2F1(const double& a, const double& b, const double& c,
                          const double& z) {
  // fmt::print("2F1({},{},{};{})\n", a, b, c, z);
  const double target_accuracy = std::numeric_limits<double>::epsilon();
  // const double target_accuracy = 1e-11;
  const size_t max_it = 1000;
  if (std::fabs(z) >= 1.) {
    throw "hypergeometric_2F1: won't converge";
  }
  // if (std::fabs(a) > std::fabs(b)) {
  //   // canonical ordering: |a| always smaller than |b|
  //   return hypergeometric_2F1(b, a, c, z);
  // }
  if (z > 0.5) {
    // (0.5,1) --> Eq.(15.3.6)
    return std::tgamma(c) * std::tgamma(c - a - b) / std::tgamma(c - a) /
               std::tgamma(c - b) *
               hypergeometric_2F1(a, b, a + b - c + 1., 1. - z) +
           std::pow(1. - z, c - a - b) * std::tgamma(c) *
               std::tgamma(a + b - c) / std::tgamma(a) / std::tgamma(b) *
               hypergeometric_2F1(c - a, c - b, c - a - b + 1., 1. - z);
  }
  if (z < 0.0) {
    // (-1,0) --> Eq.(15.3.4-5)
    if (std::fabs(c - a) > std::fabs(c - b)) {
      //  Eq.(15.3.4)
      return std::pow(1. - z, -a) *
             hypergeometric_2F1(a, c - b, c, z / (z - 1.));
    } else {
      //  Eq.(15.3.5)
      return std::pow(1. - z, -b) *
             hypergeometric_2F1(b, c - a, c, z / (z - 1.));
    }
  }
  // [0,0.5]
  if (std::max(a, b) > std::max(c - a, c - b)) {
    // fmt::print("swap!\n");
    return std::pow(1. - z, c - a - b) * hypergeometric_2F1(c - a, c - b, c, z);
  }
  double acc = 0.;
  for (auto it = 0; it < max_it; ++it) {
    double dres = 1.;
    for (auto i = 0; i < it; ++i) {
      dres *= (a + double(i)) * (b + double(i)) / (c + double(i)) * z /
              double(i + 1);
    }
    //> this one is unstable:
    // double dres = pochhammer(a, it) * pochhammer(b, it) / pochhammer(c, it) /
    //               factorial(it) * std::pow(z, it);
    // fmt::print("hypergeometric_2F1: {} <-> {}  |  {} {}\n", it, dres, acc, std::fabs(dres / acc));
    if (!std::isfinite(dres)) {
      throw "hypergeometric_2F1: failed to converge!";
      return acc;
    }
    acc += dres;
    if (std::fabs(dres / acc) <= target_accuracy) {
      // fmt::print("reached target {}\n", std::fabs(dres / acc));
      break;
    }
  }
  return acc;
}

// http://www.gasaneofisica.uns.edu.ar/papers/2001/ColavecchiaGasaneoMiragliacpc_01_138_29.pdf
double appell_F1(const double& a, const double& b1, const double& b2,
                 const double& c, const double& z1, const double& z2) {
  // const double target_accuracy = std::numeric_limits<double>::epsilon();
  const double target_accuracy = 1e-6;
  const size_t max_it = 20;
  if (std::fabs(z1) >= 1. || std::fabs(z2) >= 1.) {
    throw "appell_F1: won't converge";
  }
  // if (std::fabs(z1) > std::fabs(z2)) {
  //   // canonical ordering: |z1| always smaller than |z2|
  //   return appell_F1(a, b2, b1, c, z2, z1);
  // }
  double acc = 0.;
  for (auto it = 0; it < max_it; ++it) {
    double dres = 1.;
    for (auto i = 0; i < it; ++i) {
      dres *= (a + double(i)) * (b1 + double(i)) * (b2 + double(i)) *
              (c - a + double(i)) / (c + double(it - 1 + i)) /
              (c + double(2 * i)) / (c + double(2 * i + 1)) * (z1 * z2) /
              double(i + 1);
    }
    dres *= hypergeometric_2F1(a+it, b1+it, c+2*it, z1);
    dres *= hypergeometric_2F1(a+it, b2+it, c+2*it, z2);
    if (!std::isfinite(dres)) {
      throw "appell_F1: failed to converge!";
      return acc;
    }
    acc += dres;
    if (std::fabs(dres / acc) <= target_accuracy) {
      // fmt::print("reached target {}\n", std::fabs(dres / acc));
      break;
    }
  }
  return acc;
}

//> old and unstable
// double appell_F1(const double& a, const double& b1, const double& b2,
//                  const double& c, const double& z1, const double& z2) {
//   // const double target_accuracy = std::numeric_limits<double>::epsilon();
//   const double target_accuracy = 1e-6;
//   const size_t max_it = 20;
//   const double p_damp = 0.3;
//   if (std::fabs(z1) >= 1. || std::fabs(z2) >= 1.) {
//     throw "appell_F1: won't converge";
//   }
//   if (std::fabs(z1) > std::fabs(z2)) {
//     // ensure that |z1| always smaller than |z2|
//     return appell_F1(a, b2, b1, c, z2, z1);
//   }
//   int ratio = int(std::ceil(std::pow(std::fabs(z2 / z1), p_damp)));
//   fmt::print("---> {} {} => ratio {}\n", z1, z2, ratio);
//   // ratio = 1;
//   double acc = 0.;
//   for (auto it = 0; it < max_it; ++it) {
//     double max_dres = 0.;
//     for (auto i1 = 0; i1 <= it; ++i1) {
//       for (auto i2 = (it - i1) * ratio; i2 < (it - i1 + 1) * ratio; ++i2) {
//         double dres = pochhammer(b1, i1) * pochhammer(b2, i2) *
//                       pochhammer(a, i1 + i2) / pochhammer(c, i1 + i2) /
//                       factorial(i1) / factorial(i2) * std::pow(z1, i1) *
//                       std::pow(z2, i2);
//         if (!std::isfinite(dres)) {
//           throw "appell_F1: failed to converge!";
//           return acc;
//         }
//         fmt::print("appell_F1: {} ({}, {}) <-> {}  |  {}\n", it, i1, i2, dres,
//                    acc);
//         max_dres = std::max(max_dres, std::fabs(dres));
//         acc += dres;
//       }
//     }
//     fmt::print("---> max_dres = {}\n", max_dres);
//     if (max_dres <= target_accuracy) break;
//   }
//   return acc;
// }

}  // namespace miho
