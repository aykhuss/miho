#include "GeometricModel.h"

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace miho {

double GeometricModel::pdf_delta__mu(const std::vector<double>& delta) const {
  // Eq.(4.16) : assumes omega is int
  auto m = delta.size() - 1;  // counting starts at 0
  std::vector<double> a = a_list(delta);
  // for (auto i = 0; i < a.size(); ++i) {
  //   std::cout << "| " << i << ": " << a[i] << std::endl;
  // }
  double acc_k = 0.;
  for (auto k = 0; k <= m; ++k) {
    double acc_j = 0.;
    for (auto j = 0; j <= _omega; ++j) {
      double a_low = std::min(1., a[k + 1]);
      double a_upp = std::min(1., a[k]);
      if (is_approx(a_low, a_upp)) continue;
      if (a_low > a_upp) {
        std::cout << "invalid a list for k=" << k << ": "
                  << "a_low = " << a_low << "[" << a[k + 1]
                  << "], a_upp = " << a_upp << "[" << a[k] << "]" << std::endl;
      }
      acc_j += std::pow(-1, j) * binomial(_omega, j) *
               a_integral(a_low, a_upp, _epsilon, m, k, j);
      // std::cout << k << ":" << j << ": " << acc_j << std::endl;
    }
    double abs_del_k = std::fabs(delta[k]);
    if (abs_del_k < std::numeric_limits<double>::epsilon())
      abs_del_k = std::numeric_limits<double>::epsilon();
    acc_k += acc_j / pow(abs_del_k, m + _epsilon);
  }
  return acc_k * _epsilon * (1. + _omega) / std::pow(2., m) / (m + _epsilon);
}

double GeometricModel::a_integral(const double& a_lower, const double& a_upper,
                                  const double& epsilon, int m, int k, int j) {
  double exponent = (m + epsilon) * k - m * (m + 1) / 2. + j + 1;
  if (exponent == 0.) {
    return std::log(a_upper / a_lower);
  } else {
    return (std::pow(a_upper, exponent) - std::pow(a_lower, exponent)) /
           exponent;
  }
}

}  // namespace miho
