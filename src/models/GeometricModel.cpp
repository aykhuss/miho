#include "GeometricModel.h"

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace miho {

double GeometricModel::pdf(const double& val) const {
  // Eq.(3.14) : assumes j == 1
  return pdf_delta___delta_mu(delta_next(val)) / std::fabs(_sigma.front());
}

double GeometricModel::pdf_delta___delta_mu(const double& delta_next) const {
  // Eq.(4.10) : assumes j == 1
  double pdf_num = pdf_delta__mu(delta_next);
  if (!_q_pdf_den) {
    _pdf_den = pdf_delta__mu(_delta);
    _q_pdf_den = true;
  }
  return pdf_num / _pdf_den;
}

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

double GeometricModel::pdf_delta__mu(const double& delta_next) const {
  std::vector<double> delta = _delta;
  delta.push_back(delta_next);
  return pdf_delta__mu(delta);
}

std::vector<double> GeometricModel::a_list(const std::vector<double>& delta) {
  auto n_delta = delta.size();

  std::vector<double> a;
  a.reserve(n_delta + 1);

  // "guesses"
  a.push_back(std::numeric_limits<double>::infinity());  // a[0]
  for (auto i = 1; i < n_delta; ++i) {
    if (delta[i - 1] != 0.) {
      a.push_back(std::fabs(delta[i] / delta[i - 1]));
    } else {
      a.push_back(std::numeric_limits<double>::infinity());
    }
  }
  a.push_back(0.);  // a[n]

  // patch conflicts
  int i_low = n_delta;
  int i_upp = 0;
  // double val = -1.;
  while (true) {
    bool hit = false;
    for (auto i = 1; i < n_delta; ++i) {
      if (a[i] < a[i + 1]) {
        hit = true;
        if ((i == i_low - 1) || (i == i_upp)) {
          // we're extending a previous range
          i_low = std::min(i_low, i);
          i_upp = std::max(i_upp, i + 1);
        } else {
          // a start of a new range
          i_low = i;
          i_upp = i + 1;
        }
        // merge
        double val = 0.;
        if (delta[i_low - 1] != 0.) {
          val = std::fabs(delta[i_upp] / delta[i_low - 1]);
        } else {
          val = std::numeric_limits<double>::infinity();
        }
        val = std::pow(val, 1. / double(i_upp - i_low + 1));
        // std::cout << "merging range [" << i_low << "," << i_upp << "] = " <<
        // val
        //           << std::endl;
        for (auto j = i_low; j <= i_upp; ++j) {
          a[j] = val;
        }
        break;
      }
    }
    if (!hit) break;
  }

  return a;
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
