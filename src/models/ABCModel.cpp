#include "ABCModel.h"

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace miho {

double ABCModel::pdf(const double& val) const {
  // Eq.(3.14) : assumes j == 1
  return pdf_delta___delta_mu(delta_next(val)) / _sigma.front();
}

double ABCModel::pdf_delta___delta_mu(const double& delta_next) const {
  // Eq.(4.10) : assumes j == 1
  std::vector<double> delta = _delta;
  double pdf_den = pdf_delta__mu(delta);
  // std::cout << "pdf_den = " << pdf_den << std::endl;
  delta.push_back(delta_next);
  double pdf_num = pdf_delta__mu(delta);
  // std::cout << "pdf_num = " << pdf_num << "[" << delta_next << "]" <<
  // std::endl;
  return pdf_num / pdf_den;
}

double ABCModel::pdf_delta__mu(const std::vector<double>& delta) const {
  return a_int(delta);
}

double ABCModel::pdf_delta__mu(const double& delta_next) const {
  std::vector<double> delta = _delta;
  delta.push_back(delta_next);
  return pdf_delta__mu(delta);
}

double ABCModel::a_integrand(const double& a,
                             const std::vector<double>& delta) const {
  double a_val = a;
  if (a_val < 1e3*std::numeric_limits<double>::epsilon()) {
    a_val = 1e3*std::numeric_limits<double>::epsilon();
    return 0.;
  }
  auto m = delta.size() - 1;  // counting starts at 0
  // find \Delta^{(\pm)}
  double dp = 0.;
  double dm = 0.;
  for (auto i = 0; i < delta.size(); ++i) {
    if (delta.at(i) < 0.) {
      dm = std::min(dm, delta.at(i) / pow(a_val, i));
    } else {
      dp = std::max(dp, delta.at(i) / pow(a_val, i));
    }
  }
  // fmt::print("## > (dm,dp) = ({},{})\n", dm, dp);
  // select the correct case
  double b_int = 0.;
  if ((-1. + dp) >= (1 + dm)) {
    b_int = std::pow(2., m + 1 + _epsilon) / (m + _epsilon) /
            std::pow(dp - dm, m + _epsilon);
    // fmt::print("## > case 1: {}\n", b_int);
  } else {
    b_int = 2. / (m + _epsilon) + 2. - (dp - dm);
    // fmt::print("## > case 2: {}\n", b_int);
  }

  return b_int * std::pow(1. - a_val, _omega) /
         std::pow(a_val, m * (m + 1) / 2) / std::pow(2., m + 1) * _epsilon *
         (1. + _omega) / (m + 1 + _epsilon);
}

double ABCModel::a_int(const std::vector<double>& delta) const {
  const size_t max_a_nodes = 500;
  double result = 0.;
  double error = 0.;
  double last_integrand = 0.;
  double last_derivative = 0.;
  double curr_derivative = 0.;
  for (auto i = 0; i < max_a_nodes; ++i) {
    double a = (i + 0.5) / max_a_nodes;
    double integrand = a_integrand(a, delta);
    // fmt::print("# > x, f(x) = {}, {}\n", a, integrand);
    // compute some accumulator
    double dres = (integrand + last_integrand) / 2. / max_a_nodes;
    result += dres;
    // compute the error
    curr_derivative = (integrand - last_integrand) * max_a_nodes;
    error += std::fabs(curr_derivative - last_derivative) *
             pow(1. / max_a_nodes, 2) / 12.;
    last_derivative = curr_derivative;
  }

  // check for accuracy termination condition
  // fmt::print("# > result[{}] = {} +- {}\n", max_a_nodes, result, error);
  // std::cin.ignore();

  return result;
}

}  // namespace miho
