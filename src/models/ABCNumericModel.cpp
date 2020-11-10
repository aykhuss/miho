#include "ABCNumericModel.h"

#include <cubature.h>
#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace {

struct ABCdata {
  const std::vector<double>& delta;
  const double& epsilon;
  const double& xi;
  int omega;
};

int ab_integrand(unsigned ndim, const double* x, void* fdata, unsigned fdim,
                 double* fval) {
  if (ndim != 2) return 1;
  if (fdim != 1) return 1;
  fval[0] = 0.;
  double result = 1.;

  //> a mapping
  double a = x[0];
  // // for numerical stability
  // if (miho::is_approx(a, 0.)) a = std::numeric_limits<double>::epsilon();

  //> b-mapping
  //>  (a) x/(1-x^2)
  // // const double omx1sq = 1. - x[1] * x[1];
  // const double omx1sq = (1. - x[1]) * (1. + x[1]);  // more stable
  // const double b = x[1] / omx1sq;
  // result *= (1. + x[1] * x[1]) / (omx1sq * omx1sq);  // Jac for b-mapping
  //>  (b) arctanh(x)
  const double b = std::atanh(x[1]);
  result /= (1. - x[1]) * (1. + x[1]);  // Jac for b-mapping

  const ABCdata* abc = (ABCdata*)fdata;
  double m = abc->delta.size() - 1.;  // counting starts at 0
  //> prefactor
  result *= abc->epsilon * (1. + abc->omega) /
            (std::pow(2., m + 1.) * abc->xi * (m + 1. + abc->epsilon));
  //> a-dependent piece
  result *= std::pow(1. - a, abc->omega) / std::pow(a, m * (m + 1.) / 2.);
  //> b-dependent piece (max function)
  double max_val = std::max(1., std::fabs(b) / abc->xi);
  double pow_a = 1.;  // accumulate a-powers
  for (auto i = 1; i <= m; ++i) {
    pow_a *= a;  // pow_a = a^i
    double dtest = std::fabs(abc->delta.at(i) / pow_a - b);
    if (dtest > max_val) max_val = dtest;
  }
  result /= std::pow(max_val, m + 1. + abc->epsilon);
  //> done.
  if (!std::isfinite(result)) {
    std::cerr << "#ab_integrand: problem for a = " << a << ", b = " << b
              << std::endl;
    return 0;
  }
  fval[0] = result;
  return 0;  // success
}

}  // namespace

namespace miho {

double ABCNumericModel::pdf(const double& val) const {
  // Eq.(3.14) : assumes j == 1
  return pdf_delta___delta_mu(delta_next(val)) / std::fabs(_sigma.front());
}

double ABCNumericModel::pdf_delta___delta_mu(const double& delta_next) const {
  // Eq.(4.10) : assumes j == 1
  double pdf_num = pdf_delta__mu(delta_next);
  if (!_q_pdf_den) {
    _pdf_den = pdf_delta__mu(_delta);
    _q_pdf_den = true;
  }
  return pdf_num / _pdf_den;
}

double ABCNumericModel::pdf_delta__mu(const std::vector<double>& delta) const {
  // cubature library
  const size_t maxEval = 10000;
  const double reqAbsError = 0.;
  // const double reqRelError = 0.003;
  const double reqRelError =
      _target_accuracy;  // use same accuracy as the model

  // a = x[0]
  // b <-> x[1] mapped to [-intfy, +infty]
  const double xmin[2] = {+0. + std::numeric_limits<float>::epsilon(),
                          -1. + std::numeric_limits<float>::epsilon()};
  const double xmax[2] = {+1. - std::numeric_limits<float>::epsilon(),
                          +1. - std::numeric_limits<float>::epsilon()};
  double val[1], err[1];

  ABCdata fdata = {delta, _epsilon, _xi, _omega};

  hcubature(1, ab_integrand, &fdata, 2, xmin, xmax, maxEval, reqAbsError,
            reqRelError, ERROR_INDIVIDUAL, val, err);

  if (!std::isfinite(val[0])) {
    std::cerr << "#pdf_delta__mu: problem for delta[last] = " << delta.back()
              << std::endl;
    return 0.;
  }
  return val[0];
}

double ABCNumericModel::pdf_delta__mu(const double& delta_next) const {
  std::vector<double> delta = _delta;
  delta.push_back(delta_next);
  return pdf_delta__mu(delta);
}

}  // namespace miho
