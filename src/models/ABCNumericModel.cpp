#include "ABCNumericModel.h"

#include <cubature.h>
#include <fmt/format.h>

#include <algorithm>
// #include <boost/math/quadrature/naive_monte_carlo.hpp>
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
  /// some a-factors combined with max_val for stability
  result *= std::pow(1. - a, abc->omega) * std::pow(a, m * abc->epsilon / 2.) *
            abc->epsilon * (1. + abc->omega) /
            (std::pow(2., m + 1.) * abc->xi * (m + 1. + abc->epsilon));
  /// find max value (absorbs factor of a^(m/2) for numerical stability)
  double a_fac = std::pow(a, m / 2.);
  double b_fac = a_fac * b;
  double max_val = a_fac * std::max(1., std::fabs(b) / abc->xi);
  double del_fac = a_fac;  // accumulate a-powers
  for (auto i = 1; i < abc->delta.size(); ++i) {
    del_fac /= a;  // del_fac = a^(m/2) / a^i
    double dtest = std::fabs(abc->delta.at(i) * del_fac - b_fac);
    if (dtest > max_val) max_val = dtest;
  }
  result /= std::pow(max_val, m + 1. + abc->epsilon);
  /// done.
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
  //----- cubature library
  size_t maxEval = 100000;
  double reqAbsError = 0.;
  double reqRelError = std::min(_nint_rel_err, _target_accuracy);
  // reqRelError *= 5.;  // init value (n_loop > 1)
  // a = x[0]
  // b <-> x[1] mapped to [-intfy, +infty]
  const double xmin[2] = {+0. + std::numeric_limits<float>::epsilon(),
                          -1. + std::numeric_limits<float>::epsilon()};
  const double xmax[2] = {+1. - std::numeric_limits<float>::epsilon(),
                          +1. - std::numeric_limits<float>::epsilon()};
  double val[1], err[1];
  /// package up struct
  ABCdata fdata = {delta, _epsilon, _xi, _omega};

  double val_last = 0.;
  double err_last = -1.;
  int count = 0;
  while (true) {
    count++;
    hcubature(1, ab_integrand, &fdata, 2, xmin, xmax, maxEval, reqAbsError,
              reqRelError, ERROR_INDIVIDUAL, val, err);

    if (!std::isfinite(val[0])) {
      std::cerr << "pdf_delta__mu: problem for delta[last] = " << delta.back()
                << std::endl;
      return 0.;
    }

    if (err_last > 0. && is_approx(val[0] + val_last, 0.)) break;
    // if (err_last > 0. &&
    //     std::fabs((val[0] - val_last) / (val[0] + val_last)) < _nint_rel_err) {
    //   break;
    // }
    if ((err_last > 0.) &&
        ((std::fabs(val[0] - val_last) <
          3. * std::sqrt(err[0] * err[0] + err_last * err_last)) ||
         (std::fabs((val[0] - val_last) / (val[0] + val_last)) <
          _target_accuracy))) {
      break;
    }
    if (count >= 7) {
      // std::cerr << "pdf_delta__mu: couldn't find plateau in " << count
      //           << " steps  (acc = " << _nint_rel_err << " -> " << reqRelError
      //           << ")" << std::endl;
      // std::cerr << "last: " << val_last << " +- " << err_last << std::endl;
      // std::cerr << "curr: " << val[0] << " +- " << err[0] << std::endl;
      break;
    }

    reqRelError /= 3.;
    val_last = val[0];
    err_last = err[0];
  }
  return val[0];

  // hcubature(1, ab_integrand, &fdata, 2, xmin, xmax, maxEval, reqAbsError,
  //           reqRelError, ERROR_INDIVIDUAL, val, err);
  // if (!std::isfinite(val[0])) {
  //   std::cerr << "#pdf_delta__mu: problem for delta[last] = " << delta.back()
  //             << std::endl;
  //   return 0.;
  // }
  // return val[0];

  // //----- boost naive MC
  // using boost::math::quadrature::naive_monte_carlo;
  // auto f_ab = [&](std::vector<double> const& x) -> double {
  //   double result = 1.;
  //   /// a mapping:  linear
  //   double a = x[0];
  //   /// b-mapping:  arctanh(x)
  //   const double b = std::atanh(x[1]);
  //   result /= (1. - x[1]) * (1. + x[1]);  // Jac for b-mapping
  //   // std::cerr << "result[1] = " << result << std::endl;
  //   double m = delta.size() - 1.;  // counting starts at 0
  //   result *= std::pow(1. - a, _omega) * std::pow(a, m * _epsilon / 2.) *
  //             _epsilon * (1. + _omega) /
  //             (std::pow(2., m + 1.) * _xi * (m + 1. + _epsilon));
  //   // std::cerr << "result[3] = " << result << std::endl;
  //   //> b-dependent piece (max function)
  //   double max_val = std::max(1., std::fabs(b) / _xi);
  //   double pow_a = 1.;  // accumulate a-powers
  //   for (auto i = 1; i < delta.size(); ++i) {
  //     pow_a *= a;  // pow_a = a^i
  //     double dtest = std::fabs(delta.at(i) / pow_a - b);
  //     if (dtest > max_val) max_val = dtest;
  //   }
  //   result /= std::pow(std::sqrt(pow_a) * max_val, m + 1. + _epsilon);
  //   // std::cerr << "result[4] = " << result << std::endl;
  //   /// done.
  //   // std::cerr << "eps = " << _epsilon << std::endl;
  //   // std::cerr << "xi  = " << _xi << std::endl;
  //   // std::cerr << "omg = " << _omega << std::endl;
  //   // std::cerr << "f_ab(" << a << "," << b << ") = " << result <<
  //   std::endl;
  //   // std::cin.ignore();
  //   if (!std::isfinite(result)) {
  //     std::cerr << "#f_ab: problem for a = " << a << ", b = " << b <<
  //     std::endl; std::cerr << "m       = " << m << std::endl; std::cerr <<
  //     "eps     = " << _epsilon << std::endl; std::cerr << "xi      = " << _xi
  //     << std::endl; std::cerr << "omg     = " << _omega << std::endl;
  //     std::cerr << "max_val = " << max_val << std::endl;
  //     std::cerr << "pow_a   = " << pow_a << std::endl;
  //     std::cin.ignore();
  //     return 0.;
  //   }
  //   return result;
  // };
  // std::vector<std::pair<double, double>> bounds{{0, 1}, {-1, +1}};
  // double error_goal = 0.05;
  // boost::math::quadrature::naive_monte_carlo<double, decltype(f_ab)> mc(
  //     f_ab, bounds, error_goal, /*singular = */ true, /*threads = */ 1,
  //     /*seed = */ 0);
  // std::future<double> task = mc.integrate();
  // double res_mc = task.get();

  // // std::cerr << "cubature: " << val[0] << std::endl;
  // // std::cerr << "naive mc: " << res_mc << std::endl;
  // // std::cin.ignore();

  // return res_mc;
}

double ABCNumericModel::pdf_delta__mu(const double& delta_next) const {
  std::vector<double> delta = _delta;
  delta.push_back(delta_next);
  return pdf_delta__mu(delta);
}

}  // namespace miho
