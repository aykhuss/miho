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
  std::vector<double> delta;
  const double& epsilon;
  const double& eta;
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
  if (miho::is_approx(a, 0.)) a = std::numeric_limits<double>::epsilon();

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
  result *= (1. + abc->omega) * abc->epsilon *
            std::pow(abc->eta, abc->epsilon) / std::pow(2., m + 3.) / abc->xi;
  /// a-integrand
  result *= std::pow(1. - std::fabs(a), abc->omega) /
            std::pow(std::fabs(a), m * (m + 1.) / 2.);
  /// c integration constants
  result /= (m + 2. + abc->epsilon);
  /// find max value (absorbs factor of a^(m/2) for numerical stability)
  double max_val = std::max(abc->eta, std::fabs(b) / abc->xi);
  double del_fac = 1.;  // accumulate a-powers
  for (auto i = 0; i < abc->delta.size(); ++i) {
    double dtest = std::fabs(abc->delta.at(i) * del_fac - b);
    del_fac /= a;  // del_fac = 1 / a^i on nxt iter
    if (dtest > max_val) max_val = dtest;
  }
  result /= std::pow(max_val, m + 2. + abc->epsilon);
  /// done.
  if (!std::isfinite(result)) {
    std::cerr << "#ab_integrand: problem for a = " << a << ", b = " << b
              << std::endl;
    return 0;
  }
  fval[0] = result;
  return 0;  // success
}

int a_integrand(unsigned ndim, const double* x, void* fdata, unsigned fdim,
                 double* fval) {
  if (ndim != 1) return 1;
  if (fdim != 1) return 1;
  fval[0] = 0.;
  double result = 1.;

  //> a mapping
  double a = x[0];
  // // for numerical stability
  if (miho::is_approx(a, 0.)) a = std::numeric_limits<double>::epsilon();

  const ABCdata* abc = (ABCdata*)fdata;
  double m = abc->delta.size() - 1.;  // counting starts at 0

  /// find R & r
  double R = abc->delta[0];  // max
  double r = abc->delta[0];  // min
  double del_fac = 1.;  // accumulate a-powers
  for (auto i = 1; i < abc->delta.size(); ++i) {
    del_fac /= a;  // del_fac = 1 / a^i on nxt iter
    double dtest = abc->delta[i] * del_fac;
    if (dtest > R) R = dtest;
    if (dtest < r) r = dtest;
  }

  /// branch by xi
  /// m + 1 + abc->epsilon = 0 not relevant, as it would be non-integrable
  double acc = 0.;
  if (miho::is_approx(abc->xi, 1)) {
    R = R > 0 ? R : 0;
    r = r < 0 ? r : 0;
    const double b_cusp = (r + R) / 2.;
    const double b_1upp = 1 + r;
    const double b_1low = R - 1;
    if (b_1low > b_1upp) {
      acc = std::pow((R - r) / 2., -m - 1 - abc->epsilon) * 2. /
            (m + 1 + abc->epsilon);
    } else {
      acc = 2. / (m + 1 + abc->epsilon) + 2. - (R - r);
    }
  } else if (abc->xi < 1) {
    /*
     *  interval:  func
     *  0:         b / xi
     *  1:         (b - r)
     *  2:         1
     *  3:         (R - b)
     *  4:         -b / xi
    */
    const double mpope = m + 1 + abc->epsilon;
    double b_upp = std::numeric_limits<double>::infinity();
    double b_next;
    int i_next = -1;  // the next interval
    auto set_next = [&b_next, &i_next](int i_val, const double& b_val) {
      if ((i_next < 0) || (b_val > b_next)) {
        b_next = b_val;
        i_next = i_val;
      }
    };
    /// interval: 0
    set_next(1, (-r) * abc->xi / (1 - abc->xi));
    set_next(2, abc->xi);
    set_next(3, R * abc->xi / (1 + abc->xi));
    set_next(4, 0);
    acc += std::pow(b_next / abc->xi, -mpope) * abc->xi / mpope;
    if (i_next == 1) {
      b_upp = b_next;
      i_next = -1;
      set_next(2, 1 + r);
      set_next(3, (R + r) / 2.);
      set_next(4, abc->xi / (1 + abc->xi) * r);
      acc +=
          (std::pow(b_next - r, -mpope) - std::pow(b_upp - r, -mpope)) / mpope;
    }
    if (i_next == 2) {
      b_upp = b_next;
      i_next = -1;
      set_next(3, R - 1);
      set_next(4, -abc->xi);
      acc += b_upp - b_next;
    }
    if (i_next == 3) {
      b_upp = b_next;
      i_next = -1;
      set_next(4, -R * abc->xi / (1 - abc->xi));
      acc +=
          (std::pow(R - b_upp, -mpope) - std::pow(R - b_next, -mpope)) / mpope;
    }
    if (i_next != 4) throw "#a_integrand: xi < 1: no i_next closure";
    b_upp = b_next;
    acc += std::pow(-b_upp / abc->xi, -mpope) * abc->xi / mpope;

  } else {
    /*
     *  interval:  func
     *  0:         (b - r)
     *  1:         b / xi
     *  2:         1
     *  3:         -b / xi
     *  4:         (R - b)
    */
    // std::cout << "\n";
    // std::cout << "r = " << r << "\n";
    // std::cout << "R = " << R << "\n";
    const double mpope = m + 1 + abc->epsilon;
    double b_upp = std::numeric_limits<double>::infinity();
    double b_next;
    int i_next = -1;  // the next interval
    auto set_next = [&b_next, &i_next](int i_val, const double& b_val) {
      // std::cout << " > " << i_val << ": " << b_val << "\n";
      if ((i_next < 0) || (b_val > b_next)) {
        b_next = b_val;
        i_next = i_val;
      }
    };
    /// interval: 0
    set_next(1, (-r) * abc->xi / (1 - abc->xi));
    set_next(2, 1+r);
    set_next(3, abc->xi / (1 + abc->xi) * r);
    set_next(4, (R + r) / 2.);
    acc += std::pow(b_next - r, -mpope) / mpope;
    // std::cout << "interval 0: " << b_next<< " | " << acc << "\n";
    if (i_next == 1) {
      b_upp = b_next;
      i_next = -1;
      set_next(2, abc->xi);
      set_next(3, 0.);
      set_next(4, R * abc->xi / (1 + abc->xi));
      acc += (std::pow(b_next / abc->xi, -mpope) -
              std::pow(b_upp / abc->xi, -mpope)) *
             abc->xi / mpope;
      // std::cout << "interval 1: " << b_next<< " | " << acc << "\n";
    }
    if (i_next == 2) {
      b_upp = b_next;
      i_next = -1;
      set_next(3, -abc->xi);
      set_next(4, R - 1);
      acc += b_upp - b_next;
      // std::cout << "interval 2: " << b_next<< " | " << acc << "\n";
    }
    if (i_next == 3) {
      b_upp = b_next;
      i_next = -1;
      set_next(4, -R * abc->xi / (1 - abc->xi));
      acc += (std::pow(-b_upp / abc->xi, -mpope) -
              std::pow(-b_next / abc->xi, -mpope)) *
             abc->xi / mpope;
      // std::cout << "interval 3: " << b_next<< " | " << acc << "\n";
    }
    if (i_next != 4) throw "#a_integrand: xi > 1: no i_next closure";
    b_upp = b_next;
    acc += std::pow(R - b_upp, -mpope) / mpope;
    // std::cout << "interval 4: " << b_next<< " | " << acc << "\n";
    // std::cin.ignore();

  }

  /// combine to final result
  result = acc * std::pow(1. - std::fabs(a), abc->omega) /
           std::pow(std::fabs(a), m * (m + 1.) / 2.);

  /// done.
  if (!std::isfinite(result)) {
    std::cerr << "#a_integrand: problem for a = " << a << std::endl;
    // std::cin.ignore();
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
  if (!_q_pdf_den) {
    _pdf_den = pdf_delta__mu(_delta);
    _q_pdf_den = true;
  }
  double pdf_num = pdf_delta__mu(delta_next);
  return pdf_num / _pdf_den;
}

double ABCNumericModel::pdf_delta__mu(const std::vector<double>& delta) const {
  //----- cubature library
  size_t maxEval = 100000;
  double reqAbsError = 0.;
  double reqRelError = std::min(_nint_rel_err, _target_accuracy);
  // a = x[0] from [-1,1]
  const double xmin[1] = {-1.};
  const double xmax[1] = {+1.};
  double val[1], err[1];
  /// package up struct
  ABCdata fdata = {delta, _epsilon, _eta, _xi, _omega};
  if (!is_approx(_eta, 1)) {
    // std::cout << "rescaling!\n";
    // std::cout << "before: ";
    // for (const auto& d : fdata.delta) std::cout << d << " ";
    // std::cout << "\n";
    std::transform(fdata.delta.begin(), fdata.delta.end(), fdata.delta.begin(),
                   [eta = _eta](const double& d) { return d / eta; });
    // std::cout << "after: ";
    // for (const auto& d : fdata.delta) std::cout << d << " ";
    // std::cout << "\n";
    // std::cin.ignore();
  }

  hcubature(1, a_integrand, &fdata, 1, xmin, xmax, maxEval, reqAbsError,
            reqRelError, ERROR_INDIVIDUAL, val, err);

  if (!std::isfinite(val[0])) {
    std::cerr << "pdf_delta__mu: problem for delta[last] = " << delta.back()
              << std::endl;
    return 0.;
  }

  const auto m = delta.size() - 1;
  return val[0] / std::pow(_eta, m + 1) * (1 + _omega) * _epsilon / _xi /
         std::pow(2, m + 3) / (m + 2 + _epsilon);
}

double ABCNumericModel::pdf_delta__mu__2D(const std::vector<double>& delta) const {
  //----- cubature library
  size_t maxEval = 100000;
  double reqAbsError = 0.;
  double reqRelError = std::min(_nint_rel_err, _target_accuracy);
  // reqRelError *= 5.;  // init value (n_loop > 1)
  // a = x[0] from [-1,1]
  // b <-> x[1] from [-1,1] mapped to [-intfy, +infty]
  const double xmin[2] = {-1. + std::numeric_limits<float>::epsilon(),
                          -1. + std::numeric_limits<float>::epsilon()};
  const double xmax[2] = {+1. - std::numeric_limits<float>::epsilon(),
                          +1. - std::numeric_limits<float>::epsilon()};
  double val[1], err[1];
  /// package up struct
  ABCdata fdata = {delta, _epsilon, _eta, _xi, _omega};
  // if (!is_approx(_eta, 1.)) {
  //   std::transform(fdata.delta.begin(), fdata.delta.end(), fdata.delta.begin(),
  //                  [fdata](const double& del) { return del / fdata.eta; });
  // }

  double val_last = 0.;
  double err_last = -1.;
  int count = 0;
  while (true) {
    count++;
    hcubature(1, ab_integrand, &fdata, 2, xmin, xmax, maxEval, reqAbsError,
              reqRelError, ERROR_INDIVIDUAL, val, err);

    if (!std::isfinite(val[0])) {
      std::cerr << "pdf_delta__mu__2D: problem for delta[last] = " << delta.back()
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
      // std::cerr << "pdf_delta__mu__2D: couldn't find plateau in " << count
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

  // return res_mc;
}

double ABCNumericModel::pdf_delta__mu(const double& delta_next) const {
  std::vector<double> delta = _delta;
  delta.push_back(delta_next);
  return pdf_delta__mu(delta);
}

}  // namespace miho
