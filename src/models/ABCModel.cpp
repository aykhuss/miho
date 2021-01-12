#include "ABCModel.h"

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace {
// for a monotonic function, tell me if there's a sign change
//    => a root in the interval
// iteratively find that root
bool get_root(std::function<double(double)> func, const double& x_low,
              const double& x_upp, double& x_root) {
  // settings for root finding & termination
  const double target_accuracy = std::numeric_limits<double>::epsilon();
  const size_t max_it = 1000;
  x_root = std::numeric_limits<double>::infinity();  // no root
  double x0 = x_low;
  double y0 = func(x0);
  double x1 = x_upp;
  double y1 = func(x1);
  // has root?
  if (miho::sgn(y0) * miho::sgn(y1) >= 0) {
    // same sign @ edges (monotonic) => no root
    return false;
  }
  int it = 0;
  while (std::fabs((x1 - x0) / (x_upp - x_low)) > target_accuracy) {
    ++it;
    if (it >= max_it) {
      throw "get_root: failed to converge?!";
      return true;
    }
    // always divide by two
    x_root = (x0 + x1) / 2.;
    double y_root = func(x_root);
    if (miho::sgn(y0) * miho::sgn(y_root) >= 0) {
      // continue in [y_root,y1]
      x0 = x_root;
      y0 = y_root;
    } else {
      // continue in [y0,y_root]
      x1 = x_root;
      y1 = y_root;
    }
  }
  return true;
}
}  // namespace

namespace miho {

double ABCModel::pdf(const double& val) const {
  // Eq.(3.14) : assumes j == 1
  return pdf_delta___delta_mu(delta_next(val)) / std::fabs(_sigma.front());
}

double ABCModel::pdf_delta___delta_mu(const double& delta_next) const {
  // Eq.(4.10) : assumes j == 1
  if (!_q_pdf_den) {
    _pdf_den = pdf_delta__mu(_delta);
    _q_pdf_den = true;
  }
  double pdf_num = pdf_delta__mu(delta_next);
  return pdf_num / _pdf_den;
}

double ABCModel::pdf_delta__mu(const std::vector<double>& delta) const {
  auto m = delta.size() - 1;  // counting starts at 0

  /// delta with sign flipped for odd entries
  std::vector<double> delta_tilde(delta.size());
  std::transform(delta.begin(), delta.end(), delta_tilde.begin(),
                 [idx = 0](const double& d_idx) mutable {
                   return idx++ % 2 == 0 ? d_idx : -d_idx;
                 });

  double acc = 0.;
  for (auto j = 0; j <= _omega; ++j) {
    double acc_j = j % 2 == 0 ? 1. : -1.;
    acc_j *= binomial(_omega, j);
    acc_j *= II(j - m * (m + 1) / 2, m + 3 + _epsilon, delta) +
             II(j - m * (m + 1) / 2, m + 3 + _epsilon, delta_tilde);
    acc += acc_j;
    // std::cout << "> " << j << ": " << acc_j << "\n";
  }
  // std::cin.ignore();

  return acc * (1. + _omega) * _epsilon * std::pow(_eta, _epsilon) / _xi /
         std::pow(2., m + 3);
}

double ABCModel::pdf_delta__mu(const double& delta_next) const {
  std::vector<double> delta = _delta;
  delta.push_back(delta_next);
  return pdf_delta__mu(delta);
}

double ABCModel::II(int s, const double& alpha,
                    const std::vector<double>& delta) const {
  using Key_t = std::pair<int, int>;
  using Val_t = std::pair<double, double>;

  // std::cout << "ABCModel::II:  \n";
  // std::cout << "s     = " << s << std::endl;
  // std::cout << "alpha = " << alpha << std::endl;
  // std::cout << "delta = ";
  // for (const auto& d_i : delta) std::cout << d_i << " ";
  // std::cout << std::endl;

  /// \int_{a0}^{a1} da a^(ex)
  auto int_a_pow = [](const double& ex, const double& a0, const double& a1) {
    return is_approx(ex, -1)
               ? std::log(a1 / a0)
               : (std::pow(a1, ex + 1) - std::pow(a0, ex + 1)) / (ex + 1);
  };

  double acc = 0.;
  auto parts = min_max_partitions(delta);
  for (const std::pair<Key_t, Val_t>& p : parts) {
    const Key_t& key = p.first;
    const Val_t& val = p.second;
    const double a_low = val.first < std::numeric_limits<double>::epsilon()
                             ? std::numeric_limits<double>::epsilon()
                             : val.first;
    const double a_upp = val.second > 1. ? 1. : val.second;
    // std::cout << "[" << a_low << "," << a_upp << "] "
    //           << (a_upp - a_low) / (a_upp + a_low) << "\n";

    /// (A): numeric
    auto func_A = [&s, &alpha, &delta, eta = _eta, p = key.first,
                   q = key.second](const double& a) -> double {
      double result = 0.;
      const double R = delta[p] / std::pow(a, p);
      const double r = delta[q] / std::pow(a, q);
      if ((R - r) / 2. < eta) {
        result = std::pow(eta, 2 - alpha) / (alpha - 2) -
                 (R - r) / 2. * std::pow(eta, 1 - alpha) / (alpha - 1);
      } else {
        result = std::pow((R - r) / 2., 2 - alpha) / (alpha - 2) / (alpha - 1);
      }
      return std::pow(a, s) * result;
    };
    const double num_A = nintegrate_1D(func_A, a_low, a_upp);
    /// (A): analytic

    acc += 2. * num_A;
    // std::cout << " (A) " << num_A << "\n";

    if (is_approx(_xi, 1.)) {
      /// (B): numeric
      auto func_B = [&s, &alpha, &delta, eta = _eta, p = key.first,
                     q = key.second](const double& a) -> double {
        double result = 0.;
        const double R = delta[p] / std::pow(a, p);
        const double r = delta[q] / std::pow(a, q);
        if (r > 0.) {
          if ((R - r) / 2. < eta) {
            result = std::pow(eta, 1 - alpha);
          } else {
            result = std::pow((R - r) / 2., 1 - alpha);
          }
          result *= r / (alpha - 1);
        }
        return std::pow(a, s) * result;
      };
      const double num_B = nintegrate_1D(func_B, a_low, a_upp);
      /// (B): analytic
      double ana_B = 0.;
      if (delta[key.second] > 0.) {
        ana_B = int_a_pow(s - key.second, a_low, a_upp) * delta[key.second] *
                 std::pow(_eta, 1 - alpha) / (alpha - 1);
      }
      // acc += -ana_B;
      // std::cout << " (B) " << num_B << " v.s. " << ana_B << ": "
      //           << (num_B - ana_B) / (num_B + ana_B) << " | " << s
      //           << ", " << key.second << "\n";
      acc += -num_B;

    } else {
      /// (C): numeric
      auto func_C = [&s, &alpha, &delta, eta = _eta, xi = _xi,
                     p = key.first](const double& a) -> double {
        double result = 0.;
        const double R__oMxi = delta[p] / std::pow(a, p) / (1 - xi);
        if (R__oMxi < eta) {
          result = std::pow(eta, 2 - alpha) / (alpha - 2) -
                   R__oMxi * std::pow(eta, 1 - alpha) / (alpha - 1);
        } else {
          result = std::pow(R__oMxi, 2 - alpha) / (alpha - 2) / (alpha - 1);
        }
        return std::pow(a, s) * result;
      };
      const double num_C = nintegrate_1D(func_C, a_low, a_upp);
      /// (C): analytic

      acc += -1. / (1 - _xi) * num_C;
      std::cout << " (C) " << num_C << "\n";
      /// (D)
      auto func_D = [&s, &alpha, &delta, eta = _eta, xi = _xi,
                     q = key.second](const double& a) -> double {
        double result = 0.;
        const double r__oMxi = delta[q] / std::pow(a, q) / (1 - xi);
        if (-r__oMxi < eta) {
          result = std::pow(eta, 2 - alpha) / (alpha - 2) +
                   r__oMxi * std::pow(eta, 1 - alpha) / (alpha - 1);
        } else {
          result = std::pow(-r__oMxi, 2 - alpha) / (alpha - 2) / (alpha - 1);
        }
        return std::pow(a, s) * result;
      };
      const double res_D = nintegrate_1D(func_D, a_low, a_upp);
      acc += -1. / (1 - _xi) * res_D;
      std::cout << " (D) " << res_D << "\n";
    }
  }

  // std::cout << " acc " << acc << "\n";
  // std::cin.ignore();

  return acc;
}

}  // namespace miho
