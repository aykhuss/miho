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

  // //> new numerical integration
  double generic = nintegrate_1D(
      [&](double a) -> double { return a_integrand(a, delta); }, 0., 1.);
  // return generic;

  // //> old numerical implementation
  double numerical = a_nint(delta);

  //> analytic result for xi==1
  double result = 0.;
  auto m = delta.size() - 1;  // counting starts at 0
  double mpeps = double(m) + _epsilon;
  std::vector<double> delta_neg = delta;
  std::vector<double> delta_pos = delta;
  for (auto i = 0; i < delta.size(); ++i) {
    if (delta.at(i) < 0.) {
      delta_pos.at(i) = 0.;
    } else {
      delta_neg.at(i) = 0.;
    }
  }
  std::vector<double> a_neg = new_a_list(delta_neg);
  std::vector<double> a_pos = new_a_list(delta_pos);

  // fmt::print("---[delta:]\n");
  // for (auto i = 0; i < delta.size(); ++i) {
  //   fmt::print("{} ", delta.at(i));
  // }
  // fmt::print("\n---\n");
  // fmt::print("---[a-lists:]\n");
  // for (auto i = 0; i < a_neg.size(); ++i) {
  //   fmt::print("{:10d} ", i);
  // }
  // fmt::print("\n");
  // for (auto i = 0; i < a_neg.size(); ++i) {
  //   fmt::print("{} ", a_neg.at(i));
  // }
  // fmt::print("\n");
  // for (auto i = 0; i < a_pos.size(); ++i) {
  //   fmt::print("{} ", a_pos.at(i));
  // }
  // fmt::print("\n---\n");

  double acc = 0.;

  for (auto i_neg = 0; i_neg <= m; ++i_neg) {
    double a_neg_low = std::min(1., a_neg[i_neg + 1]);
    double a_neg_upp = std::min(1., a_neg[i_neg]);
    if (is_approx(a_neg_low, a_neg_upp)) continue;
    // fmt::print("[-] {}: [{},{}]\n", i_neg, a_neg_low, a_neg_upp);

    for (auto i_pos = 0; i_pos <= m; ++i_pos) {
      double a_pos_low = std::min(1., a_pos[i_pos + 1]);
      double a_pos_upp = std::min(1., a_pos[i_pos]);
      if (is_approx(a_pos_low, a_pos_upp)) continue;
      // fmt::print("[+] {}: [{},{}]\n", i_pos, a_pos_low, a_pos_upp);

      double a_low = std::max(a_neg_low, a_pos_low);
      double a_upp = std::min(a_neg_upp, a_pos_upp);
      if (is_approx(a_low, a_upp) || a_low > a_upp) continue;
      // we found an a-segment: `[a_low,a_upp]` where:
      //    Delta^- = delta[i_neg]/a^i_neg <= 0
      //    Delta^+ = delta[i_pos]/a^i_pos >= 0
      fmt::print("### ({},{}) : [{},{}]\n", i_neg, i_pos, a_low, a_upp);
      // for stability:
      a_low = std::max(a_low, std::numeric_limits<double>::epsilon());

      // we need to check if [a_low,a_upp] still needs to be sub-divided
      // by checking if the crossing `-1+Delta^+ == 1+Delta^-`
      // occurs within this segment
      auto func_xing = [=](double a) -> double {
        return delta[i_pos] / std::pow(a, i_pos) -
               delta[i_neg] / std::pow(a, i_neg) - 2.;
      };
      double a_xing;
      bool has_xing = get_root(func_xing, a_low, a_upp, a_xing);
      fmt::print("# > xing? {}: {}\t\t\t----------------\n", has_xing,
      a_xing);

      for (auto j = 0; j <= _omega; ++j) {
        double pre_j = std::pow(-1, j) * binomial(_omega, j);

        // first piece  [a_low,a_xing]  (a_xing==a_upp if no xing)
        if (!has_xing) a_xing = a_upp;
        if (func_xing(a_low) >= 0.) {
          // `-1+Delta^+ >= 1+Delta^-`
          acc += pre_j * std::pow(2., 1 + mpeps) / mpeps *
                 a_nintegral(a_low, a_xing, std::fabs(delta[i_neg]),
                            std::fabs(delta[i_pos]), _epsilon,
                            j - m * (m + 1) / 2, i_neg, i_pos, m);
        } else {
          // `-1+Delta^+ < 1+Delta^-`
          acc += pre_j * (2. / mpeps + 2. -
                          a_nintegral2(a_low, a_xing, std::fabs(delta[i_pos]),
                                      j - m * (m + 1) / 2, i_pos) +
                          a_nintegral2(a_low, a_xing, std::fabs(delta[i_neg]),
                                      j - m * (m + 1) / 2, i_neg));
        }
        // second piece  [a_xing,a_upp] if it exists
        if (has_xing) {
          if (func_xing(a_upp) >= 0.) {
            // `-1+Delta^+ >= 1+Delta^-`
            acc += pre_j * std::pow(2., 1 + mpeps) / mpeps *
                   a_nintegral(a_xing, a_upp, std::fabs(delta[i_neg]),
                              std::fabs(delta[i_pos]), _epsilon,
                              j - m * (m + 1) / 2, i_neg, i_pos, m);
          } else {
            // `-1+Delta^+ < 1+Delta^-`
            acc += pre_j * (2. / mpeps + 2. -
                            a_nintegral2(a_xing, a_upp, std::fabs(delta[i_pos]),
                                        j - m * (m + 1) / 2, i_pos) +
                            a_nintegral2(a_xing, a_upp, std::fabs(delta[i_neg]),
                                        j - m * (m + 1) / 2, i_neg));
          }
        }  // if has_xing
      }    // for j
    }      // for i_pos
  }        // for i_neg
  // std::cin.ignore();

  result = acc * _epsilon * (1. + _omega) / std::pow(2., m + 1) / (1 + mpeps);

  fmt::print("numerical: {}\n", numerical);
  fmt::print("generic:   {}  |  {}\n", generic, numerical / generic);
  fmt::print("analytic:  {}  |  {}\n", result, numerical / result);
  std::cin.ignore();

  return result;
}

double ABCModel::pdf_delta__mu(const double& delta_next) const {
  std::vector<double> delta = _delta;
  delta.push_back(delta_next);
  return pdf_delta__mu(delta);
}

double ABCModel::a_integral(const double& a_lower, const double& a_upper,
                            const double& alpha, const double& beta,
                            const double& epsilon, int l, int i, int j, int m) {
  // fmt::print("--- a_integral:\n");
  // fmt::print("a_lower: {}\n", a_lower);
  // fmt::print("a_upper: {}\n", a_upper);
  // fmt::print("alpha: {}\n", alpha);
  // fmt::print("beta: {}\n", beta);
  // fmt::print("epsilon: {}\n", epsilon);
  // fmt::print("l: {}\n", l);
  // fmt::print("i: {}\n", i);
  // fmt::print("j: {}\n", j);
  // fmt::print("m: {}\n", m);

  double mpeps = double(m) + epsilon;

  if ((i == j) || is_approx(alpha, 0.) || is_approx(beta, 0.)) {
    double exponent = 1. + l + i * mpeps;
    if (is_approx(exponent, 0.)) {
      // fmt::print("a_integral: i==j & exp = 0\n");
      return std::pow(alpha + beta, -mpeps) * std::log(a_upper / a_lower);
    } else {
      // fmt::print("a_integral: i==j & exp != 0\n");
      return std::pow(alpha + beta, -mpeps) / exponent *
             (std::pow(a_upper, exponent) - std::pow(a_lower, exponent));
    }
  }
  // i!=j
  if (i < j) {
    return a_integral(a_lower, a_upper, beta, alpha, epsilon, l, j, i, m);
  }
  // i>j
  // fmt::print("a_integral: default\n");
  return std::pow(alpha, -mpeps) * std::pow(a_upper, l + i * mpeps) *
         a_master(l + i * mpeps, -mpeps, i - j, a_lower / a_upper,
                  std::pow(a_upper, i - j) * beta / alpha);
}

double ABCModel::a_master(const double& p, const double& q, const double& r,
                          const double& A, const double& B) {
  double opB = 1. + B;
  double omAr = (1. - std::pow(A, r));
  return std::pow(opB, q) * omAr / r *
         appell_F1(1., (r - p - 1.) / r, -q, 2., omAr, B * omAr / opB);
  // return std::pow(opB, q) * omAr / r *
  //        rf1_adsj(1., (r - p - 1.) / r, -q, 2., omAr, B * omAr / opB);
}

double ABCModel::a_integral2(const double& a_lower, const double& a_upper,
                             const double& alpha, int l, int i) {
  double exponent = 1. + l - i;
  if (is_approx(exponent, 0.)) {
    // fmt::print("a_integral2: exp = 0\n");
    return alpha * std::log(a_upper / a_lower);
  } else {
    // fmt::print("a_integral2: exp != 0\n");
    return alpha / exponent *
           (std::pow(a_upper, exponent) - std::pow(a_lower, exponent));
  }
}

double ABCModel::a_integrand(const double& a,
                             const std::vector<double>& delta) const {
  double a_val = a;
  if (a_val < 1e3 * std::numeric_limits<double>::epsilon()) {
    a_val = 1e3 * std::numeric_limits<double>::epsilon();
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

double ABCModel::a_nint(const std::vector<double>& delta) const {
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

double ABCModel::a_nintegral(const double& a_lower, const double& a_upper,
                             const double& alpha, const double& beta,
                             const double& epsilon, int l, int i, int j,
                             int m) {
  // double analytic =
  //     a_integral(a_lower, a_upper, alpha, beta, epsilon, l, i, j, m);

  auto func_gen = [=](double a) -> double {
    return std::pow(a, l) *
           std::pow(alpha / std::pow(a, i) + beta / std::pow(a, j),
                    -(double(m) + epsilon));
  };
  double result = nintegrate_1D(func_gen, a_lower, a_upper);

  // if (std::fabs((analytic - result) / (analytic + result)) > 1e-2) {
  //   fmt::print("--- a_nintegral:\n");
  //   fmt::print("a_lower: {}\n", a_lower);
  //   fmt::print("a_upper: {}\n", a_upper);
  //   fmt::print("alpha: {}\n", alpha);
  //   fmt::print("beta: {}\n", beta);
  //   fmt::print("epsilon: {}\n", epsilon);
  //   fmt::print("l: {}\n", l);
  //   fmt::print("i: {}\n", i);
  //   fmt::print("j: {}\n", j);
  //   fmt::print("m: {}\n", m);
  //   fmt::print("---\n");
  //   fmt::print("numerical: {}\n", result);
  //   fmt::print("analytic:  {}  |  {}\n", analytic, result / analytic);
  //   fmt::print(
  //       "NIntegrate[a^({}) ({}/a^{}+{}/a^{})^(-({}+{})),{{a, {}, {}}}]\n", l,
  //       alpha, i, beta, j, m, epsilon, a_lower, a_upper);
  //   // std::cin.ignore();
  // }
  return result;
}

double ABCModel::a_nintegral2(const double& a_lower, const double& a_upper,
                              const double& alpha, int l, int i) {
  // double analytic = a_integral2(a_lower, a_upper, alpha, l, i);

  auto func_gen = [=](double a) -> double {
    return std::pow(a, l) * alpha / std::pow(a, i);
  };
  double result = nintegrate_1D(func_gen, a_lower, a_upper);

  // if (std::fabs((analytic - result) / (analytic + result)) > 1e-2) {
  //   fmt::print("--- a_nintegral2:\n");
  //   fmt::print("a_lower: {}\n", a_lower);
  //   fmt::print("a_upper: {}\n", a_upper);
  //   fmt::print("alpha: {}\n", alpha);
  //   fmt::print("l: {}\n", l);
  //   fmt::print("i: {}\n", i);
  //   fmt::print("---\n");
  //   fmt::print("numerical: {}\n", result);
  //   fmt::print("analytic:  {}  |  {}\n", analytic, result / analytic);
  //   fmt::print("NIntegrate[a^({}) * {}/a^{},{{a, {}, {}}}]\n", l, alpha, i,
  //              a_lower, a_upper);
  //   // std::cin.ignore();
  // }
  return result;
}

}  // namespace miho
