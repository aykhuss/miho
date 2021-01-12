#include "Util.h"

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <list>
#include <random>

namespace {

/// random numbers
std::mt19937_64 mt{};
std::uniform_real_distribution<double> dist(0., 1.);

}  // namespace

namespace miho {

template <typename T>
int sgn(const T& val) {
  return (T(0) < val) - (val < T(0));
}
// Explicit instantiations.
template int sgn<int>(const int&);
template int sgn<double>(const double&);

bool is_approx(const double& x, const double& y) {
  /// check the numbers are finite
  if ( !std::isfinite(x) || !std::isfinite(y) ) return false;
  // 1.0 is here because epsilon is the smallest difference
  // that can be distinguished from 1.0
  double max_val = std::max({1.0, std::fabs(x), std::fabs(y)});
  return std::fabs(x - y) <= std::numeric_limits<double>::epsilon() * max_val;
}

bool is_approx_int(const double& x) { return is_approx(std::round(x), x); }

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
  const double target_accuracy = std::numeric_limits<double>::epsilon();
  // const double target_accuracy = 1e-11;
  const size_t max_it = 1000;

  // fmt::print("2F1({},{},{};{})\n", a, b, c, z);

  // catch special cases
  if (is_approx(b, c)) return std::pow(1. - z, -a);

  if (std::fabs(z) >= 1.) {
    fmt::print("hypergeometric_2F1: won't converge\n");
    throw "hypergeometric_2F1: won't converge";
  }
  // if (std::fabs(a) > std::fabs(b)) {
  //   // canonical ordering: |a| always smaller than |b|
  //   return hypergeometric_2F1(b, a, c, z);
  // }
  auto test_non_pos_int = [](double x) -> bool {
    return is_approx_int(x) && x < 0.5;
  };
  // if (z > 0.5 && !test_non_pos_int(c - a - b) &&
  //     !test_non_pos_int(a + b - c + 1.) && !test_non_pos_int(a + b - c) &&
  //     !test_non_pos_int(c - a - b + 1.)) {
  if (z > 0.5 && !test_non_pos_int(a + b - c + 1.) &&
      !test_non_pos_int(c - a - b + 1.)) {
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
    // fmt::print("hypergeometric_2F1: {} <-> {}  |  {} {}\n", it, dres, acc,
    //   std::fabs(dres / acc));
    if (!std::isfinite(dres)) {
      fmt::print("hypergeometric_2F1: failed to converge!\n");
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
  const size_t max_it = 1000;

  // fmt::print("appell_F1:\n");
  // fmt::print("a = {}\n", a);
  // fmt::print("b1 = {}\n", b1);
  // fmt::print("b2 = {}\n", b2);
  // fmt::print("c = {}\n", c);
  // fmt::print("z1 = {}\n", z1);
  // fmt::print("z2 = {}\n", z2);

  // catch special cases
  if (is_approx(b1, 0.)) return hypergeometric_2F1(a, b2, c, z2);
  if (is_approx(b2, 0.)) return hypergeometric_2F1(a, b1, c, z1);
  if (is_approx(z1, 0.)) return hypergeometric_2F1(a, b2, c, z2);
  if (is_approx(z2, 0.)) return hypergeometric_2F1(a, b1, c, z1);

  if (is_approx(z1, 1.) && is_approx(z2, 1.))
    return std::tgamma(c) * std::tgamma(c - a - b1 - b2) / std::tgamma(c - a) /
           std::tgamma(c - b1 - b2);
  if (is_approx(z2, 1.))
    return std::tgamma(c) * std::tgamma(c - a - b2) / std::tgamma(c - a) /
           std::tgamma(c - b2) * hypergeometric_2F1(a, b1, c - b2, z1);
  if (is_approx(z1, 1.))
    return std::tgamma(c) * std::tgamma(c - a - b1) / std::tgamma(c - a) /
           std::tgamma(c - b1) * hypergeometric_2F1(a, b2, c - b1, z2);

  if (std::fabs(z1) >= 1. || std::fabs(z2) >= 1.) {
    fmt::print("appell_F1: won't converge\n");
    throw "appell_F1: won't converge";
  }
  // if (std::fabs(z1) > std::fabs(z2)) {
  //   // canonical ordering: |z1| always smaller than |z2|
  //   return appell_F1(a, b2, b1, c, z2, z1);
  // }

  // map to smaller values of |z1| & |z2|
  double max_orig = std::max(std::fabs(z1), std::fabs(z2));
  double map_z1 = z1 / (z1 - 1.);
  double map_z2 = z2 / (z2 - 1.);
  double map1_z2 = (z1 - z2) / (z1 - 1.);
  double map2_z1 = (z2 - z1) / (z2 - 1.);
  std::array<double, 4> max_arg = {
      std::max(std::fabs(z1), std::fabs(z2)),
      std::max(std::fabs(map_z1), std::fabs(map_z2)),
      std::max(std::fabs(map_z1), std::fabs(map1_z2)),
      std::max(std::fabs(map2_z1), std::fabs(map_z2))};
  auto it = std::min_element(max_arg.begin(), max_arg.end());
  int idx = 0;
  if (it != max_arg.end()) idx = std::distance(max_arg.begin(), it);
  switch (idx) {
    case 1:
      return std::pow(1. - z1, -b1) * std::pow(1. - z2, -b2) *
             appell_F1(c - a, b1, b2, c, map_z1, map_z2);
    case 2:
      return std::pow(1. - z1, -a) *
             appell_F1(a, c - b1 - b2, b2, c, map_z1, map1_z2);
    case 3:
      return std::pow(1. - z2, -a) *
             appell_F1(a, b1, c - b1 - b2, c, map2_z1, map_z2);
  }

  double acc = 0.;
  for (auto it = 0; it < max_it; ++it) {
    double dres = 1.;
    for (auto i = 0; i < it; ++i) {
      dres *= (a + double(i)) * (b1 + double(i)) * (b2 + double(i)) *
              (c - a + double(i)) / (c + double(it - 1 + i)) /
              (c + double(2 * i)) / (c + double(2 * i + 1)) * (z1 * z2) /
              double(i + 1);
    }
    dres *= hypergeometric_2F1(a + it, b1 + it, c + 2 * it, z1);
    dres *= hypergeometric_2F1(a + it, b2 + it, c + 2 * it, z2);
    // fmt::print("{}: {} {}\n",it,dres,acc);
    if (!std::isfinite(dres)) {
      fmt::print("appell_F1: failed to converge!\n");
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
//         fmt::print("appell_F1: {} ({}, {}) <-> {}  |  {}\n", it, i1, i2,
//         dres,
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

double nintegrate_1D(std::function<double(double)> func, const double& x_low,
                     const double& x_upp) {
  struct Node {
    double x;
    double y;
  };
  std::list<Node> nodes;
  const double target_accuracy = 1e-5;
  const size_t max_nodes = 5000;

  double result = 0.;
  double error = 0.;

  if (is_approx(x_low, x_upp)) return 0.;
  // if (x_low > x_upp) return -nintegrate_1D(func, x_upp, x_low);
  if (x_low > x_upp) throw "nintegrate_1D: invalid range";

  // init with two equidistant bins
  Node nlow;
  nlow.x = (x_low + x_upp) / 2. - (x_upp - x_low) / 4.;
  nlow.y = func(nlow.x);
  nodes.emplace_front(nlow);
  Node nupp;
  nupp.x = (x_low + x_upp) / 2. + (x_upp - x_low) / 4.;
  nupp.y = func(nupp.x);
  nodes.emplace_back(nupp);

  while (nodes.size() < max_nodes) {
    // fmt::print("-: ({})\n",x_low);
    // size_t cnt = 0;
    // for (auto it = nodes.begin(); it != nodes.end(); ++it) {
    //   fmt::print("{}: ({}, {})\n", cnt, it->x, it->y);
    //   cnt++;
    // }
    // fmt::print("-: ({})\n",x_upp);
    // std::cin.ignore();

    // start with the `end` node (not part of the loop below)
    auto node_pos = nodes.end();
    result = (x_upp - std::prev(node_pos)->x) * std::prev(node_pos)->y;
    error = 0.;
    double curr_derivative = 0.;
    double last_derivative = 0.;
    double jump = std::fabs(result);
    // fmt::print("> {}: {}\n", nodes.size(), jump);
    // cnt = 0;

    for (auto it = nodes.begin(); it != nodes.end(); ++it) {
      double dres = 0.;
      if (it == nodes.begin()) {
        dres = (it->x - x_low) * it->y;
      } else {
        dres = (it->x - std::prev(it)->x) * (it->y + std::prev(it)->y) / 2.;
        curr_derivative =
            (it->y - std::prev(it)->y) / (it->x - std::prev(it)->x);
        error += std::fabs(curr_derivative - last_derivative) *
                 std::pow(it->x - std::prev(it)->x, 2) / 12.;
      }
      result += dres;
      dres = std::fabs(dres);
      // fmt::print("> {}: {}\n", cnt, dres);
      // cnt++;
      if (dres > jump) {
        jump = dres;
        node_pos = it;
        // fmt::print("> NEW node_pos!\n");
      }
      last_derivative = curr_derivative;
    }

    // check for accuracy termination condition
    if ((error / result) <= target_accuracy) {
      break;
    }

    Node n_new;
    if (node_pos == nodes.begin()) {
      n_new.x = (node_pos->x + x_low) / 2.;
    } else if (node_pos == nodes.end()) {
      n_new.x = (std::prev(node_pos)->x + x_upp) / 2.;
    } else {
      n_new.x = (node_pos->x + std::prev(node_pos)->x) / 2.;
    }
    n_new.y = func(n_new.x);
    node_pos = nodes.insert(node_pos, n_new);
    // fmt::print("# new element: ({},{}) \n", n_new.x, n_new.y);
  }

  return result;
}


std::vector<double> a_list(const std::vector<double>& delta, bool invert) {
  const auto n_delta = delta.size();
  std::vector<double> a(n_delta + 1, 0.);

  /// "guesses"
  a[0] = std::numeric_limits<double>::infinity();
  a[n_delta] = 0.;
  if (invert) {
    a[0] = 0.;
    a[n_delta] = std::numeric_limits<double>::infinity();
  }
  for (auto i = 1; i < n_delta; ++i) {
    if (delta[i - 1] != 0.) {
      a[i] = std::fabs(delta[i] / delta[i - 1]);
    } else {
      a[i] = std::numeric_limits<double>::infinity();
    }
  }

  /// patch conflicts
  int i_low = n_delta;
  int i_upp = 0;
  while (true) {
    bool hit = false;
    for (auto i = 1; i < n_delta; ++i) {
      if ((!invert && a[i] < a[i + 1]) || (invert && a[i] > a[i + 1])) {
        hit = true;
        if ((i == i_low - 1) || (i == i_upp)) {
          /// we're extending a previous range
          i_low = std::min(i_low, i);
          i_upp = std::max(i_upp, i + 1);
        } else {
          /// a start of a new range
          i_low = i;
          i_upp = i + 1;
        }
        /// merge
        double val = 0.;
        if (delta[i_low - 1] != 0.) {
          val = std::fabs(delta[i_upp] / delta[i_low - 1]);
        } else {
          val = std::numeric_limits<double>::infinity();
        }
        val = std::pow(val, 1. / double(i_upp - i_low + 1));
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


std::vector<double> max_a_list(std::vector<double> delta) {
  /// drop all negative values
  auto is_negative = [](const double& d_i) { return d_i < 0.; };
  std::replace_if(delta.begin(), delta.end(), is_negative, 0.);
  return a_list(delta);
}


std::vector<double> min_a_list(std::vector<double> delta) {
  /// if there's at least one negative delta
  ///   -> run max on the flipped numbers
  auto is_negative = [](const double& d_i) { return d_i < 0.; };
  if (std::any_of(delta.begin(), delta.end(), is_negative)) {
    std::transform(delta.begin(), delta.end(), delta.begin(),
                   [](const double& d_i) { return d_i > 0. ? 0. : -d_i; });
    return a_list(delta);
  }
  /// all values are positive from here on
  ///   -> run inverted
  return a_list(delta, true);
}

std::map<std::pair<int, int>, std::pair<double, double>> min_max_partitions(
    const std::vector<double>& delta, const double& max_val) {
  std::map<std::pair<int, int>, std::pair<double, double>> result;

  std::vector<double> min_a = min_a_list(delta);
  std::vector<double> max_a = max_a_list(delta);

  // std::cerr << "min_max_partitions\n";
  // std::cerr << "min_a: [";
  // for (const auto & a : min_a) std::cerr << " " << a;
  // std::cerr << " ]\n";
  // std::cerr << "max_a: [";
  // for (const auto & a : max_a) std::cerr << " " << a;
  // std::cerr << " ]\n";

  const int m = delta.size();
  std::pair<int, int> key;

  if (min_a.back() == 0.) {
    /// min_a & max_a have same ordering
    // std::cerr << "min_a & max_a have same ordering\n";
    double a_low = 0.;
    int k_min = m;
    int k_max = m;
    while ((k_min > 0) && (k_max > 0)) {
      // std::cerr << " [" << a_low << ",";
      /// skip empty intervals
      while (k_min > 1 && min_a[k_min] < max_val &&
             is_approx(min_a[k_min], min_a[k_min - 1]))
        --k_min;
      while (k_max > 1 && max_a[k_max] < max_val &&
             is_approx(max_a[k_max], max_a[k_max - 1]))
        --k_max;
      key = std::pair<int, int>(k_min - 1, k_max - 1);
      if (min_a[k_min - 1] < max_a[k_max - 1]) {
        --k_min;
        if (min_a[k_min] >= max_val) {
          result.emplace(key, std::pair<double, double>(a_low, max_val));
          break;
        }
        result.emplace(key, std::pair<double, double>(a_low, min_a[k_min]));
        a_low = min_a[k_min];
      } else {
        --k_max;
        if (max_a[k_max] >= max_val) {
          result.emplace(key, std::pair<double, double>(a_low, max_val));
          break;
        }
        result.emplace(key, std::pair<double, double>(a_low, max_a[k_max]));
        a_low = max_a[k_max];
      }
      // std::cerr << a_low << "] " << key.first << " & " << key.second << "\n";
    }

  } else {
    /// min_a has reverted ordering w.r.t. a_max
    // std::cerr << "min_a has reverted ordering w.r.t. a_max\n";
    double a_low = 0.;
    int k_min = 0;
    int k_max = m;
    while ((k_min < m) && (k_max > 0)) {
      // std::cerr << " [" << a_low << ",";
      /// skip empty intervals
      while (k_min < m && min_a[k_min] < max_val &&
             is_approx(min_a[k_min + 1], min_a[k_min]))
        ++k_min;
      while (k_max > 1 && max_a[k_max] < max_val &&
             is_approx(max_a[k_max], max_a[k_max - 1]))
        --k_max;
      key = std::pair<int, int>(k_min, k_max - 1);
      if (min_a[k_min + 1] < max_a[k_max - 1]) {
        ++k_min;
        if (min_a[k_min] >= max_val) {
          result.emplace(key, std::pair<double, double>(a_low, max_val));
          break;
        }
        result.emplace(key, std::pair<double, double>(a_low, min_a[k_min]));
        a_low = min_a[k_min];
      } else {
        --k_max;
        if (max_a[k_max] >= max_val) {
          result.emplace(key, std::pair<double, double>(a_low, max_val));
          break;
        }
        result.emplace(key, std::pair<double, double>(a_low, max_a[k_max]));
        a_low = max_a[k_max];
      }
      // std::cerr << a_low << "] " << key.first << " & " << key.second << "\n";
    }
  }
  return result;
}

/// find min/max value for `delta_k / a^k`
std::pair<double, double> get_min_max(const double& a,
                                      const std::vector<double>& delta) {
  double max_val = 1.;  // (i = 0) case
  double min_val = 1.;
  double a_pow_i = 1.;  // accumulate a-powers
  for (auto i = 1; i < delta.size(); ++i) {
    a_pow_i *= a;  // a_pow_i = a^i
    const double dtest = delta.at(i) / a_pow_i;
    if (dtest > max_val) max_val = dtest;
    if (dtest < min_val) min_val = dtest;
  }
  return std::make_pair(min_val, max_val);
}

std::pair<double, double> get_min_max_check(const double& a,
                                            const std::vector<double>& delta) {
  auto parts =
      min_max_partitions(delta, std::numeric_limits<double>::infinity());
  for (const std::pair<std::pair<int, int>, std::pair<double, double>>& p :
       parts) {
    // std::cerr << "[" << p.first.first << "," << p.first.second << "] -> ["
    //           << p.second.first << "," << p.second.second << "]\n";
    if ((p.second.first < a) && (p.second.second >= a)) {
      // std::cerr << "  -->  " << p.second.first << " < " << a << " < "
      //           << p.second.second << std::endl;
      return std::make_pair(
          delta.at(p.first.first) / std::pow(a, p.first.first),
          delta.at(p.first.second) / std::pow(a, p.first.second));
    }
  }

  return std::make_pair(0., 0.);

  // /// find the [a_k, a_{k+1}] interval of the max piece
  // std::vector<double> a_max = miho::max_a_list(delta);
  // int k_max = -1;
  // for (auto i = 1; i < a_max.size(); ++i) {
  //   if ((a_max.at(i - 1) > a) && (a_max.at(i) <= a)) {
  //     k_max = i - 1;
  //     break;
  //   }
  // }
  // double max_val = delta.at(k_max) / std::pow(a, k_max);
  // /// find the [a_k, a_{k+1}] interval of the min piece
  // std::vector<double> a_min = miho::min_a_list(delta);
  // int k_min = -1;
  // for (auto i = 1; i < a_min.size(); ++i) {
  //   const double a_k_low = std::min(a_min.at(i - 1), a_min.at(i));
  //   const double a_k_upp = std::max(a_min.at(i - 1), a_min.at(i));
  //   if ((a_k_low < a) && (a_k_upp >= a)) {
  //     k_min = i - 1;
  //     break;
  //   }
  // }
  // double min_val = delta.at(k_min) / std::pow(a, k_min);
  // return std::make_pair(min_val, max_val);

}

void set_rand_seed(int seed_number) { mt.seed(seed_number); }
double rand() { return dist(mt); }

}  // namespace miho
