#include "Model.h"

#include <fmt/format.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <list>

// #include "ABCModel.h"
// #include "GeometricModel.h"

namespace miho {

// std::unique_ptr<Model> Model::create(ModelType model_type) {
//   switch (model_type) {
//     case ModelType::geo:
//       return std::make_unique<GeometricModel>();
//       break;
//     case ModelType::abc:
//       return std::make_unique<ABCModel>();
//       break;
//     default:
//       throw "Model::create: model not registered";
//   }
// }

const std::function<double(double)> Model::f_one = [](double x) { return 1.; };
const std::function<double(double)> Model::f_wgt = [](double x) {
  return 1. + 1e2 * std::fabs(x);
};

bool Model::node_exists(const Node& n) const {
  for (const Node& inode : _nodes) {
    if (is_approx(inode.x, n.x) && is_approx(inode.sig, n.sig) &&
        is_approx(inode.jac, n.jac) && is_approx(inode.pdf, n.pdf) &&
        is_approx(inode.func, n.func))
      return true;
  }
  return false;
}

void Model::print_nodes(const std::string& prefix) {
  for (const Node& n : _nodes) {
    fmt::print(prefix + "{:<+20.14g} {:<+20.14g}\n", n.sig, n.pdf);

    // fmt::print(prefix +
    //         "{:<+20.14g} {:<+20.14g} {:<+20.14g} [{:<+20.14g}, {:<+20.14g}]\n",
    //     n.sig, n.pdf, n.func, n.x, n.jac);

    // std::cout << n.x << " \t " << n.pdf << std::endl;
  }
}

// very naive importance sampling to integrate over the pdf.
void Model::adapt_integration(std::function<double(double)> func,
                              std::function<double(double)> rewgt) {
  // std::cout << "\n#integrate at order " << _n_orders << "\n";

  // shift by eps for numerical stability
  const double sig_ctr = sigma() * (1. + std::numeric_limits<float>::epsilon());

  // auto f_sig = [&sig_ctr](const double& x) {
  //   return sig_ctr * (1 + std::atanh(x));
  // };
  // auto f_jac = [&sig_ctr](const double& x) {
  //   return sig_ctr / (1. - x) / (1. + x);
  // };

  auto f_sig = [&sig_ctr](const double& x) {
    const volatile double oMx_sq = (1 - x) * (1 + x);
    const volatile double sgn_x = (sig_ctr>=0) ? x : -x;;
    return sig_ctr * (1 + sgn_x / oMx_sq);
  };
  auto f_jac = [&sig_ctr](const double& x) {
    const volatile double oMx_sq = (1 - x) * (1 + x);
    return sig_ctr * (1. + x * x) / (oMx_sq * oMx_sq);
  };

  /// set up central node
  Node node_i;

  // /// this is for optimisation for each integrand type...
  // node_i.x = 0;
  // node_i.sig = f_sig(node_i.x);
  // node_i.jac = f_jac(node_i.x);
  // node_i.pdf = pdf(node_i.sig);
  // node_i.func = func(node_i.sig);
  // if (node_exists(node_i)) {
  //   /// we continue with the existing list
  // } else {
  //   std::cerr << "#adapt_integration: clearing nodes cache!\n";
  //   clear();
  // }

  if (_nodes.empty()) {
    // std::cerr << "# adapt integration: start fresh..." << std::endl;
    /// central
    node_i.x = 0;
    node_i.sig = f_sig(node_i.x);
    node_i.jac = f_jac(node_i.x);
    node_i.pdf = pdf(node_i.sig);
    node_i.func = func(node_i.sig);
    _nodes.push_back(node_i);

    /// add intermediate points (rewgt stability)
    node_i.x = -0.9;
    node_i.sig = f_sig(node_i.x);
    node_i.jac = f_jac(node_i.x);
    node_i.pdf = pdf(node_i.sig);
    node_i.func = func(node_i.sig);
    _nodes.push_front(node_i);
    /// add intermediate points (rewgt stability)
    node_i.x = +0.9;
    node_i.sig = f_sig(node_i.x);
    node_i.jac = f_jac(node_i.x);
    node_i.pdf = pdf(node_i.sig);
    node_i.func = func(node_i.sig);
    _nodes.push_back(node_i);

    /// force ZERO @ -infty
    node_i.x = -1;
    node_i.sig = -std::numeric_limits<double>::infinity();
    node_i.jac = 0.;  // +std::numeric_limits<double>::infinity();
    node_i.pdf = 0.;
    node_i.func = 0.;
    _nodes.push_front(node_i);
    /// force ZERO @ +infty
    node_i.x = +1;
    node_i.sig = +std::numeric_limits<double>::infinity();
    node_i.jac = 0.;  // +std::numeric_limits<double>::infinity();
    node_i.pdf = 0.;
    node_i.func = 0.;
    _nodes.push_back(node_i);
  } else {
    /// update the function values
    // std::cerr << "# adapt integration: continue..." << std::endl;
    for (Node& node : _nodes) {
      node.func = func(node.sig);
    }
    // reset first & last entry (by definition +/- infinity)
    _nodes.front().func = 0.;
    _nodes.back().func = 0.;
  }

  /// we want a minimum # of successive node insertions with desired acc
  int nsucc = _min_nodes;
  size_t count = 0;
  while (_nodes.size() < _max_nodes) {
    count++;
    // print_nodes();
    double jump = -1.;
    double norm = 0.;
    double result = 0.;
    double error = 0.;
    double error2 = 0.;
    auto node_pos = _nodes.begin();
    for (auto it = std::next(_nodes.begin()); it != _nodes.end(); ++it) {
      /// accumulators
      auto it_prev = std::prev(it);
      auto it_next = std::next(it);
      double dx_prev = it->x - it_prev->x;
      double dx_next = it_next->x - it->x;
      norm += dx_prev * (it->pdf * it->jac + it_prev->pdf * it_prev->jac) / 2.;
      double dres = dx_prev *
                    (it->pdf * it->func * it->jac +
                     it_prev->pdf * it_prev->pdf * it_prev->jac) /
                    2.;
      result += dres;
      double df_prev = (it->pdf * it->func * it->jac -
                        it_prev->pdf * it_prev->func * it_prev->jac) /
                       dx_prev;
      double df_next = (it_next->pdf * it_next->func * it_next->jac -
                        it->pdf * it->func * it->jac) /
                       dx_next;
      double derr = std::fabs(df_prev - df_next) * pow(dx_prev, 2);  // / 12.;
      error += derr;
      double ddf = 2. * (df_next - df_prev) / (it_next->x - it_prev->x);
      double err2 = std::fabs(ddf) *
                    std::pow((it_next->x - it_prev->x) / 2., 2) * 2. / 12.;
      if (err2 > error2) error2 = err2;

      /// find the next subdivision
      double test = rewgt((it->x + it_prev->x) / 2.);
      if (count % 3 == 0) {
        // std::cout << " - importance - ";
        test *= dres;  // importance sampling
      } else {
        // std::cout << " - stratified - ";
        test *= derr;  // stratified sampling
      }
      // std::cout << "test @ " << it->x << ": " << it->sig << " ->" << test << std::endl;
      if (test > jump) {
        jump = test;
        node_pos = it;
        // std::cout << "jump!" << std::endl;
      }
    }

    /// done?
    if (std::fabs(error / result) <= _target_accuracy &&
        std::fabs(error2 / result) <= _target_accuracy &&
        std::fabs(norm - 1.) <= _target_accuracy) {
      nsucc--;
      if (nsucc < 0) break;
      // std::cerr << "reached target accuracy[" << nsucc << "]: " << result
      //           << " +/- " << error << " [" << error / result << "/"
      //           << _target_accuracy << " | " << norm << "] in " <<
      //           _nodes.size()
      //           << "steps\n";
    } else {
      nsucc = _min_nodes;
    }

    /// insert new element
    node_i.x = (node_pos->x + std::prev(node_pos)->x) / 2.;
    if (is_approx(node_i.x, -1.) || is_approx(node_i.x, -1.)) break;
    node_i.sig = f_sig(node_i.x);
    node_i.jac = f_jac(node_i.x);
    node_i.pdf = pdf(node_i.sig);
    node_i.func = func(node_i.sig);
    node_pos = _nodes.insert(node_pos, node_i);
    // fmt::print(
    //     "# new element: {:<+20.14g} {:<+20.14g} [{:<+20.14g},{:<+20.14g}]\n",
    //     node_i.sig, node_i.func, node_i.x, node_i.jac);
    // std::cin.ignore();
  }
}

double Model::integrate(std::function<double(double)> func) {
  adapt_integration(func);

  // std::cout << "\n\n#integrate: adapted nodes!\n";
  // print_nodes("#DEBUG: ");

  double result = 0.;
  for (auto it_curr = _nodes.begin(); it_curr != _nodes.end(); ++it_curr) {
    const auto it_prev = std::prev(it_curr);
    // fmt::print("{:<+20.16g} {:<+20.16g}  # = (x,y)\n", it_curr->x,
    // it_curr->y);
    if (!std::isfinite(it_curr->x) || !std::isfinite(it_prev->x) ||
        it_curr == _nodes.begin())
      continue;
    result += (it_curr->x - it_prev->x) *
              (it_curr->pdf * it_curr->func * it_curr->jac +
               it_prev->pdf * it_prev->func * it_prev->jac) /
              2.;
  }
  // std::cout << "# result = " << result << std::endl;

  return result;
}

std::pair<double, double> Model::degree_of_belief_interval(const double& p) {
  adapt_integration(f_one, f_wgt);

  double lower = 0.;
  double upper = 0.;

  double percentile = 0.;
  for (auto it = std::next(_nodes.begin()); it != _nodes.end(); ++it) {
    percentile +=
        (it->x - std::prev(it)->x) *
        (it->pdf * it->jac + std::prev(it)->pdf * std::prev(it)->jac) / 2.;
    // std::cout << "# low " << it->x << ": " << percentile << std::endl;
    if (percentile >= (1.0 - p) / 2.) {
      lower = it->sig;
      break;
    }
  }
  percentile = 0.;
  for (auto it = std::next(_nodes.rbegin()); it != _nodes.rend(); ++it) {
    percentile +=
        (std::prev(it)->x - it->x) *
        (std::prev(it)->pdf * std::prev(it)->jac + it->pdf * it->jac) / 2.;
    // std::cout << "# upp " << it->x << ": " << percentile << std::endl;
    if (percentile >= (1.0 - p) / 2.) {
      upper = it->sig;
      break;
    }
  }
  return std::pair<double, double>{lower, upper};
}

}  // namespace miho
