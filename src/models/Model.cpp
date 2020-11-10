#include "Model.h"

#include <fmt/format.h>

#include <cmath>
#include <iostream>
#include <list>

namespace miho {

const std::function<double(double)> Model::f_one = [](double x) { return 1.; };
const std::function<double(double)> Model::f_wgt = [](double x) {
  return 1. + 1e1 * std::fabs(x);
};

bool Model::node_exists(const Node& n) {
  for (const Node& inode : _nodes) {
    if (is_approx(inode.x, n.x) && is_approx(inode.y, n.y) &&
        is_approx(inode.sig, n.sig) && is_approx(inode.jac, n.jac))
      return true;
  }
  return false;
}

void Model::print_nodes(const std::string& prefix) {
  for (const Node& n : _nodes) {
    // fmt::print(prefix + "{:<+20.14g} {:<+20.14g}\n", n.sig, n.y);

    fmt::print(prefix + "{:<+20.14g} {:<+20.14g} [{:<+20.14g}, {:<+20.14g}]\n",
               n.sig, n.y, n.x, n.jac);

    // std::cout << n.x << " \t " << n.y << std::endl;
  }
}

// very naive importance sampling to integrate over the pdf.
void Model::adapt_integration(std::function<double(double)> func,
                              std::function<double(double)> rewgt) {
  // std::cout << "\n#integrate at order " << _n_orders << "\n";

  // shift by eps for numerical stability
  double sig_ctr = sigma() * (1. + std::numeric_limits<float>::epsilon());

  auto f_sig = [=](const double& x) { return sig_ctr + std::atanh(x); };
  auto f_jac = [=](const double& x) { return 1. / (1. - x) / (1. + x); };

  /// set up central node
  Node node_i;
  node_i.x = 0;
  node_i.sig = f_sig(node_i.x);
  node_i.jac = f_jac(node_i.x);
  node_i.y = pdf(node_i.sig) * func(node_i.sig);

  if (node_exists(node_i)) {
    /// we continue with the existing list
  } else {
    // std::cout << "#adapt_integration: clearing nodes cache!\n";
    _nodes.clear();
    /// central
    _nodes.push_back(node_i);

    /// add intermediate points (rewgt stability)
    node_i.x = -0.9;
    node_i.sig = f_sig(node_i.x);
    node_i.jac = f_jac(node_i.x);
    node_i.y = pdf(node_i.sig) * func(node_i.sig);
    _nodes.push_front(node_i);
    /// add intermediate points (rewgt stability)
    node_i.x = +0.9;
    node_i.sig = f_sig(node_i.x);
    node_i.jac = f_jac(node_i.x);
    node_i.y = pdf(node_i.sig) * func(node_i.sig);
    _nodes.push_back(node_i);

    /// force ZERO @ -infty
    node_i.x = -1;
    node_i.sig = -std::numeric_limits<double>::infinity();
    node_i.jac = 0.;  // +std::numeric_limits<double>::infinity();
    node_i.y = 0.;
    _nodes.push_front(node_i);
    /// force ZERO @ +infty
    node_i.x = +1;
    node_i.sig = +std::numeric_limits<double>::infinity();
    node_i.jac = 0.;  // +std::numeric_limits<double>::infinity();
    node_i.y = 0.;
    _nodes.push_back(node_i);
  }

  /// we want a minimum # of successive node insertions with desired acc
  int nsucc = _min_nodes;
  size_t count = 0;
  while (_nodes.size() < _max_nodes) {
    count++;
    // print_nodes();
    double jump = -1.;
    double result = 0.;
    double error = 0.;
    auto node_pos = _nodes.begin();
    for (auto it = std::next(_nodes.begin()); it != _nodes.end(); ++it) {
      /// accumulators
      double dres = (it->x - std::prev(it)->x) *
                    (it->y * it->jac + std::prev(it)->y * std::prev(it)->jac) /
                    2.;
      result += dres;
      double df_prev =
          (it->y * it->jac - std::prev(it)->y * std::prev(it)->jac) /
          (it->x - std::prev(it)->x);
      double df_next =
          (std::next(it)->y * std::next(it)->jac - it->y * it->jac) /
          (std::next(it)->x - it->x);
      double derr = std::fabs(df_prev - df_next) *
                    pow(it->x - std::prev(it)->x, 2);  // / 12.;
      error += derr;

      /// find the next subdivision
      double test = rewgt((it->x + std::prev(it)->x) / 2.);
      if (count % 3 == 0) {
        test *= dres;  // importance sampling
      } else {
        test *= derr;  // stratified sampling
      }
      // std::cout << "test @ " << it->x << ": " << test << std::endl;
      if (test > jump) {
        jump = test;
        node_pos = it;
        // std::cout << "jump!" << std::endl;
      }
    }

    /// done?
    if (std::fabs(error / result) <= _target_accuracy) {
      nsucc--;
      if (nsucc < 0) break;
      // std::cerr << "reached target accuracy[" << nsucc << "]: " << result
      //           << " +/- " << error << " [" << error / result << "/"
      //           << _target_accuracy << "] in " << _nodes.size() << "
      //           steps\n";
    } else {
      nsucc = _min_nodes;
    }

    /// insert new element
    node_i.x = (node_pos->x + std::prev(node_pos)->x) / 2.;
    node_i.sig = f_sig(node_i.x);
    node_i.jac = f_jac(node_i.x);
    node_i.y = pdf(node_i.sig) * func(node_i.sig);
    node_pos = _nodes.insert(node_pos, node_i);
    // fmt::print(
    //     "# new element: {:<+20.14g} {:<+20.14g} [{:<+20.14g},
    //     {:<+20.14g}]\n", node_i.sig, node_i.y, node_i.x, node_i.jac);
    // std::cin.ignore();
  }
}

double Model::integrate(std::function<double(double)> func) {
  adapt_integration(func);

  // std::cout << "\n\n#integrate: adapted nodes!\n";
  // print_nodes("#DEBUG: ");

  double result = 0.;
  for (auto it = _nodes.begin(); it != _nodes.end(); ++it) {
    // fmt::print("{:<+20.16g} {:<+20.16g}  # = (x,y)\n", it->x, it->y);
    if (it == _nodes.begin()) continue;
    result += (it->x - std::prev(it)->x) *
              (it->y * it->jac + std::prev(it)->y * std::prev(it)->jac) / 2.;
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
    percentile += (it->x - std::prev(it)->x) *
                  (it->y * it->jac + std::prev(it)->y * std::prev(it)->jac) /
                  2.;
    // std::cout << "# low " << it->x << ": " << percentile << std::endl;
    if (percentile >= (1.0 - p) / 2.) {
      lower = it->sig;
      break;
    }
  }
  percentile = 0.;
  for (auto it = std::next(_nodes.rbegin()); it != _nodes.rend(); ++it) {
    percentile += (std::prev(it)->x - it->x) *
                  (std::prev(it)->y * std::prev(it)->jac + it->y * it->jac) /
                  2.;
    // std::cout << "# upp " << it->x << ": " << percentile << std::endl;
    if (percentile >= (1.0 - p) / 2.) {
      upper = it->sig;
      break;
    }
  }
  return std::pair<double, double>{lower, upper};
}

}  // namespace miho
