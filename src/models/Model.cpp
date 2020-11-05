#include "Model.h"

#include <fmt/format.h>

#include <cmath>
#include <iostream>
#include <list>

namespace miho {

bool Model::node_exists(const Node& n) {
  for (const Node& inode : _nodes) {
    if (is_approx(inode.x, n.x) && is_approx(inode.y, n.y)) return true;
  }
  return false;
}

// very naive importance sampling to integrate over the pdf.
void Model::adapt_integration(std::function<double(double)> func) {
  // std::cout << "\n#integrate at order " << _n_orders << "\n";
  // find three nodes to seed the integration
  double min_sig = -1.;
  double max_sig = -1.;
  for (auto i_ord = 0; i_ord < _n_orders; ++i_ord) {
    double sig = sigma(i_ord);
    if ((min_sig < 0.) || (sig < min_sig)) min_sig = sig;
    if ((max_sig < 0.) || (sig > max_sig)) max_sig = sig;
  }
  double integrate_center =
      sigma() * (1. + std::numeric_limits<float>::epsilon());
  double integrate_delta = (max_sig - min_sig) / 2.;
  if (integrate_delta < 0.1 * integrate_center)
    integrate_delta = 0.1 * integrate_center;

  // std::cout << "#range: " << integrate_center << " +- " << integrate_delta
  //           << std::endl;

  Node nctr, nlow, nupp;
  // highest-order central prediction
  nctr.x = integrate_center;
  nctr.y = pdf(nctr.x) * func(nctr.x);
  // fmt::print("# > center: ({},{})\n", nctr.x, nctr.y);
  // lower variation (ordered list)
  nlow.x = integrate_center - integrate_delta;
  nlow.y = pdf(nlow.x) * func(nlow.x);
  // fmt::print("# > lower:  ({},{})\n", nlow.x, nlow.y);
  // upper variation (ordered list)
  nupp.x = integrate_center + integrate_delta;
  nupp.y = pdf(nupp.x) * func(nupp.x);
  // fmt::print("# > upper:  ({},{})\n", nupp.x, nupp.y);

  if (node_exists(nctr) && node_exists(nlow) && node_exists(nupp)) {
    // we continue with the existing list
  } else {
    // std::cout << "# clearing nodes cache!\n";
    _nodes.clear();
    _nodes.emplace_front(nctr);
    _nodes.emplace_front(nlow);
    _nodes.emplace_back(nupp);
  }

  while (_nodes.size() < _max_nodes) {
    double jump = 0.;
    double result = 0.;
    double error = 0.;
    double last_derivative = 0.;
    double curr_derivative = 0.;
    auto node_pos = _nodes.begin();
    for (auto it = _nodes.begin(); it != _nodes.end(); ++it) {
      // compute some accumulator
      if (it == _nodes.begin()) continue;
      double dres =
          (it->x - std::prev(it)->x) * (it->y + std::prev(it)->y) / 2.;
      result += dres;
      dres = std::fabs(dres);
      if (dres > jump) {
        jump = dres;
        node_pos = it;
        // std::cout << "dres: " << dres << std::endl;
      }
      // compute the error
      curr_derivative = (it->y - std::prev(it)->y) / (it->x - std::prev(it)->x);
      error += std::fabs(curr_derivative - last_derivative) *
               pow(it->x - std::prev(it)->x, 2) / 12.;
      last_derivative = curr_derivative;
    }

    // place where the termination condition (target accuracy) would go in

    Node n_new;
    n_new.x = (node_pos->x + std::prev(node_pos)->x) / 2.;
    // fmt::print("# > normal jump? {} by {}\n", n_new.x, jump);

    // possibly extend range?
    bool qextended = false;
    double x0, dres0;
    Node n1, n2;
    const double weight = 1e3;
    // front:
    n1 = *_nodes.begin();
    n2 = *std::next(_nodes.begin());
    x0 = n1.x - (n2.x - n1.x) * n1.y / (n2.y - n1.y) * weight;
    if (x0 >= n1.x) x0 = n1.x - weight * (n2.x - n1.x);
    dres0 = std::fabs((n1.x - x0) * n1.y / 2.);
    if (dres0 > jump) {
      qextended = true;
      jump = dres0;
      node_pos = _nodes.begin();
      n_new.x = x0;
      // fmt::print("# > front jump? {} by {}\n", n_new.x, jump);
    }
    // back:
    n1 = *_nodes.rbegin();
    n2 = *std::next(_nodes.rbegin());
    x0 = n1.x - (n2.x - n1.x) * n1.y / (n2.y - n1.y) * weight;
    if (x0 <= n1.x) x0 = n1.x - weight * (n2.x - n1.x);
    dres0 = std::fabs((n1.x - x0) * n1.y / 2.);
    if (dres0 > jump) {
      qextended = true;
      jump = dres0;
      node_pos = _nodes.end();
      n_new.x = x0;
      // fmt::print("# > back jump? {} by {}\n", n_new.x, jump);
    }

    // insert new element
    // fmt::print("#  >> new x = {}\n", n_new.x);
    n_new.y = pdf(n_new.x) * func(n_new.x);
    // fmt::print("#  >> new y = {}\n", n_new.y);
    node_pos = _nodes.insert(node_pos, n_new);
    // fmt::print("# new element: ({},{}) \n", n_new.x, n_new.y);

    // check for accuracy termination condition
    if (!qextended) {
      // // guaranteed to be somewhere in the "middle"
      // double region_old = (std::next(node_pos)->x - std::prev(node_pos)->x) *
      //                     (std::next(node_pos)->y + std::prev(node_pos)->y) /
      //                     2.;
      // double region_new = 0.;
      // region_new += (std::next(node_pos)->x - node_pos->x) *
      //               (std::next(node_pos)->y + node_pos->y) / 2.;
      // region_new += (node_pos->x - std::prev(node_pos)->x) *
      //               (node_pos->y + std::prev(node_pos)->y) / 2.;
      // std::cout << "region_old: " << region_old << std::endl;
      // std::cout << "region_new: " << region_new << std::endl;
      // double check = 1e4 * std::fabs(region_old - region_new) * _nodes.size();
      // double check = 1.5 * std::fabs(jump); // empirical prefactor
      // double check = 1.5 * std::fabs(jump) +
      //                std::fabs(region_old - region_new) * _nodes.size();
      double check = error;
      if (std::fabs(check / result) <= _target_accuracy) {
        // std::cout << "reached target accuracy: " << result << " +/- " << check
        //           << " [" << check / result << "/" << _target_accuracy
        //           << "] in " << _nodes.size() << " steps\n";
        break;
      }
    }
  }
}

double Model::integrate(std::function<double(double)> func) {
  adapt_integration(func);

  double result = 0.;
  for (auto it = _nodes.begin(); it != _nodes.end(); ++it) {
    // fmt::print("{:<+20.16g} {:<+20.16g}  # = (x,y)\n", it->x, it->y);
    if (it == _nodes.begin()) continue;
    result += (it->x - std::prev(it)->x) * (it->y + std::prev(it)->y) / 2.;
  }
  // std::cout << "# result = " << result << std::endl;

  return result;
}

std::pair<double, double> Model::degree_of_belief_interval(const double& p) {
  adapt_integration([](double x) { return 1.; });

  double lower = 0.;
  double upper = 0.;

  double percentile = 0.;
  for (auto it = std::next(_nodes.begin()); it != _nodes.end(); ++it) {
    percentile += (it->x - std::prev(it)->x) * (it->y + std::prev(it)->y) / 2.;
    // std::cout << "# low " << it->x << ": " << percentile << std::endl;
    if (percentile >= (1.0 - p) / 2.) {
      lower = it->x;
      break;
    }
  }
  percentile = 0.;
  for (auto it = std::next(_nodes.rbegin()); it != _nodes.rend(); ++it) {
    percentile += (std::prev(it)->x - it->x) * (it->y + std::prev(it)->y) / 2.;
    // std::cout << "# upp " << it->x << ": " << percentile << std::endl;
    if (percentile >= (1.0 - p) / 2.) {
      upper = it->x;
      break;
    }
  }
  return std::pair<double, double>{lower, upper};
}

}  // namespace miho
