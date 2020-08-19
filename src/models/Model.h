#pragma once

#include <cmath>
#include <iostream>
#include <list>
#include <stdexcept>

namespace miho {

class Model {
 public:
  Model() : _n_orders{0} {}
  virtual double sigma(int order) const = 0;
  inline double sigma() const { return sigma(_n_orders - 1); };
  inline size_t n_orders() const { return _n_orders; };
  virtual double pdf(const double& val) const = 0;
  double integrate(std::function<double(double)> func);
  // template<size_t p>
  // double moment() {
  //   return integrate([](double x) { return std::pow(x,p); });
  // }
  double moment(int p) {
    return integrate([=](double x) { return std::pow(x, p); });
  }
  double mean() {
    return integrate([](double x) { return x; });
    // return moment(1);
  }
  std::pair<double, double> degree_of_belief_interval(const double& p = 0.68);
  double median() {
    std::pair<double, double> dob = degree_of_belief_interval(0.5);
    return (dob.first + dob.second) / 2.;
  }
  void print_nodes() {
    adapt_integration([](double x) { return 1.; });
    for (const Node& n : _nodes) {
      std::cout << n.x << " \t " << n.y << std::endl;
    }
  }

 protected:
  size_t _n_orders;

  static bool is_approx(const double& x, const double& y) {
    double max_val = std::max({1.0, std::fabs(x), std::fabs(y)});
    return std::fabs(x - y) <= std::numeric_limits<double>::epsilon() * max_val;
  }

 private:
  struct Node {
    double x;
    double y;
  };

  bool node_exists(const Node& n);
  void adapt_integration(std::function<double(double)> func);

  std::list<Node> _nodes;
};

}  // namespace miho
