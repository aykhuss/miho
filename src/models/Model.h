#pragma once

#include <cmath>
#include <iostream>
#include <limits>
#include <list>
#include <memory>
#include <stdexcept>

#include "Util.h"

namespace miho {

class Model {
 public:
  Model()
      : _n_orders{0},
        _target_accuracy{0.001},
        _max_nodes{10000},
        _min_nodes{10},
        _nodes{} {}
  virtual ~Model() noexcept = default;
  virtual double sigma(int order) const = 0;
  inline double sigma() const { return sigma(_n_orders - 1); };
  inline size_t n_orders() const { return _n_orders; };
  virtual double pdf(const double& val) const = 0;

  // /// The factory design pattern
  // static std::unique_ptr<Model> create(ModelType model_type);

  void adapt_integration(std::function<double(double)> func,
                         std::function<double(double)> rewgt = f_one);
  double integrate(std::function<double(double)> func);

  inline double norm() { return integrate(f_one); }
  inline double moment(int p) {
    return integrate([p](double x) { return std::pow(x, p); });
  }
  inline double mean() {
    // double old_target_accuracy = _target_accuracy;
    // set_accuracy(_target_accuracy / 4.);
    // double result = integrate([](double x) { return x; });
    // set_accuracy(old_target_accuracy);
    // return result;
    return integrate([](double x) { return x; });
    // return moment(1);
  }
  inline double variance() {
    // this is not very efficient, as there are large cancellations!
    // return moment(2) - std::pow(mean(), 2);
    // much better convergence this way:
    double mean_val = mean();
    return integrate(
        [mean_val](double x) { return std::pow(x - mean_val, 2); });
    // double old_target_accuracy = _target_accuracy;
    // set_accuracy(_target_accuracy / 4.);
    // double result =
    //     integrate([=](double x) { return std::pow(x - mean_val, 2); });
    // set_accuracy(old_target_accuracy);
    // return result;
  }
  inline double stdev() { return std::sqrt(variance()); }
  std::pair<double, double> degree_of_belief_interval(const double& p = 0.68);
  double median() {
    std::pair<double, double> dob = degree_of_belief_interval(0.);
    return (dob.first + dob.second) / 2.;
  }
  void print_nodes(const std::string& prefix = "");
  inline void print_pdf() {
    adapt_integration(f_one, f_wgt);
    print_nodes();
  }
  inline void set_max_nodes(size_t nmax) { _max_nodes = nmax; }
  inline void set_accuracy(const double& acc) { _target_accuracy = acc; }
  inline void clear() { _nodes.clear(); }

 protected:
  size_t _n_orders;
  double _target_accuracy;
  size_t _max_nodes;
  size_t _min_nodes;

 private:
  static const std::function<double(double)> f_one;
  static const std::function<double(double)> f_wgt;
  struct Node {
    double x;
    double sig;
    double jac;
    double pdf;
    double func;
  };
  bool node_exists(const Node& n) const;
  inline double node_y(const Node& n) const { return n.jac * n.pdf * n.func; }

  std::list<Node> _nodes;
};

}  // namespace miho
