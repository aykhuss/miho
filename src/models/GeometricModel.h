#pragma once

#include <iostream>
#include <vector>

#include "Model.h"
#include "ModelPrototype.h"

namespace miho {

class GeometricModel : public ModelPrototype {
 public:
  GeometricModel()
      : _sigma(), _delta(), _q_pdf_den(false), _omega(1), _epsilon(0.1) {}
  GeometricModel(const std::vector<double>& sigma)
      : _sigma(sigma), _delta(), _q_pdf_den(false), _omega(1), _epsilon(0.1) {
    init();
  }

  inline void set_sigma(const std::vector<double>& sigma) override {
    _sigma = sigma;
    _q_pdf_den = false;
    init();
  }

  void init() {
    // std::cout << "GeometricModel - init - " << _omega << ", " << _epsilon
    //           << std::endl;
    // update everything after changing _sigma
    clear();  // clear cached nodes
    _n_orders = _sigma.size();
    _delta.clear();
    _delta.reserve(_n_orders);
    // Eq.(3.4)
    _delta.push_back(1.);  // normalised w.r.t. LO <-> _sigma[0]
    for (auto i = 1; i < _n_orders; ++i) {
      _delta.push_back((_sigma[i] - _sigma[i - 1]) / _sigma[0]);
    }
  }

  /// Model interface
  double sigma(int order) const override { return _sigma.at(order); };
  double pdf(const double& val) const override;
  /// ModelPrototype interface
  std::unique_ptr<ModelPrototype> clone() const override {
    return std::make_unique<GeometricModel>(*this);
  }

  // setters
  inline void set_epsilon(const double& epsilon) {
    _epsilon = epsilon;
    _q_pdf_den = false;
  }
  inline void set_omega(int omega) {
    _omega = omega;
    _q_pdf_den = false;
  }

  // additional public member functions
  inline double delta_next(const double& sigma_next) const {
    return (sigma_next - _sigma.back()) / _sigma.front();
  }
  double pdf_delta___delta_mu(const double& delta_next) const;
  double pdf_delta__mu(const std::vector<double>& delta) const;
  inline double pdf_delta__mu() const { return pdf_delta__mu(_delta); }
  double pdf_delta__mu(const double& delta_next) const;

  static std::vector<double> a_list(const std::vector<double>& delta);
  static double a_integral(const double& a_lower, const double& a_upper,
                           const double& epsilon, int m, int k, int j);

 protected:
  std::vector<double> _sigma;
  std::vector<double> _delta;
  // cache denominator
  mutable bool _q_pdf_den;
  mutable double _pdf_den;
  // parameters of the model
  int _omega;       // Eq.(4.9) : assumes integer
  double _epsilon;  // Eq.(4.8)
};

}  // namespace miho
