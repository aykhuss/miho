#pragma once

#include <iostream>
#include <vector>

#include "Model.h"
#include "ModelPrototype.h"

namespace miho {

class ABCNumericModel : public ModelPrototype {
 public:
  ABCNumericModel()
      : _sigma(),
        _delta(),
        _q_pdf_den(false),
        _omega(1),
        _xi(5.0),
        _epsilon(0.1),
        _eta(0.2),
        _nint_rel_err(0.01) {}
  ABCNumericModel(const std::vector<double>& sigma)
      : _sigma(sigma),
        _delta(),
        _q_pdf_den(false),
        _omega(1),
        _xi(5.0),
        _epsilon(0.1),
        _eta(0.2),
        _nint_rel_err(0.01) {
    init();
  }

  inline void set_sigma(const std::vector<double>& sigma) override {
    _sigma = sigma;
    _q_pdf_den = false;
    init();
  }

  void init() {
    // std::cout << "ABCNumericModel - init - " << _omega << ", " << _epsilon
    //           << std::endl;
    // update everything after changing _sigma
    clear();  // clear cached nodes
    _n_orders = _sigma.size();
    _delta.clear();
    _delta.reserve(_n_orders);
    // Eq.(3.4)
    _delta.push_back(1.);  // normalised w.r.t. LO <-> _sigma[0]
    // std::cerr << "ABCNumericModel::init [" << _delta.back();
    for (auto i = 1; i < _n_orders; ++i) {
      _delta.push_back((_sigma[i] - _sigma[i - 1]) / _sigma[0]);
      // std::cerr << ", " << _delta.back();
    }
    // std::cerr << "] " << std::endl;
  }

  /// Model interface
  double sigma(int order) const override { return _sigma.at(order); };
  double pdf(const double& val) const override;
  /// ModelPrototype interface
  std::unique_ptr<ModelPrototype> clone() const override {
    return std::make_unique<ABCNumericModel>(*this);
  }

  // setters
  inline void set_epsilon(const double& epsilon) {
    _epsilon = epsilon;
    _q_pdf_den = false;
  }
  inline void set_eta(const double& eta) {
    _eta = eta;
    _q_pdf_den = false;
  }
  inline void set_xi(const double& xi) {
    _xi = xi;
    _q_pdf_den = false;
  }
  inline void set_omega(int omega) {
    _omega = omega;
    _q_pdf_den = false;
  }
  inline void set_nint_rel_err(const double& acc) { _nint_rel_err = acc; }

  // additional public member functions
  inline double delta_next(const double& sigma_next) const {
    return (sigma_next - _sigma.back()) / _sigma.front();
  }
  double pdf_delta___delta_mu(const double& delta_next) const;
  double pdf_delta__mu(const std::vector<double>& delta) const;
  double pdf_delta__mu__2D(const std::vector<double>& delta) const;
  inline double pdf_delta__mu() const { return pdf_delta__mu(_delta); }
  double pdf_delta__mu(const double& delta_next) const;

 protected:
  std::vector<double> _sigma;
  std::vector<double> _delta;
  // cache denominator
  mutable bool _q_pdf_den;
  mutable double _pdf_den;
  // parameters of the model
  int _omega;
  double _xi;
  double _epsilon;
  double _eta;
  // parameter for the cubature integration
  double _nint_rel_err;  // 1%
};

}  // namespace miho
