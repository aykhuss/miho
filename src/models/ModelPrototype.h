#pragma once

#include <iostream>
#include <memory>
#include <vector>

#include "Model.h"

namespace miho {

enum class ModelType { geo, abc };

class ModelPrototype : public Model {
 public:
  using Model::sigma;  /// avoid name hiding
  virtual ~ModelPrototype() noexcept = default;
  /// The prototype design pattern
  virtual std::unique_ptr<ModelPrototype> clone() const = 0;
  /// external interface
  void set_sigma(const std::vector<double>& sigma) {
    _sigma = sigma;
    /// update everything after changing _sigma
    clear();  // clear caches
    _n_orders = _sigma.size();
    _delta.clear();
    _delta.reserve(_n_orders);
    /// Eq.(3.4)
    _delta.push_back(1.);  // normalised w.r.t. LO <-> _sigma[0]
    for (auto i = 1; i < _n_orders; ++i) {
      _delta.push_back((_sigma[i] - _sigma[i - 1]) / _sigma[0]);
    }
  }
  inline double sigma(int order) const override { return _sigma.at(order); };
  /// Model interface
  inline double pdf(const double& val) const override {
    /// Eq.(3.14) : assumes j == 1
    return pdf_delta___delta_mu(delta_next(val)) / std::fabs(_sigma.front());
  }

  /// additional internal member functions
  inline double delta_next(const double& sigma_next) const {
    return (sigma_next - _sigma.back()) / _sigma.front();
  }
  double pdf_delta___delta_mu(const double& delta_next) const {
    /// Eq.(4.10) : assumes j == 1
    if (!_q_pdf_den) {
      _pdf_den = pdf_delta__mu(_delta);
      _q_pdf_den = true;
    }
    const double pdf_num = pdf_delta__mu(delta_next);
    return pdf_num / _pdf_den;
  }
  virtual double pdf_delta__mu(const std::vector<double>& delta) const = 0;
  inline double pdf_delta__mu() const { return pdf_delta__mu(_delta); }
  double pdf_delta__mu(const double& delta_next) const {
    std::vector<double> delta = _delta;
    delta.push_back(delta_next);
    return pdf_delta__mu(delta);
  }

 protected:
  std::vector<double> _sigma;
  std::vector<double> _delta;
  /// cache denominator
  mutable bool _q_pdf_den = false;
  mutable double _pdf_den;
  virtual void clear() override {
    _q_pdf_den = false;
    Model::clear();
  }

};

}  // namespace miho
