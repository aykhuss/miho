#pragma once

#include <iostream>
#include <vector>

#include "Model.h"

namespace miho {

class ABCNumericModel : public Model {
 public:
  ABCNumericModel(const std::vector<double>& sigma)
      : _sigma(sigma), _delta(), _epsilon(0.1), _xi(1.0), _omega(1) {
    _n_orders = sigma.size();
    _delta.reserve(_n_orders);
    // Eq.(3.4)
    _delta.push_back(1.);  // normalised w.r.t. LO <-> sigma[0]
    for (auto i = 1; i < _n_orders; ++i) {
      _delta.push_back((sigma[i] - sigma[i - 1]) / sigma[0]);
    }
    // std::cout << "# ABCNumericModel: " << _n_orders << " \n";
    // for (auto i = 0; i < sigma.size(); ++i) {
    //   std::cout << "# > " << sigma[i] << ", " << _sigma[i] << ":\t" <<
    //   _delta[i]
    //             << std::endl;
    // }
    // cache the denominator (const for fixed delta's)
    _pdf_den = pdf_delta__mu(_delta);
  }

  double sigma(int order) const { return _sigma.at(order); };
  double pdf(const double& val) const;

  // setters
  inline void set_epsilon(const double& epsilon) {
    _epsilon = epsilon;
    _pdf_den = pdf_delta__mu(_delta);
  }
  inline void set_xi(const double& xi) {
    _xi = xi;
    _pdf_den = pdf_delta__mu(_delta);
  }
  inline void set_omega(int omega) {
    _omega = omega;
    _pdf_den = pdf_delta__mu(_delta);
  }

  // additional public member functions
  inline double delta_next(const double& sigma_next) const {
    return (sigma_next - _sigma.back()) / _sigma.front();
  }
  double pdf_delta___delta_mu(const double& delta_next) const;
  double pdf_delta__mu(const std::vector<double>& delta) const;
  inline double pdf_delta__mu() const { return pdf_delta__mu(_delta); }
  double pdf_delta__mu(const double& delta_next) const;

 protected:
  std::vector<double> _sigma;
  std::vector<double> _delta;
  double _pdf_den;
  // parameters of the model
  int _omega;
  double _xi;
  double _epsilon;
};

}  // namespace miho
