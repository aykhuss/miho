#pragma once

#include <iostream>
#include <vector>

#include "Model.h"

namespace miho {

class ABCModel : public Model {
 public:
  ABCModel(const std::vector<double>& sigma)
      : _sigma(sigma), _delta(), _epsilon(0.1), _xi(1.0), _omega(1) {
    _n_orders = sigma.size();
    _delta.reserve(_n_orders);
    // Eq.(3.4)
    _delta.push_back(1.);  // normalised w.r.t. LO <-> sigma[0]
    for (auto i = 1; i < _n_orders; ++i) {
      _delta.push_back((sigma[i] - sigma[i - 1]) / sigma[0]);
    }
    std::cout << "# ABCModel: " << _n_orders << " \n";
    for (auto i = 0; i < sigma.size(); ++i) {
      std::cout << "# > " << sigma[i] << ", " << _sigma[i] << ":\t" <<
      _delta[i] << std::endl;
    }
  }

  double sigma(int order) const { return _sigma.at(order); };
  double pdf(const double& val) const;

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
  // parameters of the model
  double _epsilon;
  double _xi;
  int _omega;

 private:
  // duplicate the infrastructure for the brute-force 1D integration until I
  // have the analytic result
  double a_integrand(const double& a, const std::vector<double>& delta) const;
  double a_int(const std::vector<double>& delta) const;
};

}  // namespace miho
