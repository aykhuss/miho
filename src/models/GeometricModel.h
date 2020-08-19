#pragma once

#include <iostream>
#include <vector>

#include "Model.h"

namespace miho {

class GeometricModel : public Model {
 public:
  GeometricModel(const std::vector<double>& sigma)
      : _sigma(sigma), _delta(), _epsilon(0.1), _omega(1) {
    _n_orders = sigma.size();
    _delta.reserve(_n_orders);
    // Eq.(3.4)
    _delta.push_back(1.);  // normalised w.r.t. LO <-> sigma[0]
    for (auto i = 1; i < _n_orders; ++i) {
      _delta.push_back((sigma[i] - sigma[i - 1]) / sigma[0]);
    }
    // std::cout << "# GeometricModel: " << _n_orders << " \n";
    for (auto i=0; i<sigma.size();++i) {
      // std::cout << "# > " << sigma[i] << ", " << _sigma[i] << ":\t" <<
      // _delta[i] << std::endl;
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
  inline double pdf_delta__mu() const {
    return pdf_delta__mu(_delta);
  }
  double pdf_delta__mu(const double& delta_next) const;

  static std::vector<double> a_list(const std::vector<double>& delta);
  static int binomial(int n, int k);
  static double a_integral(const double& a_lower, const double& a_upper,
                           const double& epsilon, int m, int k, int j);

 protected:
  std::vector<double> _sigma;
  std::vector<double> _delta;
  // parameters of the model
  double _epsilon;  // Eq.(4.8)
  int _omega;       // Eq.(4.9) : assumes integer
};

}  // namespace miho
