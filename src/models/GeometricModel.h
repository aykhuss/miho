#pragma once

#include <iostream>
#include <vector>

#include "ModelPrototype.h"

namespace miho {

class GeometricModel : public ModelPrototype {
 public:
  std::unique_ptr<ModelPrototype> clone() const override {
    return std::make_unique<GeometricModel>(*this);
  }

  double pdf_delta__mu(const std::vector<double>& delta) const override;

  /// setters
  inline void set_epsilon(const double& epsilon) {
    _epsilon = epsilon;
    _q_pdf_den = false;
  }
  inline void set_omega(int omega) {
    _omega = omega;
    _q_pdf_den = false;
  }

  /// additional public member functions
  static std::vector<double> a_list(const std::vector<double>& delta);
  static double a_integral(const double& a_lower, const double& a_upper,
                           const double& epsilon, int m, int k, int j);

 protected:
  /// parameters of the model
  int _omega = 1;         // Eq.(4.9) : assumes integer
  double _epsilon = 0.1;  // Eq.(4.8)
};

}  // namespace miho
