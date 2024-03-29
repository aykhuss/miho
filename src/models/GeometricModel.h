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
    clear();  // clear caches
  }
  inline void set_omega(int omega) {
    _omega = omega;
    clear();  // clear caches
  }

 protected:
  /// parameters of the model
  int _omega = 1;         // Eq.(4.9) : assumes integer
  double _epsilon = 0.1;  // Eq.(4.8)
  /// additional helper functions
  static double a_integral(const double& a_lower, const double& a_upper,
                           const double& epsilon, int m, int k, int j);
};

}  // namespace miho
