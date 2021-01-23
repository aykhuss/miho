#pragma once

#include "ModelPrototype.h"

namespace miho {

class ABCNumericModel : public ModelPrototype {
 public:
  std::unique_ptr<ModelPrototype> clone() const override {
    return std::make_unique<ABCNumericModel>(*this);
  }

  double pdf_delta__mu(const std::vector<double>& delta) const override;

  /// setters
  inline void set_epsilon(const double& epsilon) {
    _epsilon = epsilon;
    clear();  // clear caches
  }
  inline void set_eta(const double& eta) {
    _eta = eta;
    clear();  // clear caches
  }
  inline void set_xi(const double& xi) {
    _xi = xi;
    clear();  // clear caches
  }
  inline void set_omega(int omega) {
    _omega = omega;
    clear();  // clear caches
  }
  inline void set_nint_rel_err(const double& acc) { _nint_rel_err = acc; }

 protected:
  /// parameters of the model
  int _omega = 1;
  double _xi = 5.0;
  double _epsilon = 0.1;
  double _eta = 0.2;
  /// parameter for the cubature integration
  double _nint_rel_err = 0.01;  // 1%
};

}  // namespace miho
