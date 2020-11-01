#pragma once

#include <map>

#include "Model.h"

namespace miho {

class Scale1DModel : public Model {
 public:
  Scale1DModel() : _use_gauss_legendre(false), _scale_models() {}
  double sigma(int order) const;
  double pdf(const double& val) const;
  // use Move eventually!
  void add_model(
      const double& fac_mu, std::shared_ptr<miho::Model> model) {
    std::cout << "# add_model: _n_orders = " << _n_orders << std::endl;
    _scale_models.insert(
        std::pair<Scale1D, std::shared_ptr<miho::Model>>(Scale1D(fac_mu), model));
    if (_n_orders == 0) {
      _n_orders = model->n_orders();
    } else if (_n_orders != model->n_orders()) {
      throw std::runtime_error("Scale1DModel::add_model: order mismatch!");
    }
  }
  inline void use_gauss_legendre(bool flag) { _use_gauss_legendre = flag; }

 private:
  struct Scale1D {
    double fac_mu;
    Scale1D(const double& val) : fac_mu(val) {}
    bool operator<(const Scale1D& other) const { return fac_mu < other.fac_mu; }
  };

  double pdf_trapezoid(const double& val) const;
  double pdf_gauss_legendre(const double& val) const;

  bool _use_gauss_legendre;
  std::map<Scale1D, std::shared_ptr<miho::Model>> _scale_models;
};

}  // namespace miho
