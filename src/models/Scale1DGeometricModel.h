#pragma once

#include <map>

#include "GeometricModel.h"
#include "Model.h"

namespace miho {

class Scale1DGeometricModel : public Model {
 public:
  Scale1DGeometricModel() : _use_gauss_legendre(false), _scale_models() {}
  double sigma(int order) const;
  double pdf(const double& val) const;
  // use Move eventually!
  void add_model(const double& fac_mu, const GeometricModel& gm) {
    std::cout << "# add_model: _n_orders = " << _n_orders << std::endl;
    _scale_models.insert(
        std::pair<Scale1D, GeometricModel>(Scale1D(fac_mu), gm));
    if (_n_orders == 0) {
      _n_orders = gm.n_orders();
    } else if (_n_orders != gm.n_orders()) {
      throw std::runtime_error(
          "Scale1DGeometricModel::add_model: order mismatch!");
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
  std::map<Scale1D, GeometricModel> _scale_models;
};

}  // namespace miho
