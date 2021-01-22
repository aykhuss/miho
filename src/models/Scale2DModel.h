#pragma once

#include <map>

#include "Model.h"
#include "ModelPrototype.h"

namespace miho {

class Scale2DModel : public Model {
 public:
  Scale2DModel() : _use_gauss_legendre(false), _scale_models() {}
  double sigma(int order) const override;
  double pdf(const double& val) const override;
  // use Move eventually!
  void add_model(const std::pair<double, double> fac_mu,
                 std::shared_ptr<miho::ModelPrototype> model) {
    // std::cout << "# add_model: _n_orders = " << _n_orders << std::endl;
    _scale_models.insert(
        std::pair<Scale2D, std::shared_ptr<miho::ModelPrototype>>(
            Scale2D(fac_mu.first, fac_mu.second), model));
    if (_n_orders == 0) {
      _n_orders = model->n_orders();
    } else if (_n_orders != model->n_orders()) {
      throw std::runtime_error("Scale2DModel::add_model: order mismatch!");
    }
  }
  inline void use_gauss_legendre(bool flag) { _use_gauss_legendre = flag; }

 private:
  struct Scale2D {
    double fac_muR;
    double fac_muF;
    Scale2D(const double& valR, const double& valF)
        : fac_muR(valR), fac_muF(valF) {}
    bool operator<(const Scale2D& other) const {
      // first sort in muF, then in muR
      // => first loop over muR, then muF
      if (fac_muF != other.fac_muF) {
        return fac_muF < other.fac_muF;
      } else {
        return fac_muR < other.fac_muR;
      }
    }
  };

  double pdf_trapezoid(const double& val) const;
  double pdf_gauss_legendre(const double& val) const;

  bool _use_gauss_legendre;
  std::map<Scale2D, std::shared_ptr<miho::ModelPrototype>> _scale_models;
};

}  // namespace miho
