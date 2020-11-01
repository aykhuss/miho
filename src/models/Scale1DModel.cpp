#include "Scale1DModel.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace miho {

double Scale1DModel::sigma(int order) const {
  // try to return the central scale choice
  for (const auto& mod : _scale_models) {
    if (mod.first.fac_mu == 1.) {
      return mod.second->sigma(order);
    }
  }
  // no central scale found: return the first one
  return _scale_models.begin()->second->sigma(order);
}

double Scale1DModel::pdf(const double& val) const {
  if (_use_gauss_legendre) {
    return pdf_gauss_legendre(val);
  } else {
    return pdf_trapezoid(val);
  }
}

double Scale1DModel::pdf_trapezoid(const double& val) const {
  // double log_min = log(_scale_models.begin()->first.fac_mu);
  // double log_max = log(_scale_models.rbegin()->first.fac_mu);
  // std::cout << "# Scale1DModel: marginalising log(mu) [" << log_min
  //           << ", " << log_max << "]\n";

  // the entries are already sorted by key
  double result = 0.;
  double log_mu_last = 0.;
  double pdf_val_last = 0.;
  for (auto it = _scale_models.begin(); it != _scale_models.end(); ++it) {
    double log_mu_curr = log(it->first.fac_mu);
    double pdf_val_curr = it->second->pdf(val);
    if (it != _scale_models.begin()) {
      double dlog_mu = log_mu_curr - log_mu_last;
      result += dlog_mu * (pdf_val_curr + pdf_val_last) / 2.;
    }
    log_mu_last = log_mu_curr;
    pdf_val_last = pdf_val_curr;
  }
  return result;
}

double Scale1DModel::pdf_gauss_legendre(const double& val) const {
  if (_scale_models.size() != 3)
    throw std::runtime_error(
        "Scale1DModel::pdf_gauss_legendre: incompatible # of modles");
  const std::vector<std::pair<double, double>> gauss_legendre_weights = {
      {0.5, 5. / 9}, {1., 8. / 9.}, {2., 5. / 9}};
  size_t i = 0;
  double result = 0.;
  for (const auto& scl_gm : _scale_models) {
    if (is_approx(scl_gm.first.fac_mu, gauss_legendre_weights.at(i).first)) {
      result += gauss_legendre_weights.at(i).second * scl_gm.second->pdf(val);
    } else {
      throw std::runtime_error(
          "Scale1DModel::pdf_gauss_legendre: invalid scale value?!");
    }
    i++;
  }
  return result;
}

}  // namespace miho
