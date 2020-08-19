#include "Scale2DGeometricModel.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace miho {

double Scale2DGeometricModel::sigma(int order) const {
  // try to return the central scale choice
  for (const auto& mod : _scale_models) {
    if (mod.first.fac_mu == 1.) {
      return mod.second.sigma(order);
    }
  }
  // no central scale found: return the first one
  return _scale_models.begin()->second.sigma(order);
}

double Scale2DGeometricModel::pdf(const double& val) const {
  if (_use_gauss_legendre) {
    return pdf_gauss_legendre(val);
  } else {
    return pdf_trapezoid(val);
  }
}

double Scale2DGeometricModel::pdf_trapezoid(const double& val) const {
  // double log_min = log(_scale_models.begin()->first.fac_mu);
  // double log_max = log(_scale_models.rbegin()->first.fac_mu);
  // std::cout << "# Scale2DGeometricModel: marginalising log(mu) [" << log_min
  //           << ", " << log_max << "]\n";

  // the entries are already sorted by key
  double result_num = 0.;
  double result_den = 0.;
  double log_mu_last = 0.;
  double pdf_last = 0.;
  double pdf_val_last = 0.;
  for (auto it = _scale_models.begin(); it != _scale_models.end(); ++it) {
    double log_mu_curr = log(it->first.fac_mu);
    double pdf_curr = it->second.pdf_delta__mu();
    double pdf_val_curr = it->second.pdf_delta__mu(it->second.delta_next(val)) /
                          it->second.sigma(0);
    if (it != _scale_models.begin()) {
      double dlog_mu = log_mu_curr - log_mu_last;
      result_num += dlog_mu * (pdf_val_curr + pdf_val_last) / 2.;
      result_den += dlog_mu * (pdf_curr + pdf_last) / 2.;
      // std::cout << "# > " << it->first.fac_mu << ": " << log_mu_curr << " & "
      //           << log_mu_last << ": " << pdf_curr << " & " << pdf_last
      //           << std::endl;
    }
    log_mu_last = log_mu_curr;
    pdf_last = pdf_curr;
    pdf_val_last = pdf_val_curr;
  }
  return result_num / result_den;
}

double Scale2DGeometricModel::pdf_gauss_legendre(const double& val) const {
  if (_scale_models.size() != 3)
    throw std::runtime_error(
        "Scale2DGeometricModel::pdf_gauss_legendre: incompatible # of modles");
  const std::vector<std::pair<double, double>> gauss_legendre_weights = {
      {0.5, 5. / 9}, {1., 8. / 9.}, {2., 5. / 9}};
  size_t i = 0;
  double result_num = 0.;
  double result_den = 0.;
  for (const auto& scl_gm : _scale_models) {
    if (is_approx(scl_gm.first.fac_mu, gauss_legendre_weights.at(i).first)) {
      result_num += gauss_legendre_weights.at(i).second *
                    scl_gm.second.pdf_delta__mu(scl_gm.second.delta_next(val)) /
                    scl_gm.second.sigma(0);
      result_den +=
          gauss_legendre_weights.at(i).second * scl_gm.second.pdf_delta__mu();
    } else {
      throw std::runtime_error(
          "Scale2DGeometricModel::pdf_gauss_legendre: invalid scale value?!");
    }
    i++;
  }
  return result_num / result_den;
}

}  // namespace miho
