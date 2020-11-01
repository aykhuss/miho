#include "Scale2DModel.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace miho {

double Scale2DModel::sigma(int order) const {
  // try to return the central scale choice
  for (const auto& mod : _scale_models) {
    if (mod.first.fac_muR == 1. && mod.first.fac_muF == 1.) {
      return mod.second->sigma(order);
    }
  }
  // no central scale found: return the first one
  return _scale_models.begin()->second->sigma(order);
}

double Scale2DModel::pdf(const double& val) const {
  if (_use_gauss_legendre) {
    return pdf_gauss_legendre(val);
  } else {
    return pdf_trapezoid(val);
  }
}

double Scale2DModel::pdf_trapezoid(const double& val) const {
  // @todo: disgusting implementation; generalise to arbitrary N dimensions
  // with generic templated class...

  double last_logR = 0.;
  double curr_logR = 0.;
  double last_logF = 0.;
  double curr_logF = 0.;
  bool is_head_R = true;
  bool is_head_F = true;

  double last_pdf_num = 0.;
  double curr_pdf_num = 0.;

  double last_intR_num = 0.;
  double curr_intR_num = 0.;

  double result = 0.;

  for (auto it = _scale_models.begin(); it != _scale_models.end(); ++it) {
    curr_logR = log(it->first.fac_muR);
    curr_logF = log(it->first.fac_muF);
    curr_pdf_num = it->second->pdf(val);
    // std::cout << is_head_R << is_head_F << "muR: " << it->first.fac_muR << ",
    // muF: " << it->first.fac_muF
    //           << "  [" << curr_pdf_num << "," << curr_pdf_den << "]\n";
    if (!is_head_R) {
      double dlogR = curr_logR - last_logR;
      curr_intR_num += dlogR * (curr_pdf_num + last_pdf_num) / 2.;
      // std::cout << "---" << curr_intR_den << std::endl;
    }
    is_head_R = false;
    // reached the end of a muF slice
    if (!is_approx(std::next(it)->first.fac_muF, it->first.fac_muF)) {
      // accumulate
      if (!is_head_F) {
        double dlogF = curr_logF - last_logF;
        result += dlogF * (curr_intR_num + last_intR_num) / 2.;
      }
      is_head_F = false;
      // prepare for next muF slice...
      is_head_R = true;
      last_logF = curr_logF;
      last_intR_num = curr_intR_num;
      curr_intR_num = 0.;
      curr_pdf_num = 0.;  // for safe measure
    }
    last_logR = curr_logR;
    last_pdf_num = curr_pdf_num;
  }
  // std::cout << "result = " << result << std::endl;
  return result;
}

double Scale2DModel::pdf_gauss_legendre(const double& val) const {
  if (_scale_models.size() != 9) {
    std::cerr << "Scale2DModel::pdf_gauss_legendre: incompatible # of modles\n";
    return 0.;
  }
  const std::vector<std::pair<double, double>> GL_weights = {
      {0.5, 5. / 9}, {1., 8. / 9.}, {2., 5. / 9}};
  size_t iR = 0;
  size_t iF = 0;
  double result = 0.;
  for (const auto& scl_gm : _scale_models) {
    // std::cout << "Scale2DModel::GL " << iR << " & " << iF << std::endl;
    if (is_approx(scl_gm.first.fac_muR, GL_weights.at(iR).first) &&
        is_approx(scl_gm.first.fac_muF, GL_weights.at(iF).first)) {
      result += GL_weights.at(iR).second * GL_weights.at(iF).second *
                scl_gm.second->pdf(val);
    } else {
      throw std::runtime_error(
          "Scale2DModel::pdf_gauss_legendre: invalid scale value?!");
    }
    iR++;
    if (iR > 2) {
      iR = 0;
      iF++;
    }
  }
  return result;
}

}  // namespace miho
