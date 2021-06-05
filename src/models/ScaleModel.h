#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>

#include "Model.h"
#include "ModelPrototype.h"

namespace miho {

enum class ScalePrescription { scale_average, scale_marginalisation };

/// would this be a better abstraction? maybe some other time...
// template <std::size_t N>
// class ScaleType {
//  public:
//   ScaleType() : _fac() {}
//   /// needed for map sorting
//   bool operator<(const ScaleType<N>& other) const;
//  private:
//   std::array<double, N> _fac;
// };

template <std::size_t N>
class ScaleModel : public Model {
  using Scale_t = std::array<double, N>;
  using MPFun_t = std::function<double(const ModelPrototype&)>;

 public:
  using Model::sigma;  /// avoid name hiding

  void add_model(const Scale_t& fac_mu, std::shared_ptr<ModelPrototype> model) {
    // std::cout << "# add_model: _n_orders = " << _n_orders << std::endl;
    _scale_models.insert(std::pair<Scale_t, std::shared_ptr<ModelPrototype>>(
        Scale_t(fac_mu), model));
    if (_n_orders == 0) {
      _n_orders = model->n_orders();
    } else if (_n_orders != model->n_orders()) {
      throw std::runtime_error("ScaleModel::add_model: order mismatch!");
    }
    clear();
  }

  void init() const {
    // std::cerr << "=== init ===\n";
    int count = 0;
    for (const auto& mod : _scale_models) {
      // ///===
      // std::cerr << count << ": ";
      // for (const auto& val : mod.first) {
      //   std::cerr << val << " ";
      // }
      // std::cerr << "\n";
      // count++;
      // ///===
      for (auto idim = 0; idim < mod.first.size(); ++idim) {
        const double& val = mod.first.at(idim);
        auto& ivec_mu = _fac.at(idim);
        if (!ivec_mu.empty()) {
          if (is_approx(ivec_mu.back(), val) || (ivec_mu.back() > val))
            continue;
        }
        ivec_mu.push_back(val);
      }
      ///===
    }
    // std::cerr << "=== reconstructed ===\n";
    // for (const auto& fac : _fac) {
    //   ///===
    //   std::cerr << "> ";
    //   for (const auto& val : fac) {
    //     std::cerr << val << " ";
    //   }
    //   std::cerr << "\n";
    // }
  }

  double sigma(int order) const override;
  double pdf(const double& val) const override;

  /// setters
  inline void use_gauss_legendre(bool flag) {
    _use_gauss_legendre = flag;
    clear();  // clear caches
    ///@todo if GL check for compatibility or "complete the the square"
  }
  inline void set_prescription(ScalePrescription prescr) {
    _prescription = prescr;
    clear();  // clear caches
  }

 protected:
  /// cache denominator
  mutable bool _q_init = false;
  mutable std::array<std::vector<double>, N> _fac;
  mutable bool _q_pdf_den = false;
  mutable double _pdf_den;
  virtual void clear() override {
    _q_init = false;
    _q_pdf_den = false;
    Model::clear();
  }

 private:
  bool _use_gauss_legendre = false;
  ScalePrescription _prescription = ScalePrescription::scale_average;

  /// from C++20 on, no operator< so would need to pass to the map:
  /// std::lexicographical_compare
  std::map<Scale_t, std::shared_ptr<ModelPrototype>> _scale_models;

  double pdf_scale_average(const double& val) const;
  double pdf_scale_marginalisation(const double& val) const;

  /// integration over the scales(s)
  double int_mu(MPFun_t fun) const;

  double int_mu_trapezoid(
      MPFun_t fun,
      typename std::map<Scale_t,
                        std::shared_ptr<ModelPrototype>>::const_iterator& it,
      int idim = 0) const;

  double int_mu_gauss_legendre(
      std::function<double(const ModelPrototype&)> fun,
      typename std::map<Scale_t,
                        std::shared_ptr<ModelPrototype>>::const_iterator& it,
      int idim = 0) const;
};

}  // namespace miho

#include "ScaleModel.hpp"
