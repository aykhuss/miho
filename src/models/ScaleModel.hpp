
namespace miho {

template <std::size_t N>
double ScaleModel<N>::sigma(int order) const {
  /// try to return the central scale choice
  for (const auto& mod : _scale_models) {
    if (std::all_of(mod.first.begin(), mod.first.end(),
                    [](const double& fac) { return is_approx(fac, 1.); })) {
      return mod.second->sigma(order);
    }
  }
  /// no central scale found: return the first one
  return _scale_models.begin()->second->sigma(order);
}

template <std::size_t N>
double ScaleModel<N>::pdf(const double& val) const {
  if (!_q_init) {
    init();
    _q_init = true;
  }
  switch (_prescription) {
    case ScalePrescription::scale_average:
      // std::cerr << "PDF: scale_average\n";
      return pdf_scale_average(val);
    case ScalePrescription::scale_marginalisation:
      // std::cerr << "PDF: scale_marginalisation\n";
      return pdf_scale_marginalisation(val);
    default:
      throw "ScaleModel::pdf: unknown prescription";
  }
}

template <std::size_t N>
double ScaleModel<N>::pdf_scale_average(const double& val) const {
  // const auto norm = int_mu([&val](const ModelPrototype& mod) { return 1.; });
  // std::cerr << "norm = " << norm << std::endl;
  return int_mu([&val](const ModelPrototype& mod) { return mod.pdf(val); });
}

template <std::size_t N>
double ScaleModel<N>::pdf_scale_marginalisation(const double& val) const {
  if (!_q_pdf_den) {
    _pdf_den =
        int_mu([](const ModelPrototype& mod) { return mod.pdf_delta__mu(); });
    _q_pdf_den = true;
  }
  const auto pdf_num = int_mu([&val](const ModelPrototype& mod) {
    return mod.pdf_delta__mu(mod.delta_next(val)) / std::fabs(mod.sigma(0));
  });
  return pdf_num / _pdf_den;
}

template <std::size_t N>
double ScaleModel<N>::int_mu(MPFun_t fun) const {
  auto it_map = _scale_models.begin();
  if (_use_gauss_legendre) {
    return int_mu_gauss_legendre(fun, it_map);
  } else {
    return int_mu_trapezoid(fun, it_map);
  }
}

template <std::size_t N>
double ScaleModel<N>::int_mu_trapezoid(
    MPFun_t fun,
    typename std::map<Scale_t, std::shared_ptr<ModelPrototype>>::const_iterator&
        it,
    int idim) const {
  /// ended up at a "leaf": return the function evaluation
  if (idim >= N) {
    const double ret = fun(*(it->second));
    ++it;
    return ret;
  }
  /// accumulate integral for dimension `idim`
  double result = 0.;
  double log_mu_last = 0.;
  double fun_val_last = 0.;
  for (auto iscl = 0; iscl < _fac.at(idim).size(); ++iscl) {
    // std::cerr << "> " << idim << "[" << iscl << "] " << it->first.at(idim)
    //           << " vs " << _fac.at(idim).at(iscl) << "\t";
    // std::cerr << "# ";
    // for (const auto& val : it->first) std::cerr << val << " ";
    // std::cerr << "\n";
    const double log_mu_curr = log(it->first.at(idim));
    const double fun_val_curr = int_mu_trapezoid(fun, it, idim + 1);
    if (iscl > 0) {
      const double dlog_mu = log_mu_curr - log_mu_last;
      result += dlog_mu * (fun_val_curr + fun_val_last) / 2.;
    }
    log_mu_last = log_mu_curr;
    fun_val_last = fun_val_curr;
  }
  return result / (log(_fac.at(idim).back()) - log(_fac.at(idim).front()));
}

template <std::size_t N>
double ScaleModel<N>::int_mu_gauss_legendre(
    MPFun_t fun,
    typename std::map<Scale_t, std::shared_ptr<ModelPrototype>>::const_iterator&
        it,
    int idim) const {
  constexpr std::array<double, 3> _weights_gauss_legendre = {5. / 9., 8. / 9.,
                                                             5. / 9.};
  /// ended up at a "leaf": return the function evaluation
  if (idim >= N) {
    const double ret = fun(*(it->second));
    ++it;
    return ret;
  }
  ///@todo: check _fac.at(idim).size() == 3 here?
  /// accumulate integral for dimension `idim`
  double result = 0.;
  for (auto iscl = 0; iscl < _fac.at(idim).size(); ++iscl) {
    // std::cerr << "> " << idim << "[" << iscl << "] " << it->first.at(idim)
    //           << " vs " << _fac.at(idim).at(iscl) << "\t";
    // std::cerr << "# ";
    // for (const auto& val : it->first) std::cerr << val << " ";
    // std::cerr << "\n";
    result += _weights_gauss_legendre[iscl] *
              int_mu_gauss_legendre(fun, it, idim + 1);
  }
  return result / 2.;  // 2 = sum(_weights_gauss_legendre): GL range was [-1,+1]
}

}  // namespace miho
