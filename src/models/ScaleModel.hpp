
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
  std::cerr << "\n##### pdf[" << val << "]\n";
  auto it_map = _scale_models.begin();
  auto res =
      int_mu_trapezoid([](const ModelPrototype& mod) { return mod.sigma(); }, it_map);
  std::cerr << "##### result = " << res << "\n\n";
  return res;
  if (_use_gauss_legendre) {
    return pdf_gauss_legendre(val);
  } else {
    return pdf_trapezoid(val);
  }
}

template <std::size_t N>
double ScaleModel<N>::int_mu_trapezoid(
    std::function<double(const ModelPrototype&)> fun,
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
  return result;
}


template <std::size_t N>
double ScaleModel<N>::int_gauss_legendre(
    std::function<double(const ModelPrototype&)> fun,
    typename std::map<Scale_t, std::shared_ptr<ModelPrototype>>::const_iterator&
        it,
    int idim) const {
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
    result +=
        _weights_gauss_legendre[iscl] * int_gauss_legendre(fun, it, idim + 1);
  }
  return result;
}


template <std::size_t N>
double ScaleModel<N>::pdf_trapezoid(const double& val) const {
  return 0.;
}

template <std::size_t N>
double ScaleModel<N>::pdf_gauss_legendre(const double& val) const {
  return 0.;
}

}  // namespace miho
