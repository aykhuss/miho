#pragma once

#include <functional>
#include <map>
#include <vector>

// extern "C" double r2f1_adsj(const double& a, const double& b, const double& c,
//                             const double& z);
// extern "C" double rf1_adsj(const double& a, const double& b1, const double& b2,
//                            const double& c, const double& x, const double& y);

namespace miho {

template <typename T>
int sgn(const T& val);

bool is_approx(const double& x, const double& y);

bool is_approx_int(const double& x, const double& y);

int binomial(int n, int k);

int factorial(int n);

template <typename T>
T pochhammer(const T& val, int n);

double hypergeometric_2F1(const double& a, const double& b, const double& c,
                          const double& z);

double appell_F1(const double& a, const double& b1, const double& b2,
                 const double& c, const double& z1, const double& z2);

double nintegrate_1D(std::function<double(double)> func, const double& x_low,
                     const double& x_upp);

double erf_inv(const double& x);


/// input: list of delta: d[0], ..., d[n]
/// output: list of transitions: infty == a[0] >= ... >= a[n+1] == 0
///         for max{|d[0]|, |d[1]|/a, ..., |d[n]|/a^n}
///   for inverted case:
///       : list of transitions: 0 == a[0] >= ... >= a[n+1] == infty
///         for min{|d[0]|, |d[1]|/a, ..., |d[n]|/a^n}
std::vector<double> a_list(const std::vector<double>& delta,
                           bool invert = false);
/// variants without the absolute value for min/max cases
std::vector<double> max_a_list(std::vector<double> delta);
std::vector<double> min_a_list(std::vector<double> delta);
/// get the powers and ranges for min/max
std::map<std::pair<int, int>, std::pair<double, double>> min_max_partitions(
    const std::vector<double>& delta, const double& max_val = 1.);
/// find min/max value for `delta_k / a^k`
std::pair<double, double> get_min_max(const double& a,
                                      const std::vector<double>& delta);
std::pair<double, double> get_min_max_check(const double& a,
                                            const std::vector<double>& delta);

/// simple accessor routines for random numbers to do some testing
/// not threadsafe!
void set_rand_seed(int seed_number);
double rand();

}  // namespace miho
