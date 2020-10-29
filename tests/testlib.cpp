#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <limits>

#include "Util.h"
#include "GeometricModel.h"
#include "Model.h"

miho::GeometricModel gm({1., 2., 2.5, 2.6});

TEST_CASE("binomial", "[Util]") {
  REQUIRE(miho::binomial(6, 3) == 20);
  REQUIRE(miho::binomial(6, 4) == 15);
  REQUIRE(miho::binomial(7, 2) == 21);
}

TEST_CASE("factorial", "[Util]") {
  REQUIRE(miho::factorial(0) == 1);
  REQUIRE(miho::factorial(1) == 1);
  REQUIRE(miho::factorial(7) == 5040);
}

TEST_CASE("pochhammer", "[Util]") {
  int output_1 = miho::pochhammer<int>(42, 5);
  int check_1 = 164490480;
  REQUIRE(output_1 == check_1);

  double output_2 = miho::pochhammer<double>(4.2, 7);
  double check_2 = 750203.;
  REQUIRE(output_2 == Approx(check_2));

  double output_3 = miho::pochhammer<double>(-1.1, 3);
  double check_3 = 0.099;
  REQUIRE(output_3 == Approx(check_3));
}

TEST_CASE("hypergeometric_2F1", "[Util]") {
  double output_1 = miho::hypergeometric_2F1(1, 2, 3, 0.1);
  double check_1 = 1.0721;
  REQUIRE(output_1 == Approx(check_1));

  // output_1 = r2f1_adsj(1, 2, 3, 0.1);
  // REQUIRE(output_1 == Approx(check_1));

  double output_2 = miho::hypergeometric_2F1(0.1, 3.3, 4.2, 0.9);
  double check_2 = 1.14889;
  REQUIRE(output_2 == Approx(check_2));

  // output_2 = r2f1_adsj(0.1, 3.3, 4.2, 0.9);
  // REQUIRE(output_2 == Approx(check_2));


  double output_3 = miho::hypergeometric_2F1(11, 7, 5.3, -0.4);
  double check_3 = 0.0048317;
  REQUIRE(output_3 == Approx(check_3));

  // output_3 = r2f1_adsj(11, 7, 5.3, -0.4);
  // REQUIRE(output_3 == Approx(check_3));


  double output_4 = miho::hypergeometric_2F1(3, 8, 5.2, -0.3);
  double check_4 = 0.307245;
  REQUIRE(output_4 == Approx(check_4));

  // output_4 = r2f1_adsj(3, 8, 5.2, -0.3);
  // REQUIRE(output_4 == Approx(check_4));


  double output_5 = miho::hypergeometric_2F1(2, 5, 9, -0.8);
  double check_5 = 0.490718;
  REQUIRE(output_5 == Approx(check_5));

  // output_5 = r2f1_adsj(2, 5, 9, -0.8);
  // REQUIRE(output_5 == Approx(check_5));

}

TEST_CASE("appell_f1", "[Util]") {
  double output_1 = miho::appell_F1(1, 2, 3, 4, 0.1, 0.2);
  double check_1 = 1.24827;
  REQUIRE(output_1 == Approx(check_1));

  // output_1 = rf1_adsj(1, 2, 3, 4, 0.1, 0.2);
  // REQUIRE(output_1 == Approx(check_1));

  double output_2 = miho::appell_F1(0.1, 3.3, 4.2, 44, 0.9, 0.01);
  double check_2 = 1.0072;
  REQUIRE(output_2 == Approx(check_2));

  // output_2 = rf1_adsj(0.1, 3.3, 4.2, 44, 0.9, 0.01);
  // REQUIRE(output_2 == Approx(check_2));

  double output_3 = miho::appell_F1(11, 7, 5.3, 2.1, -0.5, 0.75);
  double check_3 = 5.65353e+8;
  REQUIRE(output_3 == Approx(check_3));

  // output_3 = rf1_adsj(11, 7, 5.3, 2.1, -0.5, 0.75);
  // REQUIRE(output_3 == Approx(check_3));
}

TEST_CASE("limits", "[Model]") {
  REQUIRE( (std::numeric_limits<double>::infinity() > 42.e+42 ) == true);
}


TEST_CASE("a list", "[GeometricModel]") {
  double inf = std::numeric_limits<double>::infinity();

  auto output_1 = miho::GeometricModel::a_list({1., 0.5, 0.1, 0.01});
  REQUIRE_THAT(output_1, Catch::Approx<double>({inf, 0.5, 0.2, 0.1, 0.}));

  auto output_2 = miho::GeometricModel::a_list({1., 0.01, 0.5, 0.1});
  REQUIRE_THAT(output_2,
               Catch::Approx<double>({inf, 0.707107, 0.707107, 0.2, 0.}));

  auto output_3 = miho::GeometricModel::a_list({1., 0.01, 0.1, 0.5});
  REQUIRE_THAT(output_3,
               Catch::Approx<double>({inf, 0.793701, 0.793701, 0.793701, 0.}));
}
