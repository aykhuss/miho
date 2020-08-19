#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <limits>

#include "GeometricModel.h"

miho::GeometricModel gm({1., 2., 2.5, 2.6});

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

TEST_CASE("binomial", "[GeometricModel]") {
  REQUIRE(miho::GeometricModel::binomial(6, 3) == 20);
  REQUIRE(miho::GeometricModel::binomial(6, 4) == 15);
  REQUIRE(miho::GeometricModel::binomial(7, 2) == 21);
}
