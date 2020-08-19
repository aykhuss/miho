#include <fmt/format.h>

#include <CLI/CLI.hpp>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <regex>

#include "Model.h"
#include "GeometricModel.h"
#include "Scale1DGeometricModel.h"

int main(int argc, char const* argv[]) {
  CLI::App app{"ミホ (miho) - theory uncertainties from MIssing Higher Orders"};

  app.require_subcommand(/* min */ 0, /* max */ 1);

  CLI::App* app_gm = app.add_subcommand("geometric", "Use the GeometricModel.");

  CLI::App* app_scl1gm = app.add_subcommand("scale1d", "Use the Scale1DGeometricModel.");
  std::vector<double> scl1_vec;
  app_scl1gm->add_option("--scales", scl1_vec, "The scale factors.");
  bool scl1_gl = false;
  app_scl1gm->add_flag("--gauss-legendre", scl1_gl, "Use the 3-point Gauss–Legendre quadrature.");

  bool flag_bool;
  app.add_flag("--bool,-b", flag_bool, "This is a bool flag");

  CLI11_PARSE(app, argc, argv);

  std::shared_ptr<miho::Model> cli_model = nullptr;

  if (app.got_subcommand(app_gm)) {
    fmt::print("# GeometricModel\n");

    std::string line;
    double sig;
    std::vector<double> sigma;

    fmt::print("# Enter XS values separated by spaces: \n");
    std::getline(std::cin, line);
    std::istringstream stream(line);
    std::copy(std::istream_iterator<double>(stream), std::istream_iterator<double>(),
              std::back_inserter(sigma));

    std::cout << "# Numbers you entered: ";
    std::copy(sigma.begin(), sigma.end(),
              std::ostream_iterator<double>(std::cout, " "));
    std::cout << '\n';

    miho::GeometricModel test_gm(sigma);
    // double val_low = 0.;
    // double val_upp = 70.;
    // size_t n_steps = 700;
    // for (auto i = 0; i <= n_steps; ++i) {
    //   double val = val_low + i * (val_upp - val_low) / n_steps;
    //   double xs = test_gm.pdf(val);
    //   std::cout << val << "\t" << xs << std::endl;
    // }
    //double norm = test_gm.integrate([](double) { return 1.; });

    cli_model = std::shared_ptr<miho::Model>(new miho::GeometricModel(sigma));
  }

  if (app.got_subcommand(app_scl1gm)) {
    fmt::print("# Scale1DGeometricModel\n");

    std::cout << "# Scales you entered: ";
    std::copy(scl1_vec.begin(), scl1_vec.end(),
              std::ostream_iterator<double>(std::cout, " "));
    std::cout << '\n';

    // cli_model = std::shared_ptr<miho::Scale1DGeometricModel>(
    //     new miho::Scale1DGeometricModel());
    std::shared_ptr<miho::Scale1DGeometricModel> scl1gm =
        std::shared_ptr<miho::Scale1DGeometricModel>(
            new miho::Scale1DGeometricModel());

    scl1gm->use_gauss_legendre(scl1_gl);

    // std::string lline;
    // while (std::getline(std::cin, lline)) {
    //   if (std::regex_match(lline, std::regex("^[[:space:]]*#.*"))) continue;
    //   if (std::regex_match(lline, std::regex("^[[:space:]]*$"))) break;
    //   std::cout << "$ |" << lline << "|" << std::endl;
    // }
    // return 0;


    for (const auto& scl : scl1_vec) {
      std::string line;
      double sig;
      std::vector<double> sigma;
      fmt::print("# Enter XS[{}×μ₀] values separated by spaces: \n", scl);
      //line.clear();  // not needed
      std::getline(std::cin, line);
      //std::cout << "line = " << line << std::endl;
      std::istringstream stream(line);
      //std::cout << "stream = " << stream.str() << std::endl;
      std::copy(std::istream_iterator<double>(stream),
                std::istream_iterator<double>(), std::back_inserter(sigma));

      std::cout << "# Numbers you entered: ";
      std::copy(sigma.begin(), sigma.end(),
                std::ostream_iterator<double>(std::cout, " "));
      std::cout << '\n';

      scl1gm->add_model(scl, miho::GeometricModel(sigma));
    }

    cli_model = scl1gm;
  }

  //>>> do something with cli_model <<<

  // double val_low = 450.;
  // double val_upp = 550.;
  // size_t n_steps = 400;
  // // double val_low = 0.;
  // // double val_upp = 70.;
  // // size_t n_steps = 700;
  // for (auto i = 0; i <= n_steps; ++i) {
  //   double val = val_low + i * (val_upp - val_low) / n_steps;
  //   double xs = cli_model->pdf(val);
  //   std::cout << val << "\t" << xs << std::endl;
  // }

  cli_model->print_nodes();

  double norm = cli_model->integrate([](double) { return 1.; });
  // fmt::print("# normalisation: {}", norm);

  double median = cli_model->median();
  std::pair<double, double> DoB68 = cli_model->degree_of_belief_interval();
  std::pair<double, double> DoB95 = cli_model->degree_of_belief_interval(0.95);
  double mean   = cli_model->mean();
  fmt::print("#result: {}  {}  {} {}  {} {} \n", mean, median, DoB68.first,DoB68.second, DoB95.first,DoB95.second);



  return 0;


//
//   // miho::GeometricModel gm({12.9953, 30.705, 41.8182, 46.2879});
//   // miho::GeometricModel gm0({12.9953});
//   // miho::GeometricModel gm1({12.9953, 30.705});
//   // miho::GeometricModel gm2({12.9953, 30.705, 41.8182});
//   // miho::GeometricModel gm3({12.9953, 30.705, 41.8182, 46.2879});
//
//   //!--- extract from Fig. 4.4 left
//   //! -0.07523526321881069, 16.04463991595518
//   //!  0.9262960621789333, 36.734923418891704  -> 68%
//   //!  upper: 1.06598498944988,   47.1450871460665
//   //!  1.9284817620925736, 46.3248338838389    -> 68%
//   //!  upper: 2.0685134569670853, 50.92051649587643
//   //!  2.927008306549977,  47.98621807142033   -> 68%
//   //!  upper: 3.0672447716811635, 49.10831456272649
//   miho::GeometricModel gm({16.04463991595518, 36.734923418891704,
//                            46.3248338838389, 47.98621807142033});
//   miho::GeometricModel gm0({16.04463991595518});
//   miho::GeometricModel gm1({16.04463991595518, 36.734923418891704});
//   miho::GeometricModel gm2(
//       {16.04463991595518, 36.734923418891704, 46.3248338838389});
//   miho::GeometricModel gm3({16.04463991595518, 36.734923418891704,
//                             46.3248338838389, 47.98621807142033});
//
//   // double val_low = 0.;
//   // double val_upp = 70.;
//   // size_t n_steps = 700;
//   // for (auto i = 0; i <= n_steps; ++i) {
//   //   double val = val_low + i * (val_upp - val_low) / n_steps;
//   //   // double xs = gm.pdf(val);
//   //   // std::cout << val << "\t" << xs << std::endl;
//   //   std::cout << val << "\t" << gm0.pdf(val) << "\t" << gm1.pdf(val) << "\t"
//   //             << gm2.pdf(val) << "\t" << gm3.pdf(val) << std::endl;
//   // }
//
//   //double norm = gm.integrate([](double) { return 1.; });
//   // double E = gm.integrate([](double x) { return x; });
//   // double E2 = gm.integrate([](double x) { return x*x; });
//
//   return 0;
}
