#include <fmt/format.h>

#include <CLI/CLI.hpp>
#include <cmath>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

#include "GeometricModel.h"
#include "Model.h"
#include "Scale1DGeometricModel.h"
#include "Scale2DGeometricModel.h"


void print_format(std::string, std::shared_ptr<miho::Model>);

int main(int argc, char const* argv[]) {
  CLI::App app{"ミホ (miho) - theory uncertainties from MIssing Higher Orders"};
  app.footer("[powered by CERN QCD coffee]");

  // global settings for all models
  bool flag_pdf = false;
  app.add_flag("--pdf", flag_pdf, "Print out the probability distribution.");
  std::string format_string{"{median} {dob68_low} {dob68_upp} \n"};
  app.add_option("--format", format_string, "Output formatting string.");
  size_t nmax = 1000;
  app.add_option("--nmax", nmax, "Set the maximum number of PDF evaliations.");
  double accuracy = 0.001;  // default: 0.1%
  app.add_option("--accuracy", accuracy, "Set the target accuracy of the integration.");

  app.require_subcommand(/* min */ 0, /* max */ 1);

  CLI::App* app_gm = app.add_subcommand("geometric", "Use the GeometricModel.");

  CLI::App* app_scl1gm = app.add_subcommand(
      "scale1d", "Use the GeometricModel marginalising over one scale.");
  std::vector<double> scl1_vec;
  app_scl1gm->add_option("--scales", scl1_vec, "The scale factors.")
      ->required();
  bool scl1_gl = false;
  app_scl1gm->add_flag("--gauss-legendre", scl1_gl,
                       "Use the 3-point Gauss–Legendre quadrature.");

  CLI::App* app_scl2gm = app.add_subcommand(
      "scale2d", "Use the GeometricModel marginalising over two scales.");
  std::vector<std::pair<double, double>> scl2_vec;
  app_scl2gm
      ->add_option("--scales", scl2_vec,
                   "The scale factors as pairs.")
      ->required();
  bool scl2_gl = false;
  app_scl2gm->add_flag("--gauss-legendre", scl2_gl,
                       "Use the 3-point Gauss–Legendre quadrature.");

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
    std::copy(std::istream_iterator<double>(stream),
              std::istream_iterator<double>(), std::back_inserter(sigma));

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
    // double norm = test_gm.integrate([](double) { return 1.; });

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
      // line.clear();  // not needed
      std::getline(std::cin, line);
      // std::cout << "line = " << line << std::endl;
      std::istringstream stream(line);
      // std::cout << "stream = " << stream.str() << std::endl;
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

  if (app.got_subcommand(app_scl2gm)) {
    fmt::print("# Scale2DGeometricModel\n");

    std::cout << "# Scales you entered: ";
    for (const auto& scl : scl2_vec) {
      fmt::print("[{},{}] ", scl.first, scl.second);
    }
    fmt::print("\n");

    std::shared_ptr<miho::Scale2DGeometricModel> scl2gm =
        std::shared_ptr<miho::Scale2DGeometricModel>(
            new miho::Scale2DGeometricModel());

    scl2gm->use_gauss_legendre(scl2_gl);

    for (const auto& scl : scl2_vec) {
      std::string line;
      double sig;
      std::vector<double> sigma;
      fmt::print("# Enter XS[{}×μ₀,{}×μ₀] values separated by spaces: \n",
                 scl.first, scl.second);
      // line.clear();  // not needed
      std::getline(std::cin, line);
      // std::cout << "line = " << line << std::endl;
      std::istringstream stream(line);
      // std::cout << "stream = " << stream.str() << std::endl;
      std::copy(std::istream_iterator<double>(stream),
                std::istream_iterator<double>(), std::back_inserter(sigma));

      std::cout << "# Numbers you entered: ";
      std::copy(sigma.begin(), sigma.end(),
                std::ostream_iterator<double>(std::cout, " "));
      std::cout << '\n';

      scl2gm->add_model(scl, miho::GeometricModel(sigma));
    }

    cli_model = scl2gm;
  }

  //>>> evaluate cli_model <<<

  cli_model->set_max_nodes(nmax);
  cli_model->set_accuracy(accuracy);

  if (flag_pdf) {
    cli_model->print_nodes();
  }

  print_format(format_string, cli_model);



  return 0;

  double norm = cli_model->integrate([](double) { return 1.; });
  // fmt::print("# normalisation: {}", norm);

  double median2 = cli_model->median();
  std::pair<double, double> DoB682 = cli_model->degree_of_belief_interval();
  std::pair<double, double> DoB95 = cli_model->degree_of_belief_interval(0.95);
  double mean = cli_model->mean();
  double stdev = cli_model->stdev();
  double m1 = cli_model->moment(1);
  double m2 = cli_model->moment(2);
  double m3 = cli_model->moment(3);
  double variance =
      cli_model->integrate([=](double x) { return pow(x - mean, 2); });
  fmt::print("#result: {}  {}  {} {}  {} {}  {} {}  {} {} {} \n", mean, median2,
             DoB682.first, DoB682.second, DoB95.first, DoB95.second, stdev,
             std::sqrt(variance), m1, m2, m3);

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
  //   //   std::cout << val << "\t" << gm0.pdf(val) << "\t" << gm1.pdf(val) <<
  //   "\t"
  //   //             << gm2.pdf(val) << "\t" << gm3.pdf(val) << std::endl;
  //   // }
  //
  //   //double norm = gm.integrate([](double) { return 1.; });
  //   // double E = gm.integrate([](double x) { return x; });
  //   // double E2 = gm.integrate([](double x) { return x*x; });
  //
  //   return 0;
}


void print_format(std::string format_string, std::shared_ptr<miho::Model> model) {
  fmt::dynamic_format_arg_store<fmt::format_context> format_arg_list;
  if (std::regex_search(format_string, std::regex("\\{norm\\}"))) {
    double norm = model->norm();
    format_arg_list.push_back<double>(fmt::arg("norm", norm));
  }
  if (std::regex_search(format_string, std::regex("\\{median\\}"))) {
    double median = model->median();
    format_arg_list.push_back<double>(fmt::arg("median", median));
  }
  if (std::regex_search(format_string, std::regex("\\{dob68_low\\}")) ||
      std::regex_search(format_string, std::regex("\\{dob68_upp\\}"))) {
    std::pair<double, double> DoB68 = model->degree_of_belief_interval();
    format_arg_list.push_back(fmt::arg("dob68_low", DoB68.first));
    format_arg_list.push_back(fmt::arg("dob68_upp", DoB68.second));
  }
  if (std::regex_search(format_string, std::regex("\\{dob95_low\\}")) ||
      std::regex_search(format_string, std::regex("\\{dob95_upp\\}"))) {
    std::pair<double, double> DoB95 = model->degree_of_belief_interval(0.95);
    format_arg_list.push_back(fmt::arg("dob95_low", DoB95.first));
    format_arg_list.push_back(fmt::arg("dob95_upp", DoB95.second));
  }
  if (std::regex_search(format_string, std::regex("\\{mean\\}")) ||
      std::regex_search(format_string, std::regex("\\{stdev\\}"))) {
    double mean = model->mean();  // first compute the mean o recycle adaption
    double stdev = model->stdev();
    format_arg_list.push_back(fmt::arg("mean", mean));
    format_arg_list.push_back(fmt::arg("stdev", stdev));
  }
  fmt::vprint(format_string.c_str(), format_arg_list);
}

