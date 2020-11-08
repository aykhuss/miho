#include <fmt/color.h>
#include <fmt/format.h>

#include <CLI/CLI.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <regex>
#include <string>
#include <vector>

#include "ABCModel.h"
#include "ABCNumericModel.h"
#include "GeometricModel.h"
#include "Model.h"
#include "Scale1DGeometricModel.h"
#include "Scale1DModel.h"
#include "Scale2DGeometricModel.h"
#include "Scale2DModel.h"

// fwd declare some functions
std::vector<double> parse_input(std::istream& data_stream);
std::vector<std::vector<double>> parse_data_file(const std::string& file_name);
void print_format(std::string, std::shared_ptr<miho::Model>);

int main(int argc, char const* argv[]) {
  CLI::App app{"ミホ (miho) - theory uncertainties from MIssing Higher Orders"};
  app.footer("[powered by CERN QCD coffee]");

  //----- global settings for all models
  //> output format
  std::string format_string{"{median} {dob68_low} {dob68_upp} \n"};
  CLI::Option* app_format =
      app.add_option("--format", format_string, "Output formatting string.");
  //> PDF settings
  bool flag_pdf = false;
  app.add_flag("--pdf", flag_pdf, "Print out the probability distribution.");
  size_t nmax = 2000;
  app.add_option("--nmax", nmax,
                 "Set the maximum number of PDF evaluations (default: 5000).");
  double accuracy = 0.005;  // default: 0.5%
  app.add_option(
      "--accuracy", accuracy,
      "Set the relative target accuracy of the integration (default: 0.5%).");
  //> scale marginalisation settings
  std::vector<double> vec_scl1d;
  CLI::Option* app_scl1d =
      app.add_option("--scl1d", vec_scl1d, "A list of scale factors.");
  std::vector<std::pair<double, double>> vec_scl2d;
  CLI::Option* app_scl2d =
      app.add_option("--scl2d", vec_scl2d, "A list of pairs of scale factors.");
  bool scl_gl = false;
  app.add_flag("--gl", scl_gl, "Gauss–Legendre");
  // app_scl_gl->multi_option_policy(CLI::MultiOptionPolicy::Throw);
  app_scl1d->excludes(app_scl2d);
  app_scl2d->excludes(app_scl1d);
  //> file as input
  std::string file_name;
  CLI::Option* app_file =
      app.add_option("--file,-f", file_name, "Provide an input data file.");

  //----- subcommands to choose the model
  app.require_subcommand(/* min */ 1, /* max */ 1);
  //> standard geometric
  CLI::App* app_gm = app.add_subcommand("geo", "Use the GeometricModel.");
  int gm_omg = 1;
  app_gm->add_option("--omega", gm_omg,
                     "Model parameter for the prior of 'a' (default: 1).");
  double gm_eps = 0.1;
  app_gm->add_option("--epsilon", gm_eps,
                     "Model parameter for the prior of 'c' (default: 0.1).");
  //> ABC model
  CLI::App* app_abc = app.add_subcommand("abc", "Use the ABCModel.");
  int abc_omg = 1;
  app_abc->add_option("--omega", abc_omg,
                      "Model parameter for the prior of 'a' (default: 1).");
  double abc_xi = 1.;
  app_abc->add_option("--xi", abc_xi,
                      "Model parameter for the prior of 'b' (default: 1.).");
  double abc_eps = 0.1;
  app_abc->add_option("--epsilon", abc_eps,
                      "Model parameter for the prior of 'c' (default: 0.1).");
  // 1D scale
  CLI::App* app_scl1gm = app.add_subcommand(
      "scale1d", "Use the GeometricModel marginalising over one scale.");
  std::vector<double> scl1_vec;
  app_scl1gm->add_option("--scales", scl1_vec, "The scale factors.")
      ->excludes(app_file);
  bool scl1_gl = false;
  app_scl1gm->add_flag("--gauss-legendre", scl1_gl,
                       "Use the 3-point Gauss–Legendre quadrature.");
  // 2D scale
  CLI::App* app_scl2gm = app.add_subcommand(
      "scale2d", "Use the GeometricModel marginalising over two scales.");
  std::vector<std::pair<double, double>> scl2_vec;
  app_scl2gm->add_option("--scales", scl2_vec, "The scale factors as pairs.")
      ->excludes(app_file);
  bool scl2_gl = false;
  app_scl2gm->add_flag("--gauss-legendre", scl2_gl,
                       "Use the 3-point Gauss–Legendre quadrature.");

  CLI11_PARSE(app, argc, argv);

  std::shared_ptr<miho::Model> cli_model = nullptr;

  //----- GeometricModel
  if (app.got_subcommand(app_gm)) {
    fmt::print("# GeometricModel\n");
    //> create a bare Model object with the correct model settings
    miho::GeometricModel gm;
    gm.set_omega(gm_omg);
    gm.set_epsilon(gm_eps);

    if (*app_scl1d) {
      //### 1D scale
      // fmt::print("entered app_scl1d\n");
      std::shared_ptr<miho::Scale1DModel> scl1d =
          std::shared_ptr<miho::Scale1DModel>(new miho::Scale1DModel());
      scl1d->use_gauss_legendre(scl_gl);
      if (*app_file) {
        std::vector<std::vector<double>> data = parse_data_file(file_name);
        for (std::vector<double> record : data) {
          double scl = record.front();
          std::vector<double> sigma(record.begin() + 1, record.end());
          gm.set_sigma(sigma);
          scl1d->add_model(scl, std::make_shared<miho::GeometricModel>(gm));
        }
      } else {
        for (const auto& scl : vec_scl1d) {
          fmt::print(
              "# Enter XS[{}×μ₀] values @ LO NLO ... separated by spaces: \n",
              scl);
          std::vector<double> sigma = parse_input(std::cin);
          gm.set_sigma(sigma);
          scl1d->add_model(scl, std::make_shared<miho::GeometricModel>(gm));
        }
      }
      cli_model = scl1d;
    } else if (*app_scl2d) {
      //### 2D scale
      // fmt::print("entered app_scl2d\n");
      std::shared_ptr<miho::Scale2DModel> scl2d =
          std::shared_ptr<miho::Scale2DModel>(new miho::Scale2DModel());
      scl2d->use_gauss_legendre(scl_gl);
      if (*app_file) {
        std::vector<std::vector<double>> data = parse_data_file(file_name);
        for (std::vector<double> record : data) {
          std::pair<double, double> scl{record.at(0), record.at(1)};
          std::vector<double> sigma(record.begin() + 2, record.end());
          gm.set_sigma(sigma);
          scl2d->add_model(scl, std::make_shared<miho::GeometricModel>(gm));
        }
      } else {
        for (const auto& scl : vec_scl2d) {
          fmt::print(
              "# Enter XS[{}×μ₀,{}×μ₀] values @ LO NLO ... separated by "
              "spaces: "
              "\n",
              scl.first, scl.second);
          std::vector<double> sigma = parse_input(std::cin);
          gm.set_sigma(sigma);
          scl2d->add_model(scl, std::make_shared<miho::GeometricModel>(gm));
        }
      }
      cli_model = scl2d;
    } else {
      //### just numbers
      std::vector<double> sigma;
      if (*app_file) {
        std::vector<std::vector<double>> data = parse_data_file(file_name);
        if (data.size() != 1) {
          throw std::runtime_error(file_name + " contains != 1 record(s)");
        }
        sigma = data.front();
      } else {
        fmt::print("# Enter XS values @ LO NLO ... separated by spaces: \n");
        sigma = parse_input(std::cin);
      }
      gm.set_sigma(sigma);
      cli_model = std::make_shared<miho::GeometricModel>(gm);
    }
  }

  //----- ABCModel
  if (app.got_subcommand(app_abc)) {
    fmt::print("# ABCNumericModel\n");
    //> create a bare Model object with the correct model settings
    miho::ABCNumericModel abc;
    abc.set_omega(abc_omg);
    abc.set_xi(abc_xi);
    abc.set_epsilon(abc_eps);

    if (*app_scl1d) {
      //### 1D scale
      // fmt::print("entered app_scl1d\n");
      std::shared_ptr<miho::Scale1DModel> scl1d =
          std::shared_ptr<miho::Scale1DModel>(new miho::Scale1DModel());
      scl1d->use_gauss_legendre(scl_gl);
      if (*app_file) {
        std::vector<std::vector<double>> data = parse_data_file(file_name);
        for (std::vector<double> record : data) {
          double scl = record.front();
          std::vector<double> sigma(record.begin() + 1, record.end());
          abc.set_sigma(sigma);
          scl1d->add_model(scl, std::make_shared<miho::ABCNumericModel>(abc));
        }
      } else {
        for (const auto& scl : vec_scl1d) {
          fmt::print(
              "# Enter XS[{}×μ₀] values @ LO NLO ... separated by spaces: \n",
              scl);
          std::vector<double> sigma = parse_input(std::cin);
          abc.set_sigma(sigma);
          scl1d->add_model(scl, std::make_shared<miho::ABCNumericModel>(abc));
        }
      }
      cli_model = scl1d;
    } else if (*app_scl2d) {
      //### 2D scale
      // fmt::print("entered app_scl2d\n");
      std::shared_ptr<miho::Scale2DModel> scl2d =
          std::shared_ptr<miho::Scale2DModel>(new miho::Scale2DModel());
      scl2d->use_gauss_legendre(scl_gl);
      if (*app_file) {
        std::vector<std::vector<double>> data = parse_data_file(file_name);
        for (std::vector<double> record : data) {
          std::pair<double, double> scl{record.at(0), record.at(1)};
          std::vector<double> sigma(record.begin() + 2, record.end());
          abc.set_sigma(sigma);
          scl2d->add_model(scl, std::make_shared<miho::ABCNumericModel>(abc));
        }
      } else {
        for (const auto& scl : vec_scl2d) {
          fmt::print(
              "# Enter XS[{}×μ₀,{}×μ₀] values @ LO NLO ... separated by "
              "spaces: "
              "\n",
              scl.first, scl.second);
          std::vector<double> sigma = parse_input(std::cin);
          abc.set_sigma(sigma);
          scl2d->add_model(scl, std::make_shared<miho::ABCNumericModel>(abc));
        }
      }
      cli_model = scl2d;
    } else {
      //### just numbers
      std::vector<double> sigma;
      if (*app_file) {
        std::vector<std::vector<double>> data = parse_data_file(file_name);
        if (data.size() != 1) {
          throw std::runtime_error(file_name + " contains != 1 record(s)");
        }
        sigma = data.front();
      } else {
        fmt::print("# Enter XS values @ LO NLO ... separated by spaces: \n");
        sigma = parse_input(std::cin);
      }
      abc.set_sigma(sigma);
      cli_model = std::make_shared<miho::ABCNumericModel>(abc);
    }
  }

  //----- Scale1DGeometricModel
  if (app.got_subcommand(app_scl1gm)) {
    fmt::print("# Scale1DGeometricModel\n");

    // std::cout << "# Scales you entered: ";
    // std::copy(scl1_vec.begin(), scl1_vec.end(),
    //           std::ostream_iterator<double>(std::cout, " "));
    // std::cout << '\n';

    std::shared_ptr<miho::Scale1DGeometricModel> scl1gm =
        std::shared_ptr<miho::Scale1DGeometricModel>(
            new miho::Scale1DGeometricModel());

    scl1gm->use_gauss_legendre(scl1_gl);

    if (*app_file) {
      std::vector<std::vector<double>> data = parse_data_file(file_name);
      for (std::vector<double> record : data) {
        double scl = record.front();
        std::vector<double> sigma(record.begin() + 1, record.end());
        scl1gm->add_model(scl, miho::GeometricModel(sigma));
      }
    } else {
      for (const auto& scl : scl1_vec) {
        fmt::print(
            "# Enter XS[{}×μ₀] values @ LO NLO ... separated by spaces: \n",
            scl);
        std::vector<double> sigma = parse_input(std::cin);
        // std::cout << "# Numbers you entered: ";
        // std::copy(sigma.begin(), sigma.end(),
        //           std::ostream_iterator<double>(std::cout, " "));
        // std::cout << '\n';
        scl1gm->add_model(scl, miho::GeometricModel(sigma));
      }
    }

    cli_model = scl1gm;
  }

  //----- Scale2DGeometricModel
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

    if (*app_file) {
      std::vector<std::vector<double>> data = parse_data_file(file_name);
      for (std::vector<double> record : data) {
        std::pair<double, double> scl{record.at(0), record.at(1)};
        std::vector<double> sigma(record.begin() + 2, record.end());
        scl2gm->add_model(scl, miho::GeometricModel(sigma));
      }
    } else {
      for (const auto& scl : scl2_vec) {
        fmt::print(
            "# Enter XS[{}×μ₀,{}×μ₀] values @ LO NLO ... separated by spaces: "
            "\n",
            scl.first, scl.second);
        std::vector<double> sigma = parse_input(std::cin);
        // std::cout << "# Numbers you entered: ";
        // std::copy(sigma.begin(), sigma.end(),
        //           std::ostream_iterator<double>(std::cout, " "));
        // std::cout << '\n';
        scl2gm->add_model(scl, miho::GeometricModel(sigma));
      }
    }

    cli_model = scl2gm;
  }

  //> check that cli_model is set, otherwise scream for help & die
  if (cli_model == nullptr) {
    fmt::print(fg(fmt::color::crimson) | fmt::emphasis::bold,
               "\nno model selected!\n\n");
    std::cout << app.help() << std::endl;
    return 1;
  }

  //>>> evaluate cli_model <<<

  cli_model->set_max_nodes(nmax);
  cli_model->set_accuracy(accuracy);

  if (flag_pdf) {
    cli_model->print_pdf();
    if (*app_format) print_format(format_string, cli_model);
  } else {
    print_format(format_string, cli_model);
  }

  return 0;
}

std::vector<double> parse_input(std::istream& data_stream) {
  std::vector<double> record_vec;
  std::string line;
  std::getline(std::cin, line);
  std::istringstream input_stream(line);
  std::copy(std::istream_iterator<double>(input_stream),
            std::istream_iterator<double>(), std::back_inserter(record_vec));
  return record_vec;
}

std::vector<std::vector<double>> parse_data_file(const std::string& file_name) {
  // fmt::print("You gave an input file {}\n", file_name);
  std::ifstream data_file(file_name, std::ios_base::in);

  if (data_file.is_open()) {
    std::vector<std::vector<double>> data;

    std::string record;
    while (std::getline(data_file, record)) {
      if (std::regex_match(record, std::regex("^[[:space:]]*#.*")) ||
          std::regex_match(record, std::regex("^[[:space:]]*$")))
        continue;
      std::istringstream record_stream(record);
      std::vector<double> record_vec;
      double record_entry;
      while (record_stream >> record_entry) record_vec.push_back(record_entry);
      data.emplace_back(record_vec);
    }
    data_file.close();
    // fmt::print("read vector:\n");
    // for (const std::vector<double>& line : data) {
    //   for (const double& entry : line) fmt::print("{} ", entry);
    //   fmt::print("\n");
    // }
    return data;
  } else {
    throw std::runtime_error("couldn't open data file " + file_name);
  }
}

void print_format(std::string format_string,
                  std::shared_ptr<miho::Model> model) {
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
    format_arg_list.push_back(fmt::arg("mean", mean));
  }
  if (std::regex_search(format_string, std::regex("\\{stdev\\}"))) {
    double stdev = model->stdev();
    format_arg_list.push_back(fmt::arg("stdev", stdev));
  }
  //> maybe cache mean to save time for stdev?
  // if (std::regex_search(format_string, std::regex("\\{mean\\}")) ||
  //     std::regex_search(format_string, std::regex("\\{stdev\\}"))) {
  //   double mean = model->mean();  // first compute the mean o recycle adaption
  //   double stdev = model->stdev();
  //   format_arg_list.push_back(fmt::arg("mean", mean));
  //   format_arg_list.push_back(fmt::arg("stdev", stdev));
  // }
  fmt::vprint(format_string + "\n", format_arg_list);
}
