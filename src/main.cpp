#include <fmt/args.h>
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
#include "ScaleModel.h"
#include "Scale1DModel.h"
#include "Scale2DModel.h"
#include "Util.h"

// fwd declare some functions
std::vector<double> parse_input(std::istream& data_stream);
std::vector<std::vector<double>> parse_data_file(const std::string& file_name);
void print_format(std::string, std::unique_ptr<miho::Model>&);

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
  size_t nmax = 5000;
  app.add_option("--nmax", nmax,
                 "Set the maximum number of PDF evaluations (default: 5000).");
  double accuracy = 0.005;  // default: 0.5%
  app.add_option(
      "--accuracy", accuracy,
      "Set the relative target accuracy of the integration (default: 0.5%).");
  //> scale marginalisation settings
  std::size_t sclNd;
  CLI::Option* app_sclNd =
      app.add_option("--sclNd", sclNd, "The # of independent scales.");
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
  CLI::App* app_geo = app.add_subcommand("geo", "Use the GeometricModel.");
  int geo_omg = 1;
  app_geo->add_option("--omega", geo_omg,
                      "Model parameter for the prior of 'a' (default: 1).");
  double geo_eps = 0.1;
  app_geo->add_option("--epsilon", geo_eps,
                      "Model parameter for the prior of 'c' (default: 0.1).");
  //> ABC model
  CLI::App* app_abc = app.add_subcommand("abc", "Use the ABCModel.");
  int abc_omg = 1;
  app_abc->add_option("--omega", abc_omg,
                      "Model parameter for the prior of 'a' (default: 1).");
  double abc_xi = 2.;
  app_abc->add_option("--xi", abc_xi,
                      "Model parameter for the prior of 'b' (default: 2.).");
  double abc_eps = 0.1;
  app_abc->add_option("--epsilon", abc_eps,
                      "Model parameter for the prior of 'c' (default: 0.1).");
  double abc_eta = 0.1;
  app_abc->add_option("--eta", abc_eta,
                      "Model parameter for the prior of 'c' (default: 0.1).");

  CLI11_PARSE(app, argc, argv);

  /// set the underlying prototype model with the correct settings
  std::unique_ptr<miho::ModelPrototype> proto_model = nullptr;
  if (app.got_subcommand(app_geo)) {
    fmt::print("# GeometricModel\n");
    auto gm = std::make_unique<miho::GeometricModel>();
    gm->set_omega(geo_omg);
    gm->set_epsilon(geo_eps);
    proto_model = std::move(gm);
  } else if (app.got_subcommand(app_abc)) {
    fmt::print("# ABCModel\n");
    auto abc = std::make_unique<miho::ABCModel>();
    abc->set_omega(abc_omg);
    abc->set_xi(abc_xi);
    abc->set_epsilon(abc_eps);
    abc->set_eta(abc_eta);
    proto_model = std::move(abc);
  } else {
    fmt::print(fg(fmt::color::crimson) | fmt::emphasis::bold,
               "\nno model selected!\n\n");
    std::cout << app.help() << std::endl;
    return 1;
  }

  /// set the final model (possible scale marginalisations, etc...)
  std::unique_ptr<miho::Model> cli_model = nullptr;

if (*app_sclNd) {
    //### ND scale
    // fmt::print("entered app_sclNd\n");
    auto sclNd = std::make_unique<miho::ScaleModel<2>>();
    sclNd->use_gauss_legendre(scl_gl);
    if (*app_file) {
      std::vector<std::vector<double>> data = parse_data_file(file_name);
      for (std::vector<double> record : data) {
        std::array<double, 2> scl({record.at(0), record.at(1)});
        std::vector<double> sigma(record.begin() + 2, record.end());
        auto clone_model = proto_model->clone();
        clone_model->set_sigma(sigma);
        sclNd->add_model(scl, std::move(clone_model));
      }
    }
    cli_model = std::move(sclNd);

  } else if (*app_scl1d) {
    //### 1D scale
    // fmt::print("entered app_scl1d\n");
    auto scl1d = std::make_unique<miho::Scale1DModel>();
    scl1d->use_gauss_legendre(scl_gl);
    if (*app_file) {
      std::vector<std::vector<double>> data = parse_data_file(file_name);
      for (std::vector<double> record : data) {
        double scl = record.front();
        std::vector<double> sigma(record.begin() + 1, record.end());
        auto clone_model = proto_model->clone();
        clone_model->set_sigma(sigma);
        scl1d->add_model(scl, std::move(clone_model));
      }
    } else {
      for (const auto& scl : vec_scl1d) {
        fmt::print(
            "# Enter XS[{}×μ₀] values @ LO NLO ... separated by spaces: \n",
            scl);
        std::vector<double> sigma = parse_input(std::cin);
        auto clone_model = proto_model->clone();
        clone_model->set_sigma(sigma);
        scl1d->add_model(scl, std::move(clone_model));
      }
    }
    cli_model = std::move(scl1d);
  } else if (*app_scl2d) {
    //### 2D scale
    // fmt::print("entered app_scl2d\n");
    auto scl2d = std::make_unique<miho::Scale2DModel>();
    scl2d->use_gauss_legendre(scl_gl);
    if (*app_file) {
      std::vector<std::vector<double>> data = parse_data_file(file_name);
      for (std::vector<double> record : data) {
        std::pair<double, double> scl{record.at(0), record.at(1)};
        std::vector<double> sigma(record.begin() + 2, record.end());
        auto clone_model = proto_model->clone();
        clone_model->set_sigma(sigma);
        scl2d->add_model(scl, std::move(clone_model));
      }
    } else {
      for (const auto& scl : vec_scl2d) {
        fmt::print(
            "# Enter XS[{}×μ₀,{}×μ₀] values @ LO NLO ... separated by "
            "spaces: "
            "\n",
            scl.first, scl.second);
        std::vector<double> sigma = parse_input(std::cin);
        auto clone_model = proto_model->clone();
        clone_model->set_sigma(sigma);
        scl2d->add_model(scl, std::move(clone_model));
      }
    }
    cli_model = std::move(scl2d);
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
    proto_model->set_sigma(sigma);
    cli_model = std::move(proto_model);
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

  /// run the adaption and get the 68% DoB interval
  /// override the target accuracy if distribution is too narrow
  std::pair<double, double> dob68 = cli_model->degree_of_belief_interval();
  double rel_dob =
      std::fabs((dob68.first - dob68.second) / (dob68.first + dob68.second));
  if (rel_dob < accuracy) {
    std::cerr << "overriding target accuracy to: " << rel_dob
              << "[input: " << accuracy << "]" << std::endl;
    // cli_model->clear();
    cli_model->set_accuracy(rel_dob);
  }

  if (flag_pdf) {
    // fmt::print("\nprinting PDF...\n");
    cli_model->print_pdf();
    if (*app_format) print_format(format_string, cli_model);
  } else {
    // fmt::print("\nprinting format...\n");
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
                  std::unique_ptr<miho::Model>& model) {
  fmt::dynamic_format_arg_store<fmt::format_context> format_arg_list;
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
  if (std::regex_search(format_string, std::regex("\\{norm\\}"))) {
    double norm = model->norm();
    format_arg_list.push_back<double>(fmt::arg("norm", norm));
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
  //   double mean = model->mean();  // first compute the mean o recycle
  //   adaption double stdev = model->stdev();
  //   format_arg_list.push_back(fmt::arg("mean", mean));
  //   format_arg_list.push_back(fmt::arg("stdev", stdev));
  // }
  fmt::vprint(format_string + "\n", format_arg_list);
}
