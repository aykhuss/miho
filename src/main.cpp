#include <fmt/format.h>

#include <CLI/CLI.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

#include "ABCModel.h"
#include "GeometricModel.h"
#include "Model.h"
#include "Scale1DGeometricModel.h"
#include "Scale2DGeometricModel.h"

std::vector<double> parse_input(std::istream& data_stream);
std::vector<std::vector<double>> parse_data_file(const std::string& file_name);
void print_format(std::string, std::shared_ptr<miho::Model>);

int main(int argc, char const* argv[]) {
  CLI::App app{"ミホ (miho) - theory uncertainties from MIssing Higher Orders"};
  app.footer("[powered by CERN QCD coffee]");

  //----- global settings for all models
  bool flag_pdf = false;
  app.add_flag("--pdf", flag_pdf, "Print out the probability distribution.");
  std::string format_string{"{median} {dob68_low} {dob68_upp} \n"};
  CLI::Option* app_format =
      app.add_option("--format", format_string, "Output formatting string.");
  size_t nmax = 10000;
  app.add_option("--nmax", nmax, "Set the maximum number of PDF evaliations.");
  double accuracy = 0.001;  // default: 0.1%
  app.add_option("--accuracy", accuracy,
                 "Set the target accuracy of the integration.");
  std::string file_name;
  CLI::Option* app_file =
      app.add_option("--file,-f", file_name, "Provide an input data file.");

  //----- subcommands to choose the model
  app.require_subcommand(/* min */ 0, /* max */ 1);
  // standard geometric
  CLI::App* app_gm = app.add_subcommand("geometric", "Use the GeometricModel.");
  // ABC model
  CLI::App* app_abc = app.add_subcommand("abc", "Use the ABCModel.");
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

    std::vector<double> sigma;
    if (*app_file) {
      std::vector<std::vector<double>> data = parse_data_file(file_name);
      if (data.size()!=1) {
        throw std::runtime_error(file_name + " contains != 1 record(s)");
      }
      sigma = data.front();
    } else {
      fmt::print("# Enter XS values @ LO NLO ... separated by spaces: \n");
      sigma = parse_input(std::cin);

      // std::cout << "# Numbers you entered: ";
      // std::copy(sigma.begin(), sigma.end(),
      //           std::ostream_iterator<double>(std::cout, " "));
      // std::cout << '\n';
    }

    cli_model = std::shared_ptr<miho::Model>(new miho::GeometricModel(sigma));
  }

  //----- ABCModel
  if (app.got_subcommand(app_abc)) {
    fmt::print("# ABCModel\n");

    std::vector<double> sigma;
    if (*app_file) {
      std::vector<std::vector<double>> data = parse_data_file(file_name);
      if (data.size()!=1) {
        throw std::runtime_error(file_name + " contains != 1 record(s)");
      }
      sigma = data.front();
    } else {
      fmt::print("# Enter XS values @ LO NLO ... separated by spaces: \n");
      sigma = parse_input(std::cin);

      // std::cout << "# Numbers you entered: ";
      // std::copy(sigma.begin(), sigma.end(),
      //           std::ostream_iterator<double>(std::cout, " "));
      // std::cout << '\n';
    }

    cli_model = std::shared_ptr<miho::Model>(new miho::ABCModel(sigma));
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
            "# Enter XS[{}×μ₀,{}×μ₀] values @ LO NLO ... separated by spaces: \n",
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

  //>>> evaluate cli_model <<<

  cli_model->set_max_nodes(nmax);
  cli_model->set_accuracy(accuracy);

  if (flag_pdf) {
    cli_model->print_nodes();
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
    double stdev = model->stdev();
    format_arg_list.push_back(fmt::arg("mean", mean));
    format_arg_list.push_back(fmt::arg("stdev", stdev));
  }
  fmt::vprint(format_string + "\n", format_arg_list);
}
