#include <CLI/CLI.hpp>
#include <filesystem>
#include <spdlog/fmt/bundled/ranges.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>

void shell_pipeline(std::string filename, std::string ser_file,
                    const std::map<std::string, bool> &controls,
                    double edge_bb_ratio, double, double, double);

void section_pipeline(std::string filename, std::string sec_file,
                      std::string ser_file, double target_edge = 0.1);

void transfer_pipeline(std::string shellfile, std::string section_file,
                       int upsample_level);

void boundary_perservation(std::string filename);

int main(int argc, char **argv) {
  CLI::App program{"Robustness Test for Bijective Remeshing"};
  program.require_subcommand(1);

  std::string filename, output_dir = "./", input_file;
  bool serial = false;
  program.add_option("-i,--input", input_file, "input name")
      ->required()
      ->check(CLI::ExistingFile)
      ->each([&filename](const std::string &s) {
        filename = std::filesystem::path(s).filename();
      });
  program.add_option("-o,--output", output_dir, "output dir")
      ->default_str("./");
  program.add_option_function<std::string>(
      "-l,--logdir",
      [&filename](const std::string &logdir) {
        auto file_logger = spdlog::basic_logger_mt(
            "vanilla", logdir + "/" + filename + ".log");
        spdlog::set_default_logger(file_logger);
        spdlog::flush_on(spdlog::level::info);
      },
      "log dir");
  program.add_option_function<int>(
      "--loglevel",
      [](const int &l) {
        spdlog::set_level(static_cast<spdlog::level::level_enum>(l));
      },
      "log level");
  program.add_flag("--serial", serial);

  ///// Shell construct
  auto shell = program.add_subcommand("shell");
  double edge_bb_ratio = 0.2;

  std::map<std::string, bool> pipeline_controls = {
      {"simplify-shell", false}, {"extract-stage", false},
      {"volume-centric", false}, {"allow-intersect", false},
      {"Parallel", true},        {"clean-shell", false},
  };

  auto enable_flag = [&shell, &pipeline_controls](std::string s,
                                                  std::string explain = "") {
    shell->add_flag("--" + s, pipeline_controls[s], explain);
  };
  enable_flag("simplify-shell");
  enable_flag("extract-stage");
  enable_flag("volume-centric");
  enable_flag("allow-intersect",
              "Allow the input **surface** to have self intersection.");
  enable_flag("clean-shell",
              "maintain a self-intersetion free throughout the process.");

  double angle_distortion = 1e-3;
  double target_thickness = 5e-3;
  double init_thick = 1e-4;
  std::string suffix = "";
  shell->add_option("--edge", edge_bb_ratio,
                    "target edge length / bounding box");
  shell->add_option("--angle", angle_distortion, "angle distortion (cosine)");
  shell->add_option("--thick", target_thickness,
                    "target thickness wrt bounding box");
  shell->add_option("--init-thick", init_thick,
                    "(Internal), thickness for initialization", true);
  shell->add_option("--suffix", suffix, "suffix identifier");
  shell->callback([&]() {
    filename = std::filesystem::path(input_file).filename();
    pipeline_controls["Parallel"] = !serial;
    init_thick = std::min(target_thickness / 2, init_thick);
    spdlog::info("{}", pipeline_controls);
    spdlog::info("angle {} target thick {} init thick {}", angle_distortion,
                 target_thickness, init_thick);
    shell_pipeline(input_file, output_dir + "/" + filename + suffix + ".h5",
                   pipeline_controls, edge_bb_ratio, angle_distortion,
                   target_thickness, init_thick);
  });

  ///// Section remesh
  auto section = program.add_subcommand("section");
  double edge = 0.3;
  std::string sec_file;
  section->add_option("--edge", edge, "target edge length / bounding box");
  section->add_option("--section-input", sec_file);
  section->callback([&]() {
    // double edge = 0.3;
    filename = std::filesystem::path(input_file).filename();
    section_pipeline(input_file, sec_file, output_dir + "/" + filename + ".ply",
                     edge);
  });

  ///// Transfer utils
  auto transfer = program.add_subcommand("transfer");
  int sample_level = 0;
  bool inverse = false;
  transfer->add_option("--section-file", sec_file);
  transfer->add_option("--sample-level", sample_level);
  transfer->callback(
      [&]() { transfer_pipeline(input_file, sec_file, sample_level); });

  ///// Experimental area
  auto experimental = program.add_subcommand("experimental");
  experimental->callback([&]() { boundary_perservation(input_file); });

  CLI11_PARSE(program, argc, argv);
}
