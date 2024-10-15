#include <ekf/ekf.hpp>
#include <format>
#include <interface/interface.hpp>
#include <iostream>
#include <ranges>
#include <stdio.h>
#include <string>
#include <unordered_map>

#ifdef DEBUG
#define LOG(message) std::cout << message << std::endl;
#else
#define LOG(message)
#endif

std::string in_dir = DATA_INPUT_DIR;
std::string out_dir = DATA_OUTPUT_DIR;

using json = nlohmann::json;
using SIMAP = std::unordered_map<std::string, size_t>;
using ISMAP = std::unordered_map<size_t, std::string>;
using SCMAP = std::unordered_map<std::string, std::complex<double>>;
using ICMAP = std::unordered_map<size_t, std::complex<double>>;
using SDMAP = std::unordered_map<std::string, double>;
using IDMAP = std::unordered_map<size_t, double>;

struct Info {
  SIMAP index_lookup;
  ISMAP name_lookup;
  SDMAP base_voltages;
  SCMAP base_powers;
  SCMAP voltages;
  SCMAP powers;
};

SDMAP extractRealInjections(const json &topology) {
  json obj = topology["injections"]["power_real"];
  json ids = obj["ids"];
  json values = obj["values"];

  SDMAP mapping;
  for (const auto [id, value] : std::views::zip(ids, values)) {
    mapping.insert({id, value});
  }
  return mapping;
};

SDMAP extractImagInjections(const json &topology) {
  json obj = topology["injections"]["power_imaginary"];
  json ids = obj["ids"];
  json values = obj["values"];

  SDMAP mapping;
  for (const auto [id, value] : std::views::zip(ids, values)) {
    mapping.insert({id, value});
  }
  return mapping;
};

SDMAP extractBaseVoltages(const json &topology) {
  json obj = topology["base_voltage_magnitudes"];
  json ids = obj["ids"];
  json values = obj["values"];

  SDMAP mapping;
  for (const auto [id, value] : std::views::zip(ids, values)) {
    mapping.insert({id, value});
  }
  return mapping;
};

SIMAP extractNodes(const json &topology) {
  json obj = topology["base_voltage_magnitudes"];
  json ids = obj["ids"];

  SIMAP mapping;
  for (const auto [idx, id] : std::views::enumerate(ids)) {
    mapping.insert({id, idx});
    LOG(std::format("{} : {}", (std::string)id, idx));
  }
  return mapping;
};

ISMAP reflect(const SIMAP &lookup) {
  ISMAP mapping;
  for (const auto &[key, value] : lookup) {
    mapping.insert({value, key});
  }
  return mapping;
}

int main(int argc, char *argv[]) {
  std::cout << "State Estimator...\n";
  json topology = loadTopology();
  SIMAP index_lookup = extractNodes(topology);
  ISMAP name_lookup = reflect(index_lookup);
  return 0;
}
