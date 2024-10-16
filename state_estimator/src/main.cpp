#include <ekf/ekf.hpp>
#include <format>
#include <interface/interface.hpp>
#include <iostream>
#include <ranges>
#include <set>
#include <stdio.h>
#include <string>
#include <unordered_map>
#include <vector>

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
using CVEC = std::vector<std::complex<double>>;
using CMAT = std::vector<CVEC>;

struct Info {
  json topology;
  SIMAP index_lookup;
  ISMAP name_lookup;
  SCMAP base_voltages;
  SCMAP base_powers;
  SCMAP voltages;
  SCMAP powers;
};

CMAT extractAdmittance(const json &topology, const SIMAP &lookup) {
  json obj = topology["admittance"];
  json from = obj["from_equipment"];
  json to = obj["to_equipment"];
  json values = obj["admittance_list"];
  std::set<std::string> nodes(from);

  CMAT mat(nodes.size(), CVEC(nodes.size()));
  for (const auto [src, dst, value] : std::views::zip(from, to, values)) {
    int i, j = 0;
    if (lookup.contains(src))
      i = lookup.at(src);
    if (lookup.contains(dst))
      j = lookup.at(dst);
    if (i != 0 and j != 0)
      mat[i][j] = std::complex<double>(value[0], value[1]);
  }
  return mat;
};

SDMAP extractRealInjections(const json &injections) {
  json obj = injections["power_real"];
  json ids = obj["ids"];
  json values = obj["values"];
  json units = obj["units"];

  double scale = 1.0;
  if (units == "kW")
    scale = 1000.0;

  SDMAP mapping;
  for (const auto [id, value] : std::views::zip(ids, values)) {
    double x = scale * (double)value;
    if (mapping.contains(id)) {
      mapping[id] += x;
    } else {
      mapping.insert({id, x});
    }
  }
  return mapping;
};

SDMAP extractImaginaryInjections(const json &injections) {
  json obj = injections["power_imaginary"];
  json ids = obj["ids"];
  json values = obj["values"];
  json units = obj["units"];

  double scale = 1.0;
  if (units == "kVAR")
    scale = 1000.0;

  SDMAP mapping;
  for (const auto [id, value] : std::views::zip(ids, values)) {
    double x = scale * (double)value;
    if (mapping.contains(id)) {
      mapping[id] += x;
    } else {
      mapping.insert({id, x});
    }
  }
  return mapping;
};

SDMAP extractBaseVoltagesMagnitudes(const json &topology) {
  json obj = topology["base_voltage_magnitudes"];
  json ids = obj["ids"];
  json values = obj["values"];

  SDMAP mapping;
  for (const auto [id, value] : std::views::zip(ids, values)) {
    mapping.insert({id, (double)value});
  }
  return mapping;
};

SDMAP extractBaseVoltagesAngles(const json &topology) {
  json obj = topology["base_voltage_angles"];
  json ids = obj["ids"];
  json values = obj["values"];

  SDMAP mapping;
  for (const auto [id, value] : std::views::zip(ids, values)) {
    mapping.insert({id, (double)value});
  }
  return mapping;
};

SIMAP extractNodes(const json &topology) {
  json obj = topology["base_voltage_magnitudes"];
  json ids = obj["ids"];

  SIMAP mapping;
  for (const auto [idx, id] : std::views::enumerate(ids)) {
    mapping.insert({id, idx});
  }
  return mapping;
};

ISMAP reflectMap(const SIMAP &lookup) {
  ISMAP mapping;
  for (const auto &[key, value] : lookup) {
    mapping.insert({value, key});
  }
  return mapping;
}

SCMAP cartMap(const SDMAP &real, const SDMAP &imag) {
  SCMAP mapping;
  for (const auto &e : real) {
    std::string key = e.first;
    mapping.insert({key, std::complex<double>(real.at(key), imag.at(key))});
  }
  return mapping;
}

SCMAP polarMap(const SDMAP &mag, const SDMAP &ang) {
  SCMAP mapping;
  for (const auto &e : mag) {
    std::string key = e.first;
    std::complex<double> value = std::polar(mag.at(key), ang.at(key));
    mapping.insert({key, value});
  }
  return mapping;
}

void printSCMap(const SCMAP &map) {
  for (const auto &e : map) {
    LOG(std::format("{} = {},{}", e.first, e.second.real(), e.second.imag()));
  }
}

int main(int argc, char *argv[]) {
  std::cout << "State Estimator...\n";
  json topology = loadTopology();
  SIMAP index_lookup = extractNodes(topology);
  ISMAP name_lookup = reflectMap(index_lookup);
  LOG(std::format("lookup: {}", index_lookup.size()));
  SDMAP base_voltage_mag = extractBaseVoltagesMagnitudes(topology);
  SDMAP base_voltage_ang = extractBaseVoltagesAngles(topology);
  SCMAP base_voltages = polarMap(base_voltage_mag, base_voltage_ang);
  LOG(std::format("base_voltages: {}", base_voltages.size()));
  SDMAP base_power_real = extractRealInjections(topology["injections"]);
  SDMAP base_power_imag = extractImaginaryInjections(topology["injections"]);
  SCMAP base_powers = cartMap(base_power_real, base_power_imag);
  LOG(std::format("powers: {}", base_powers.size()));
  CMAT admittance = extractAdmittance(topology, index_lookup);
  LOG(std::format("admittance: {}, {}", admittance.size(),
                  admittance[0].size()));
  return 0;
}
