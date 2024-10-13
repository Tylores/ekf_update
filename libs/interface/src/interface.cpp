#include <fstream>
#include <interface/interface.hpp>

using namespace nlohmann;

std::string in_dir = DATA_INPUT_DIR;
std::string out_dir = DATA_OUTPUT_DIR;

json loadTopology() {
  std::ifstream f(in_dir + "/topology.json");
  return json::parse(f);
}
