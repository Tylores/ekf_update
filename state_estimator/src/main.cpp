#include <interface/interface.hpp>
#include <iostream>

#ifdef DEBUG
#define LOG(message) std::cout << message << std::endl;
#else
#define LOG(message)
#endif

int main(int argc, char *argv[]) {
  std::cout << "State Estimator...\n";
  nlohmann::json topology = loadTopology();
  LOG(topology["admittance"])
  return 0;
}
