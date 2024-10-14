# ekf_update
Trying to rewrite ekf with more modern c++

## Setup

```shell
sudo apt update
sudo apt upgrade -y
sudo apt install -y liblapack-dev libblas-dev libmpfr-dev libgmp-dev
```

## Build

```shell
cmake -S . -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=1
cmake --build build
```

## Run

```shell
./build/bin/state_estimator
```
