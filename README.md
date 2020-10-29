# ミホ (miho) 
> theory uncertainties from **MI**ssing **H**igher **O**rders

## Clone the repository

External libraries are included in this repository through submodules. In addition to cloning the repository, you therefore also need to initialise and update those modules:
```bash
git clone git@github.com:aykhuss/miho.git
cd miho
git submodule init
git submodule update
```

or just:
```bash
git clone --recurse-submodules git@github.com:aykhuss/miho.git
```

## Install

We use [CMake](https://cmake.org/) to build our project. 
On macOS using [homebrew](https://brew.sh/), install it like so: `brew install cmake`.

To build and install the project, run:
```bash
mkdir build
cd build
cmake ..
make
make install
```
that should create a `bin` folder containing the `miho` executable in the main project folder. 


## Usage

Run `miho -h` to get a help display. 

maybe use logging at some point:
* https://github.com/gabime/spdlog
* https://github.com/SergiusTheBest/plog
