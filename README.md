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

* probability distribution: get the x-y value pairs for the PDF by adding `--pdf` anywhere before the model specification
* use the format string passed via `--format 'the format string'` to control the information you want and how it is displayed, e.g. by default the format is `'{median} {dob68_low} {dob68_upp}'`. The list of available placeholders are:
    - `{median}`: the median of the distribution
    - `{dob68_low}` & `{dob68_upp}`: the lower and upper edges of the 68% degree-of-belief interval
    - `{dob95_low}` & `{dob95_upp}`: the lower and upper edges of the 95% degree-of-belief interval
    - `{mean}`: the mean of the distribution
    - `{stdev}`: the standard deviation of the distribution
    - `{norm}`: the normalisation of the PDF (should be compatible with `1`)

A step-by-step example is given [here](examples/gamma_N3LO/README.md) for the Drell-Yan cross section.

## Citation

[![arXiv](https://img.shields.io/badge/arXiv-2106.04585%20-orange)](https://arxiv.org/abs/2106.04585)
[![DOI](https://img.shields.io/badge/DOI-10.1007/JHEP09(2021)122-B31B1B.svg)](https://doi.org/10.1007/JHEP09(2021)122)

If you use miho, please cite the following reference:

```
@article{Duhr:2021mfd,
    author = "Duhr, Claude and Huss, Alexander and Mazeliauskas, Aleksas and Szafron, Robert",
    title = "{An analysis of Bayesian estimates for missing higher orders in perturbative calculations}",
    eprint = "2106.04585",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "CERN-TH-2021-058",
    doi = "10.1007/JHEP09(2021)122",
    journal = "JHEP",
    volume = "09",
    pages = "122",
    year = "2021"
}
```
