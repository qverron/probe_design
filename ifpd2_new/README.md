# iFISH Probe Design (II)

[![DOI](https://zenodo.org/badge/210782569.svg)](https://zenodo.org/badge/latestdoi/210782569) ![](https://img.shields.io/github/license/ggirelli/ifpd2.svg?style=flat) ![](https://github.com/ggirelli/ifpd2/workflows/Python%20package/badge.svg?branch=main&event=push)  
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/ifpd2) ![PyPI - Format](https://img.shields.io/pypi/format/ifpd2) ![PyPI - Status](https://img.shields.io/pypi/status/ifpd2)  
![](https://img.shields.io/github/release/ggirelli/ifpd2.svg?style=flat) ![](https://img.shields.io/github/release-date/ggirelli/ifpd2.svg?style=flat) ![](https://img.shields.io/github/languages/code-size/ggirelli/ifpd2.svg?style=flat)  
![](https://img.shields.io/github/watchers/ggirelli/ifpd2.svg?label=Watch&style=social) ![](https://img.shields.io/github/stars/ggirelli/ifpd2.svg?style=social)

[PyPi](https://pypi.org/project/ifpd2/) | [docs](https://ggirelli.github.io/ifpd2/)

`ifpd2` is a Python3.8+ package containing tools for selection of complementary oligonucleotides to build iFISH probes. It is based on our previous `ifpd` package, but works with a different and more detailed database format, allowing for more precise control on the probe design process. Read the online [documentation](https://ggirelli.github.io/ifpd2/) for more details.

## Requirements

`ifpd2` is fully implemented in Python3.8+, thus you need the corresponding Python version to run it. Check out [here](https://realpython.com/installing-python/) how to install Python+ on your machine if you don't have it yet.

`ifpd2` has been tested with Python 3.8 and 3.9. We recommend installing it using `pipx` (see [below](https://github.com/ggirelli/ifpd2#installation)) to avoid dependency conflicts with other packages. The packages it depends on are listed in our [dependency graph](https://github.com/ggirelli/ifpd2/network/dependencies). We use [`poetry`](https://github.com/python-poetry/poetry) to handle our dependencies.

## Installation

We recommend installing `ifpd2` using [`pipx`](https://github.com/pipxproject/pipx). Check how to install `pipx` [here](https://github.com/pipxproject/pipx#install-pipx) if you don't have it yet!

Once you have `pipx` ready on your system, install the latest stable release of `ifpd2` by running: `pipx install ifpd2`. If you see the stars (✨ 🌟 ✨), then the installation went well!

## Usage

All `ifpd2` commands are accessible via the `ifpd2` keyword on the terminal. For each command, you can access its help page by using the `-h` option. More details on how to run `ifpd2` are available in the online [documentation](https://ggirelli.github.io/ifpd2).

## Contributing

We welcome any contributions to `ifpd2`. In short, we use [`black`](https://github.com/psf/black) to standardize code format. Any code change also needs to pass `mypy` checks. For more details, please refer to our [contribution guidelines](https://github.com/ggirelli/ifpd2/blob/main/CONTRIBUTING.md) if this is your first time contributing! Also, check out our [code of conduct](https://github.com/ggirelli/ifpd2/blob/main/CODE_OF_CONDUCT.md).

## License

`MIT License - Copyright (c) 2021 Gabriele Girelli`
