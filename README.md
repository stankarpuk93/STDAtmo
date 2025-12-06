# STDAtmo
A lightweight desktop calculator of 1976 U.S. Standard Atmosphere 


# Standard Atmosphere Calculator

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Qt](https://img.shields.io/badge/Qt-6.10-green.svg)](https://www.qt.io/)
[![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey.svg)]()

A lightweight Qt-based calculator for standard atmospheric properties with an US ISA 1976 implementation.

## Features

- **Complete ISA 1976 Model** - 8 atmospheric layers (0-85 km)
- **Multi-Unit Support** - SI and US Customary units
- **Temperature Deviation** - Calculate with non-standard conditions
- **Export Results** - Save to CSV format
- **Real-Time Calculations** - Instant parameter updates

## Installation

### For End Users
Download the latest release for your platform from the [Releases](../../releases) page.

### For Developers

#### Prerequisites
- Qt Creator 10+ with Qt 6.10+
- MinGW 11.2+ compiler (included with Qt installer)

#### Build Steps
1. **Open Qt Creator**
2. **File → Open File or Project**
3. Select `CMakeLists.txt`
4. Click **Configure Project**
5. Click the hammer icon **Build** (Ctrl+B)
6. Click the play button **Run** (Ctrl+R)

## Usage

1. Select unit system (SI or US)
2. Enter geometric altitude
3. Set temperature deviation (optional)
4. Click **Calculate**

### Output Parameters
- Temperature (K/°C or °R/°F)
- Pressure (Pa/hPa or psf/psi)
- Density (kg/m³ or slug/ft³)
- Speed of Sound (m/s or ft/s)
- Dynamic Viscosity (Pa·s or lb·s/ft²)

## Project Structure