# URENIMOD Radiation Force Simulation

**Ultrasonic Membrane Capacitance Simulation for Neural Stimulation**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## Overview

This repository implements the ultrasonic membrane capacitance simulation described in the URENIMOD radiation force proposal. The simulation models how focused ultrasound pressure waves can modulate neural membrane capacitance to achieve non-invasive neural stimulation.

### Key Features

- **Physics-Based Model**: Implements pressure → membrane deformation → capacitance → neural response pathway
- **Parameterized Design**: Easily modify all simulation parameters in centralized locations
- **Comprehensive Analysis**: Includes pressure-response curves, pulse-width analysis, and bipolar mechanism studies
- **Publication-Ready Figures**: Generates high-quality plots for research publications

## Quick Start

### Prerequisites

```bash
pip install numpy scipy matplotlib
```

### Basic Simulation

```python
python membrane_pressure_simulation.py
```

### Comprehensive Analysis

```python
python test_membrane_simulation.py
```

### Parameter Demo

```python
python quick_parameter_demo.py
```

## Core Physics Model

The simulation implements the following key equations from the proposal:

1. **Gaussian Pressure Profile**: P_US(t) = P_max × exp(-(t - t_centre)²/(2×σ_t²))
2. **Membrane Deformation**: d(t) = d₀ × (1 - P_US(t)/E)
3. **Time-Varying Capacitance**: C_m(t) = C_m0 / (1 - P_US(t)/E)
4. **Modified Hodgkin-Huxley**: dV/dt = (1/C_m(t)) × (-I_ion - V×dC_m/dt + I_ext)
5. **Capacitive Current**: I_cap(t) = V(t) × dC_m/dt

## Easy Parameter Modification

All parameters are centralized in the `main()` functions for easy experimentation:

```python
# ========== MAIN PARAMETERS - MODIFY HERE FOR EXPERIMENTS ==========

# Ultrasonic Pressure Parameters
P_max = 50e3         # Pa - Peak pressure (main variable to explore)
t_centre = 50e-3     # s - Pulse center time (50 ms)
sigma_t = 10e-3      # s - Pulse width standard deviation (10 ms)

# Axon Parameters
diameter = 1e-6      # m - Axon diameter (1 μm)
E = 1e6             # Pa - Young's modulus (1 MPa)

# ... other parameters ...
```

## File Structure

```
URENIMOD-RadiationForce/
├── membrane_pressure_simulation.py    # Main simulator with centralized params
├── test_membrane_simulation.py        # Comprehensive analysis suite
├── quick_parameter_demo.py            # Interactive parameter demo
├── implementation_summary.md          # Detailed technical documentation
├── requirements.txt                   # Python dependencies
├── figures/                           # Generated analysis plots
└── README.md                          # This file
```

## Results

The simulation generates multiple analysis figures:

- **Basic Mechanism**: Pressure → Current → Voltage response chains
- **Pressure-Response Curves**: Threshold analysis across pressure ranges
- **Pulse-Width Analysis**: Effect of temporal pulse characteristics
- **Bipolar Mechanism**: Demonstration of depolarizing/hyperpolarizing phases

## Research Applications

This simulation is designed for:

- **Parameter Optimization**: Finding optimal ultrasonic parameters for neural stimulation
- **Safety Analysis**: Studying pressure thresholds and membrane effects
- **Pulse Design**: Optimizing temporal characteristics for selective stimulation
- **Mechanism Understanding**: Analyzing the capacitive current pathway

## Citation

If you use this simulation in your research, please cite:

```
Segawa, R. (2025). URENIMOD Ultrasonic Membrane Capacitance Simulation. 
GitHub repository: https://github.com/whizzkid-ox/URENIMOD-RadiationForce
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

**Ryo Segawa** (whizznihil.kid@gmail.com)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Acknowledgments

- Based on the URENIMOD radiation force simulation proposal
- Implements modified Hodgkin-Huxley neural dynamics
- Uses numerical integration for accurate capacitive current modeling 