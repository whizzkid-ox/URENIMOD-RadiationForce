# URENIMOD Ultrasonic Membrane Capacitance Simulation

**Author**: Ryo Segawa (whizznihil.kid@gmail.com)  

## Quick Start for Parameter Modification

All simulation parameters are now defined in the `main()` function for easy experimentation. Simply edit the parameters in the marked section:

```python
def main():
    # ========== MAIN PARAMETERS - MODIFY HERE FOR EXPERIMENTS ==========
    
    # Ultrasonic Pressure Parameters
    P_max = 50e3         # Pa - Peak pressure (main variable to explore)
    t_centre = 50e-3     # s - Pulse center time (50 ms)
    sigma_t = 10e-3      # s - Pulse width standard deviation (10 ms)
    
    # Axon Parameters
    length = 3e-3        # m - Axon length (3 mm)
    diameter = 1e-6      # m - Axon diameter (1 μm)
    d0 = 1.4e-9         # m - Static membrane thickness (1.4 nm)
    E = 1e6             # Pa - Young's modulus (1 MPa)
    
    # Hodgkin-Huxley Parameters  
    C_m0 = 1e-2         # F/m² - Static membrane capacitance (1 μF/cm²)
    g_Na = 1.2          # S/m² - Sodium conductance (120 mS/cm²)
    g_K = 0.36          # S/m² - Potassium conductance (36 mS/cm²)
    g_L = 0.003         # S/m² - Leak conductance (0.3 mS/cm²)
    
    # Nernst Potentials
    E_Na = 50e-3        # V - Sodium reversal potential (50 mV)
    E_K = -77e-3        # V - Potassium reversal potential (-77 mV)
    E_L = -54.4e-3      # V - Leak reversal potential (-54.4 mV)
    V_rest = -65e-3     # V - Resting potential (-65 mV)
    
    # ====================================================================
```

### Common Parameter Modifications:

1. **Higher pressure for stimulation**:
   ```python
   P_max = 150e3  # 150 kPa instead of 50 kPa
   ```

2. **Different pulse characteristics**:
   ```python
   t_centre = 25e-3   # Earlier pulse at 25 ms
   sigma_t = 5e-3     # Narrower pulse (5 ms width)
   ```

3. **Different axon properties**:
   ```python
   diameter = 10e-6   # Larger axon (10 μm)
   E = 0.5e6         # Softer membrane (0.5 MPa)
   ```

### 🎯 **Core Physics Model**

**Pressure → Membrane Deformation → Capacitance Change → Neural Response**

1. **Gaussian Pressure Profile**:
   ```
   P_US(t) = P_max × exp(-(t - t_centre)²/(2×σ_t²))
   ```

2. **Membrane Thickness Change**:
   ```
   d(t) = d₀ × (1 - P_US(t)/E)
   ```

3. **Time-Varying Capacitance**:
   ```
   C_m(t) = C_m0 / (1 - P_US(t)/E)
   ```

4. **Modified Hodgkin-Huxley Equation**:
   ```
   dV/dt = (1/C_m(t)) × (-I_ion - V×dC_m/dt + I_ext)
   ```

5. **Capacitive Current** (the key mechanism):
   ```
   I_cap(t) = V(t) × dC_m/dt
   ```

### 📊 **Implemented Analysis (Per PDF Specifications)**

#### **Step 4.1: Key Time-Series Plotting**
✅ **Causal Chain Visualisation**:
- Input pressure profile P_US(t)
- Resulting capacitive current I_cap(t) = V(t) × dC_m/dt  
- Membrane potential response V(t)
- Time-varying capacitance C_m(t)

#### **Step 4.2: Parameter Dependencies**
✅ **Pressure-Response Curve**:
- Systematic variation of P_max (10-300 kPa)
- Binary firing response (yes/no)
- Threshold pressure determination
- Peak voltage vs pressure analysis

✅ **Pulse-Width-Response Curve**:
- Variation of σ_t (pulse width) 
- Effect on dC_m/dt magnitude
- Relationship between pulse steepness and threshold

#### **Step 4.3: Mechanism Interpretation**
✅ **Bipolar Stimulus Analysis**:
- Rising phase (t < t_centre): depolarizing current
- Falling phase (t > t_centre): hyperpolarizing current  
- Potential anode-break excitation patterns

### 🗂️ **Generated Files**

#### **Core Implementation**
- `membrane_pressure_simulation.py` - Main simulator class
- `test_membrane_simulation.py` - Comprehensive analysis suite

#### **Generated Analysis Figures**
- `mechanism_analysis_*.png` - Multi-pressure mechanism visualization
- `pressure_response_detailed_*.png` - Comprehensive pressure-response analysis
- `pulse_width_analysis_*.png` - Pulse width dependency curves
- `bipolar_mechanism_*.png` - Bipolar stimulus mechanism analysis
- `membrane_simulation_*.png` - Individual simulation results

### 📋 **Key Parameters (From PDF)**

#### **Ultrasonic Pressure**
- Peak Pressure (P_max): 50 kPa (adjustable main variable)
- Pulse Centre Time (t_centre): 50 ms
- Pulse Width (σ_t): 10 ms

#### **Axon Properties**
- Length: 3 mm
- Diameter: 1 μm
- Membrane thickness (d₀): 1.4 nm
- Young's modulus (E): 1 MPa

#### **Hodgkin-Huxley Parameters**
- Static capacitance (C_m0): 1 μF/cm²
- Standard HH conductances and potentials
- Resting potential: -65 mV

### 🔬 **Analysis Results**

#### **Mechanism Confirmation**
✅ **Smooth pressure profile successfully elicits membrane response via capacitive mechanism**

#### **Bipolar Nature**
✅ **Rising phase**: Creates depolarizing current (positive dC_m/dt)
✅ **Falling phase**: Creates hyperpolarizing current (negative dC_m/dt)

#### **Parameter Dependencies**
- **Firing threshold**: Dependent on peak pressure
- **Pulse width effect**: Smaller σ_t → larger dC_m/dt → lower pressure threshold
- **Current magnitude**: Proportional to V(t) × dC_m/dt

### 🚀 **Usage Examples**

#### **Basic Simulation**
```python
from membrane_pressure_simulation import UltrasonicMembraneSimulator

sim = UltrasonicMembraneSimulator()
results = sim.simulate(t_span=(0, 100e-3))
fig = sim.plot_results(results, save_path='figures/simulation.png')
```

#### **Pressure-Response Analysis**
```python
pressure_range, responses, fig = sim.pressure_response_curve()
```

#### **Comprehensive Analysis**
```bash
python test_membrane_simulation.py
```

### 📈 **Key Findings**

1. **Capacitive current mechanism works**: V(t) × dC_m/dt provides effective neural stimulation
2. **Bipolar nature confirmed**: Both depolarizing and hyperpolarizing phases present
3. **Parameter sensitivity**: Pulse width significantly affects threshold pressure
4. **Smooth pressure profiles**: More physiologically realistic than square waves

## Comprehensive Testing with Custom Parameters

The test file (`test_membrane_simulation.py`) also has centralized parameters for all analysis tests:

```python
def main():
    # ========== MAIN PARAMETERS FOR ALL TESTS - MODIFY HERE ==========
    
    # Base Simulation Parameters
    P_max_base = 50e3        # Pa - Base pressure for standard tests
    t_centre = 50e-3         # s - Pulse center time
    sigma_t_base = 10e-3     # s - Base pulse width
    
    # Test-Specific Parameters
    mechanism_pressures = [25e3, 50e3, 100e3, 150e3]    # kPa for mechanism test
    pressure_response_range = (10e3, 300e3, 30)         # min, max, points
    pulse_width_pressures = [75e3, 100e3, 125e3]        # kPa for pulse width test
    pulse_width_range = (2e-3, 30e-3, 15)               # min, max, points
    bipolar_pressure = 100e3                             # kPa for bipolar demo
    
    # Simulation Time Spans
    basic_time_span = (0, 100e-3)      # 100 ms for basic tests
    extended_time_span = (0, 150e-3)   # 150 ms for response curves
    
    # ================================================================
```

### Example Test Modifications:

1. **Test higher pressure ranges**:
   ```python
   pressure_response_range = (50e3, 500e3, 40)  # 50-500 kPa, 40 points
   mechanism_pressures = [100e3, 200e3, 300e3, 400e3]  # Higher pressures
   ```

2. **Different pulse width analysis**:
   ```python
   pulse_width_range = (1e-3, 50e-3, 25)  # 1-50 ms, 25 points
   pulse_width_pressures = [150e3, 200e3, 250e3]  # Higher test pressures
   ```

3. **Extended time analysis**:
   ```python
   extended_time_span = (0, 200e-3)  # 200 ms simulation time
   ```
