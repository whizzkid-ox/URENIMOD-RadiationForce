"""
Author: Ryo Segawa (whizznihil.kid@gmail.com)

Test script for Ultrasonic Membrane Capacitance Simulation

This script implements the analysis and visualization specified in the PDF:
- Pressure-response curves
- Pulse-width-response curves  
- Mechanism visualization
- Parameter dependency analysis
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime
from membrane_pressure_simulation import UltrasonicMembraneSimulator

def test_basic_mechanism():
    """Test the basic pressure → capacitance → neural response mechanism"""
    print("Testing basic mechanism...")
    
    # Use default parameters for basic test
    sim = UltrasonicMembraneSimulator()
    
    # Test with different pressure levels
    pressures = [25e3, 50e3, 100e3, 150e3]  # 25, 50, 100, 150 kPa

def test_basic_mechanism_with_params(base_params, pressures, time_span):
    """Test the basic pressure → capacitance → neural response mechanism"""
    print("Testing basic mechanism...")
    
    # Create simulator with base parameters
    sim = UltrasonicMembraneSimulator(**base_params)
    
    fig, axes = plt.subplots(len(pressures), 4, figsize=(16, 12))
    
    for i, P_max in enumerate(pressures):
        sim.P_max = P_max
        results = sim.simulate(t_span=time_span)
        
        if results['success']:
            t_ms = results['time'] * 1000
            
            # Pressure profile
            axes[i, 0].plot(t_ms, results['pressure']/1000, 'b-', linewidth=2)
            axes[i, 0].set_ylabel(f'{P_max/1000:.0f} kPa\nPressure (kPa)')
            axes[i, 0].grid(True, alpha=0.3)
            if i == 0:
                axes[i, 0].set_title('Pressure P_US(t)')
            
            # Capacitive current
            axes[i, 1].plot(t_ms, results['capacitive_current']*1e6, 'r-', linewidth=2)
            axes[i, 1].set_ylabel('I_cap (μA/m²)')
            axes[i, 1].grid(True, alpha=0.3)
            if i == 0:
                axes[i, 1].set_title('Capacitive Current')
            
            # Membrane potential
            axes[i, 2].plot(t_ms, results['voltage']*1000, 'g-', linewidth=2)
            axes[i, 2].set_ylabel('V (mV)')
            axes[i, 2].grid(True, alpha=0.3)
            if i == 0:
                axes[i, 2].set_title('Membrane Potential')
            
            # Capacitance
            axes[i, 3].plot(t_ms, results['capacitance']*1e6, 'purple', linewidth=2)
            axes[i, 3].set_ylabel('C_m (μF/cm²)')
            axes[i, 3].grid(True, alpha=0.3)
            if i == 0:
                axes[i, 3].set_title('Membrane Capacitance')
            
            if i == len(pressures) - 1:
                for ax in axes[i, :]:
                    ax.set_xlabel('Time (ms)')
    
    plt.tight_layout()
    return fig

def test_pressure_response_curve():
    """Generate detailed pressure-response curve"""
    print("Generating pressure-response curve...")
    
    # Use default parameters for pressure response analysis
    sim = UltrasonicMembraneSimulator()
    
    # Wide range of pressures
    pressure_range = np.linspace(10e3, 300e3, 30)  # 10-300 kPa

def test_pressure_response_with_params(base_params, pressure_range_params, time_span):
    """Generate detailed pressure-response curve"""
    print("Generating pressure-response curve...")
    
    # Create simulator with base parameters
    sim = UltrasonicMembraneSimulator(**base_params)
    
    # Unpack pressure range parameters
    p_min, p_max, n_points = pressure_range_params
    pressure_range = np.linspace(p_min, p_max, n_points)
    responses = []
    max_voltages = []
    max_currents = []
    
    original_P_max = sim.P_max
    
    for P_max in pressure_range:
        sim.P_max = P_max
        results = sim.simulate(t_span=time_span)
        
        if results['success']:
            V_max = np.max(results['voltage'])
            I_cap_max = np.max(np.abs(results['capacitive_current']))
            
            max_voltages.append(V_max * 1000)  # Convert to mV
            max_currents.append(I_cap_max * 1e6)  # Convert to μA/m²
            
            # Check for action potential
            fired = V_max > -20e-3  # -20 mV threshold
            responses.append(1 if fired else 0)
        else:
            max_voltages.append(sim.V_rest * 1000)
            max_currents.append(0)
            responses.append(0)
    
    sim.P_max = original_P_max
    
    # Create comprehensive plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Response curve (binary)
    axes[0, 0].plot(pressure_range/1000, responses, 'bo-', linewidth=2, markersize=6)
    axes[0, 0].set_xlabel('Peak Pressure (kPa)')
    axes[0, 0].set_ylabel('Action Potential Fired (0/1)')
    axes[0, 0].set_title('Pressure-Response Curve (Binary)')
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].set_ylim(-0.1, 1.1)
    
    # Maximum voltage vs pressure
    axes[0, 1].plot(pressure_range/1000, max_voltages, 'g-', linewidth=2)
    axes[0, 1].axhline(y=-20, color='r', linestyle='--', label='AP Threshold')
    axes[0, 1].set_xlabel('Peak Pressure (kPa)')
    axes[0, 1].set_ylabel('Maximum Voltage (mV)')
    axes[0, 1].set_title('Peak Membrane Potential vs Pressure')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].legend()
    
    # Maximum capacitive current vs pressure
    axes[1, 0].plot(pressure_range/1000, max_currents, 'r-', linewidth=2)
    axes[1, 0].set_xlabel('Peak Pressure (kPa)')
    axes[1, 0].set_ylabel('Max Capacitive Current (μA/m²)')
    axes[1, 0].set_title('Peak Capacitive Current vs Pressure')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Current vs voltage (stimulus-response relationship)
    axes[1, 1].plot(max_currents, max_voltages, 'ko-', linewidth=2, markersize=4)
    axes[1, 1].axhline(y=-20, color='r', linestyle='--', label='AP Threshold')
    axes[1, 1].set_xlabel('Max Capacitive Current (μA/m²)')
    axes[1, 1].set_ylabel('Maximum Voltage (mV)')
    axes[1, 1].set_title('Voltage vs Current Relationship')
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].legend()
    
    plt.tight_layout()
    return fig, pressure_range, responses

def test_pulse_width_response():
    """Test pulse width effects"""
    print("Testing pulse width effects...")
    
    # Use default parameters for pulse width analysis
    sim = UltrasonicMembraneSimulator()
    
    # Range of pulse widths (σ_t)
    sigma_range = np.linspace(2e-3, 30e-3, 15)  # 2-30 ms
    pressure_levels = [75e3, 100e3, 125e3]  # Different pressure levels

def test_pulse_width_with_params(base_params, pressure_levels, pulse_width_range_params, time_span):
    """Test pulse width effects"""
    print("Testing pulse width effects...")
    
    # Create simulator with base parameters
    sim = UltrasonicMembraneSimulator(**base_params)
    
    # Unpack pulse width range parameters
    sigma_min, sigma_max, n_points = pulse_width_range_params
    sigma_range = np.linspace(sigma_min, sigma_max, n_points)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    for P_max in pressure_levels:
        responses = []
        max_dCdt = []
        
        sim.P_max = P_max
        original_sigma = sim.sigma_t
        
        for sigma_t in sigma_range:
            sim.sigma_t = sigma_t
            results = sim.simulate(t_span=time_span)
            
            if results['success']:
                V_max = np.max(results['voltage'])
                fired = V_max > -20e-3
                responses.append(1 if fired else 0)
                
                # Calculate max dC_m/dt
                dCdt_max = np.max(np.abs([sim.capacitance_derivative(t) for t in results['time']]))
                max_dCdt.append(dCdt_max * 1e9)  # Convert to nF/s per m²
            else:
                responses.append(0)
                max_dCdt.append(0)
        
        sim.sigma_t = original_sigma
        
        # Plot response vs pulse width
        axes[0].plot(sigma_range*1000, responses, 'o-', linewidth=2, 
                    label=f'{P_max/1000:.0f} kPa', markersize=6)
        
        # Plot max dC/dt vs pulse width
        axes[1].plot(sigma_range*1000, max_dCdt, 's-', linewidth=2,
                    label=f'{P_max/1000:.0f} kPa', markersize=6)
    
    axes[0].set_xlabel('Pulse Width σ_t (ms)')
    axes[0].set_ylabel('Action Potential Fired (0/1)')
    axes[0].set_title('Pulse-Width-Response Curves')
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()
    axes[0].set_ylim(-0.1, 1.1)
    
    axes[1].set_xlabel('Pulse Width σ_t (ms)')
    axes[1].set_ylabel('Max |dC_m/dt| (nF·s⁻¹·m⁻²)')
    axes[1].set_title('Peak Capacitance Derivative vs Pulse Width')
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()
    
    plt.tight_layout()
    return fig

def test_bipolar_mechanism():
    """Analyse the bipolar nature of the capacitive stimulus"""
    print("Analysing bipolar stimulus mechanism...")
    
    # Use higher pressure for clearer bipolar demonstration
    sim = UltrasonicMembraneSimulator(P_max=100e3)  # 100 kPa
    
    results = sim.simulate(t_span=(0, 100e-3))

def test_bipolar_with_params(base_params, bipolar_pressure, time_span):
    """Analyse the bipolar nature of the capacitive stimulus"""
    print("Analysing bipolar stimulus mechanism...")
    
    # Create simulator with base parameters and set pressure for bipolar demo
    params = base_params.copy()
    params['P_max'] = bipolar_pressure
    sim = UltrasonicMembraneSimulator(**params)
    
    results = sim.simulate(t_span=time_span)
    
    if not results['success']:
        return None
        
    t_ms = results['time'] * 1000
    
    # Calculate components
    dCdt = np.array([sim.capacitance_derivative(t) for t in results['time']])
    I_cap = results['capacitive_current']
    
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    
    # Pressure and its derivative (to show bipolar nature)
    axes[0].plot(t_ms, results['pressure']/1000, 'b-', linewidth=2, label='P_US(t)')
    ax0_twin = axes[0].twinx()
    dPdt = np.gradient(results['pressure'], results['time'])
    ax0_twin.plot(t_ms, dPdt/1e6, 'b--', linewidth=2, label='dP/dt', alpha=0.7)
    axes[0].set_ylabel('Pressure (kPa)', color='b')
    ax0_twin.set_ylabel('dP/dt (MPa/s)', color='b')
    axes[0].set_title('Pressure Profile and Derivative')
    axes[0].grid(True, alpha=0.3)
    axes[0].axvline(x=sim.t_centre*1000, color='k', linestyle=':', alpha=0.5, label='Peak time')
    axes[0].legend(loc='upper left')
    ax0_twin.legend(loc='upper right')
    
    # Capacitance derivative (shows bipolar current source)
    axes[1].plot(t_ms, dCdt*1e9, 'r-', linewidth=2, label='dC_m/dt')
    axes[1].axhline(y=0, color='k', linestyle='-', alpha=0.3)
    axes[1].axvline(x=sim.t_centre*1000, color='k', linestyle=':', alpha=0.5)
    axes[1].set_ylabel('dC_m/dt (nF·s⁻¹·m⁻²)')
    axes[1].set_title('Capacitance Derivative (Bipolar Nature)')
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()
    
    # Resulting capacitive current and membrane potential
    axes[2].plot(t_ms, I_cap*1e6, 'r-', linewidth=2, label='I_cap = V × dC_m/dt')
    ax2_twin = axes[2].twinx()
    ax2_twin.plot(t_ms, results['voltage']*1000, 'g-', linewidth=2, label='V(t)')
    axes[2].axhline(y=0, color='k', linestyle='-', alpha=0.3)
    axes[2].axvline(x=sim.t_centre*1000, color='k', linestyle=':', alpha=0.5)
    axes[2].set_ylabel('Capacitive Current (μA/m²)', color='r')
    ax2_twin.set_ylabel('Membrane Potential (mV)', color='g')
    axes[2].set_xlabel('Time (ms)')
    axes[2].set_title('Capacitive Current and Membrane Response')
    axes[2].grid(True, alpha=0.3)
    axes[2].legend(loc='upper left')
    ax2_twin.legend(loc='upper right')
    
    plt.tight_layout()
    return fig

def main():
    """Main test function"""
    print("=" * 70)
    print("COMPREHENSIVE MEMBRANE CAPACITANCE SIMULATION ANALYSIS")
    print("=" * 70)
    
    # ========== MAIN PARAMETERS FOR ALL TESTS - MODIFY HERE ==========
    
    # Base Simulation Parameters (used across all tests)
    P_max_base = 50e3        # Pa - Base peak pressure for standard tests
    t_centre = 50e-3         # s - Pulse center time (50 ms)
    sigma_t_base = 10e-3     # s - Base pulse width (10 ms)
    
    # Axon Parameters
    length = 3e-3            # m - Axon length (3 mm)
    diameter = 1e-6          # m - Axon diameter (1 μm)
    d0 = 1.4e-9             # m - Static membrane thickness (1.4 nm)
    E = 1e6                 # Pa - Young's modulus (1 MPa)
    
    # Hodgkin-Huxley Parameters
    C_m0 = 1e-2             # F/m² - Static membrane capacitance
    g_Na = 1.2              # S/m² - Sodium conductance
    g_K = 0.36              # S/m² - Potassium conductance
    g_L = 0.003             # S/m² - Leak conductance
    E_Na = 50e-3            # V - Sodium reversal potential
    E_K = -77e-3            # V - Potassium reversal potential
    E_L = -54.4e-3          # V - Leak reversal potential
    V_rest = -65e-3         # V - Resting potential
    
    # Test-Specific Parameters
    mechanism_pressures = [0, 50e3, 100e3, 500e3, 1000e3]    # kPa for mechanism test
    pressure_response_range = (0, 1000e3, 30)         # min, max, points for pressure curve
    pulse_width_pressures = [0, 50e3, 100e3, 500e3, 1000e3]        # kPa for pulse width test
    pulse_width_range = (0, 30e-3, 15)               # min, max, points for sigma test
    bipolar_pressure = 1000e3                             # kPa for bipolar demonstration
    
    # Simulation Time Spans
    basic_time_span = (0, 100e-3)      # 100 ms for basic tests
    extended_time_span = (0, 150e-3)   # 150 ms for response curves
    
    # ================================================================
    
    print(f"Test Parameters:")
    print(f"  Base pressure: {P_max_base/1000:.1f} kPa")
    print(f"  Pulse center: {t_centre*1000:.1f} ms")
    print(f"  Base pulse width: {sigma_t_base*1000:.1f} ms")
    print(f"  Pressure range: {pressure_response_range[0]/1000:.0f}-{pressure_response_range[1]/1000:.0f} kPa")
    print(f"  Pulse width range: {pulse_width_range[0]*1000:.0f}-{pulse_width_range[1]*1000:.0f} ms")
    print("")
    
    # Create base parameter set for creating simulators
    base_params = {
        'P_max': P_max_base, 't_centre': t_centre, 'sigma_t': sigma_t_base,
        'length': length, 'diameter': diameter, 'd0': d0, 'E': E,
        'C_m0': C_m0, 'g_Na': g_Na, 'g_K': g_K, 'g_L': g_L,
        'E_Na': E_Na, 'E_K': E_K, 'E_L': E_L, 'V_rest': V_rest
    }
    
    # Create figures directory
    os.makedirs('figures', exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Test 1: Basic mechanism visualisation
    try:
        fig1 = test_basic_mechanism_with_params(base_params, mechanism_pressures, basic_time_span)
        path1 = os.path.join('figures', f'mechanism_analysis_{timestamp}.png')
        fig1.savefig(path1, dpi=300, bbox_inches='tight')
        print(f"✓ Mechanism analysis saved to {path1}")
    except Exception as e:
        print(f"✗ Mechanism analysis failed: {e}")
    
    # Test 2: Pressure-response curve (detailed)
    try:
        fig2, pressure_range, responses = test_pressure_response_with_params(
            base_params, pressure_response_range, extended_time_span)
        path2 = os.path.join('figures', f'pressure_response_detailed_{timestamp}.png')
        fig2.savefig(path2, dpi=300, bbox_inches='tight')
        print(f"✓ Detailed pressure-response analysis saved to {path2}")
        
        # Find threshold
        firing_indices = np.where(np.array(responses) == 1)[0]
        if len(firing_indices) > 0:
            threshold_pressure = pressure_range[firing_indices[0]]
            print(f"  → Firing threshold: {threshold_pressure/1000:.1f} kPa")
        else:
            print("  → No firing observed in tested range")
            
    except Exception as e:
        print(f"✗ Pressure-response analysis failed: {e}")
    
    # Test 3: Pulse-width-response curves
    try:
        fig3 = test_pulse_width_with_params(
            base_params, pulse_width_pressures, pulse_width_range, extended_time_span)
        path3 = os.path.join('figures', f'pulse_width_analysis_{timestamp}.png')
        fig3.savefig(path3, dpi=300, bbox_inches='tight')
        print(f"✓ Pulse-width analysis saved to {path3}")
    except Exception as e:
        print(f"✗ Pulse-width analysis failed: {e}")
    
    # Test 4: Bipolar mechanism analysis
    try:
        fig4 = test_bipolar_with_params(base_params, bipolar_pressure, basic_time_span)
        if fig4:
            path4 = os.path.join('figures', f'bipolar_mechanism_{timestamp}.png')
            fig4.savefig(path4, dpi=300, bbox_inches='tight')
            print(f"✓ Bipolar mechanism analysis saved to {path4}")
    except Exception as e:
        print(f"✗ Bipolar mechanism analysis failed: {e}")
    
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE!")
    print("All figures saved in the 'figures' directory.")
    print("=" * 70)
    
    plt.show()

if __name__ == "__main__":
    main() 