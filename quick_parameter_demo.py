#!/usr/bin/env python3
"""
Quick Parameter Modification Demo for URENIMOD Simulation

Author: Ryo Segawa (whizznihil.kid@gmail.com)

This script demonstrates how easy it is to modify parameters 
in the centralized parameter system.
"""

import numpy as np
import matplotlib.pyplot as plt
from membrane_pressure_simulation import UltrasonicMembraneSimulator

def demo_parameter_modification():
    """Demonstrate how to easily modify simulation parameters"""
    
    print("=" * 60)
    print("PARAMETER MODIFICATION DEMO")
    print("=" * 60)
    
    # ========== EASILY MODIFY THESE PARAMETERS ==========
    
    # Scenario 1: Default parameters (from PDF)
    scenario1_params = {
        'P_max': 50e3,       # 50 kPa
        't_centre': 50e-3,   # 50 ms
        'sigma_t': 10e-3,    # 10 ms
        'diameter': 1e-6,    # 1 μm
        'E': 1e6            # 1 MPa
    }
    
    # Scenario 2: Higher pressure, softer membrane
    scenario2_params = {
        'P_max': 150e3,      # 150 kPa (3x higher)
        't_centre': 50e-3,   # 50 ms
        'sigma_t': 10e-3,    # 10 ms
        'diameter': 1e-6,    # 1 μm
        'E': 0.5e6          # 0.5 MPa (softer)
    }
    
    # Scenario 3: Narrow pulse, larger axon
    scenario3_params = {
        'P_max': 100e3,      # 100 kPa
        't_centre': 30e-3,   # 30 ms (earlier)
        'sigma_t': 5e-3,     # 5 ms (narrower)
        'diameter': 5e-6,    # 5 μm (larger)
        'E': 1e6            # 1 MPa
    }
    
    # ====================================================
    
    scenarios = [
        ("Default (PDF)", scenario1_params),
        ("High Pressure + Soft", scenario2_params), 
        ("Narrow Pulse + Large Axon", scenario3_params)
    ]
    
    fig, axes = plt.subplots(3, 3, figsize=(15, 12))
    
    for i, (name, params) in enumerate(scenarios):
        print(f"\nRunning {name} scenario...")
        print(f"  Pressure: {params['P_max']/1000:.0f} kPa")
        print(f"  Pulse width: {params['sigma_t']*1000:.0f} ms")
        print(f"  Axon diameter: {params['diameter']*1e6:.0f} μm")
        print(f"  Young's modulus: {params['E']/1e6:.1f} MPa")
        
        # Create simulator with these parameters
        sim = UltrasonicMembraneSimulator(**params)
        results = sim.simulate(t_span=(0, 100e-3))
        
        if results['success']:
            t_ms = results['time'] * 1000
            
            # Plot pressure
            axes[i, 0].plot(t_ms, results['pressure']/1000, 'b-', linewidth=2)
            axes[i, 0].set_ylabel('Pressure (kPa)')
            axes[i, 0].set_title(f'{name}\nPressure Profile')
            axes[i, 0].grid(True, alpha=0.3)
            
            # Plot capacitive current
            axes[i, 1].plot(t_ms, results['capacitive_current']*1e6, 'r-', linewidth=2)
            axes[i, 1].set_ylabel('I_cap (μA/m²)')
            axes[i, 1].set_title('Capacitive Current')
            axes[i, 1].grid(True, alpha=0.3)
            
            # Plot membrane potential
            axes[i, 2].plot(t_ms, results['voltage']*1000, 'g-', linewidth=2)
            axes[i, 2].set_ylabel('Voltage (mV)')
            axes[i, 2].set_title('Membrane Potential')
            axes[i, 2].grid(True, alpha=0.3)
            
            if i == 2:  # Last row
                for j in range(3):
                    axes[i, j].set_xlabel('Time (ms)')
            
            # Print results
            V_max = np.max(results['voltage']) * 1000
            I_max = np.max(np.abs(results['capacitive_current'])) * 1e6
            print(f"  Max voltage: {V_max:.1f} mV")
            print(f"  Max capacitive current: {I_max:.1f} μA/m²")
        
        else:
            print(f"  Simulation failed!")
    
    plt.tight_layout()
    plt.savefig('figures/parameter_demo.png', dpi=300, bbox_inches='tight')
    print(f"\nDemo figure saved to figures/parameter_demo.png")
    print("\n" + "=" * 60)
    print("DEMO COMPLETE!")
    print("See how easy it is to modify parameters!")
    print("=" * 60)
    
    return fig

if __name__ == "__main__":
    demo_parameter_modification()
    plt.show() 