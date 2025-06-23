"""
Author: Ryo Segawa (whizznihil.kid@gmail.com)


Ultrasonic Pressure-Induced Membrane Capacitance Simulation

This module implements the simulation plan from URENIMOD_simulation_RadiationForce_v1.pdf:
- Gaussian pressure pulse → membrane deformation → capacitance change → modified HH model

The core mechanism: P(t) → Δd(t) → C_m(t) → modified dV/dt with capacitive current
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from neuron import h
import logging
import os
from datetime import datetime

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class UltrasonicMembraneSimulator:
    """
    Simulator for ultrasonic pressure effects on axon membrane capacitance
    
    Implements the pressure → membrane deformation → capacitance → neural response model
    """
    
    def __init__(self, P_max=50e3, t_centre=50e-3, sigma_t=10e-3, 
                 length=3e-3, diameter=1e-6, d0=1.4e-9, E=1e6,
                 C_m0=1e-2, g_Na=1.2, g_K=0.36, g_L=0.003,
                 E_Na=50e-3, E_K=-77e-3, E_L=-54.4e-3, V_rest=-65e-3):
        """
        Initialize simulator with parameters from the PDF simulation plan
        
        Parameters:
        -----------
        Ultrasonic Pressure Parameters:
        - P_max : float, Peak pressure (Pa) - main variable to explore
        - t_centre : float, Pulse center time (s)
        - sigma_t : float, Pulse width standard deviation (s)
        
        Axon Parameters:
        - length : float, Axon length (m)
        - diameter : float, Axon diameter (m)
        - d0 : float, Static membrane thickness (m)
        - E : float, Young's modulus (Pa)
        
        Hodgkin-Huxley Parameters:
        - C_m0 : float, Static membrane capacitance (F/m²)
        - g_Na, g_K, g_L : float, Ionic conductances (S/m²)
        - E_Na, E_K, E_L : float, Nernst potentials (V)
        - V_rest : float, Resting potential (V)
        """
        
        # Ultrasonic Pressure Parameters
        self.P_max = P_max        
        self.t_centre = t_centre   
        self.sigma_t = sigma_t     
        
        # Axon Parameters (Single Compartment)
        self.length = length       
        self.diameter = diameter   
        self.area = np.pi * self.diameter * self.length  # m² - Membrane area
        
        # Physical Properties
        self.d0 = d0              
        self.E = E                
        
        # Hodgkin-Huxley Parameters
        self.C_m0 = C_m0          
        self.g_Na = g_Na          
        self.g_K = g_K            
        self.g_L = g_L            
        
        # Nernst potentials (V)
        self.E_Na = E_Na          
        self.E_K = E_K            
        self.E_L = E_L            
        
        self.V_rest = V_rest      
        
        logger.info("UltrasonicMembraneSimulator initialised")
        logger.info(f"Peak pressure: {self.P_max/1000:.1f} kPa")
        logger.info(f"Pulse center: {self.t_centre*1000:.1f} ms")
        logger.info(f"Pulse width: {self.sigma_t*1000:.1f} ms")
    
    def pressure_profile(self, t):
        """Calculate temporal pressure profile P_US(t) = P_max * exp(-(t - t_centre)²/(2*σ_t²))"""
        return self.P_max * np.exp(-((t - self.t_centre)**2) / (2 * self.sigma_t**2))
    
    def membrane_thickness(self, t):
        """Calculate time-varying membrane thickness d(t) = d0 * (1 - P_US(t) / E)"""
        P_t = self.pressure_profile(t)
        return self.d0 * (1 - P_t / self.E)
    
    def membrane_capacitance(self, t):
        """Calculate time-varying membrane capacitance C_m(t) = C_m0 / (1 - P_US(t) / E)"""
        P_t = self.pressure_profile(t)
        return self.C_m0 / (1 - P_t / self.E)
    
    def capacitance_derivative(self, t):
        """Calculate dC_m/dt analytically"""
        P_t = self.pressure_profile(t)
        dP_dt = P_t * (-(t - self.t_centre) / self.sigma_t**2)
        denominator = (1 - P_t / self.E)**2
        return (self.C_m0 / self.E) * dP_dt / denominator
    
    def alpha_m(self, V):
        """Sodium activation rate constant"""
        V_mv = V * 1000  # Convert to mV for HH equations
        if abs(V_mv + 40) < 1e-6:
            return 1.0
        return 0.1 * (V_mv + 40) / (1 - np.exp(-(V_mv + 40) / 10))
    
    def beta_m(self, V):
        """Sodium inactivation rate constant"""
        V_mv = V * 1000
        return 4.0 * np.exp(-(V_mv + 65) / 18)
    
    def alpha_h(self, V):
        """Sodium inactivation rate constant"""
        V_mv = V * 1000
        return 0.07 * np.exp(-(V_mv + 65) / 20)
    
    def beta_h(self, V):
        """Sodium inactivation rate constant"""
        V_mv = V * 1000
        return 1.0 / (1 + np.exp(-(V_mv + 35) / 10))
    
    def alpha_n(self, V):
        """Potassium activation rate constant"""
        V_mv = V * 1000
        if abs(V_mv + 55) < 1e-6:
            return 0.1
        return 0.01 * (V_mv + 55) / (1 - np.exp(-(V_mv + 55) / 10))
    
    def beta_n(self, V):
        """Potassium inactivation rate constant"""
        V_mv = V * 1000
        return 0.125 * np.exp(-(V_mv + 65) / 80)
    
    def ionic_currents(self, V, m, h, n):
        """Calculate total ionic current I_ion"""
        I_Na = self.g_Na * m**3 * h * (V - self.E_Na)
        I_K = self.g_K * n**4 * (V - self.E_K)  
        I_L = self.g_L * (V - self.E_L)
        return I_Na + I_K + I_L
    
    def capacitive_current(self, V, t):
        """Calculate capacitive current I_cap = V * dC_m/dt"""
        return V * self.capacitance_derivative(t)
    
    def hh_system(self, t, y, I_ext=0):
        """Modified Hodgkin-Huxley system with time-varying capacitance"""
        V, m, h, n = y
        
        # Gating variable derivatives
        dm_dt = self.alpha_m(V) * (1 - m) - self.beta_m(V) * m
        dh_dt = self.alpha_h(V) * (1 - h) - self.beta_h(V) * h  
        dn_dt = self.alpha_n(V) * (1 - n) - self.beta_n(V) * n
        
        # Ionic current
        I_ion = self.ionic_currents(V, m, h, n)
        
        # Capacitive current (the key ultrasonic stimulus)
        I_cap = self.capacitive_current(V, t)
        
        # Time-varying capacitance
        C_m_t = self.membrane_capacitance(t)
        
        # Modified membrane potential equation: dV/dt = (1/C_m(t)) * (-I_ion - V*dC_m/dt + I_ext)
        dV_dt = (1 / C_m_t) * (-I_ion - I_cap + I_ext)
        
        return [dV_dt, dm_dt, dh_dt, dn_dt]
    
    def steady_state_values(self, V):
        """Calculate steady-state gating variable values"""
        m_inf = self.alpha_m(V) / (self.alpha_m(V) + self.beta_m(V))
        h_inf = self.alpha_h(V) / (self.alpha_h(V) + self.beta_h(V))
        n_inf = self.alpha_n(V) / (self.alpha_n(V) + self.beta_n(V))
        return m_inf, h_inf, n_inf
    
    def simulate(self, t_span=(0, 100e-3), dt=1e-5, I_ext=0):
        """Run the complete pressure-capacitance-neural simulation"""
        
        # Time vector
        t_eval = np.arange(t_span[0], t_span[1], dt)
        
        # Initial conditions (resting state)
        V0 = self.V_rest
        m0, h0, n0 = self.steady_state_values(V0)
        y0 = [V0, m0, h0, n0]
        
        # Solve the modified HH system
        logger.info(f"Simulating from {t_span[0]*1000:.1f} to {t_span[1]*1000:.1f} ms...")
        
        sol = solve_ivp(
            lambda t, y: self.hh_system(t, y, I_ext),
            t_span, y0, t_eval=t_eval,
            method='RK45', rtol=1e-8, atol=1e-10
        )
        
        if not sol.success:
            logger.error(f"Integration failed: {sol.message}")
            return None
        
        # Extract solution
        t = sol.t
        V, m, h, n = sol.y
        
        # Calculate derived quantities
        P_t = np.array([self.pressure_profile(ti) for ti in t])
        C_m_t = np.array([self.membrane_capacitance(ti) for ti in t])
        I_cap_t = np.array([self.capacitive_current(V[i], t[i]) for i in range(len(t))])
        d_t = np.array([self.membrane_thickness(ti) for ti in t])
        
        results = {
            'time': t,
            'voltage': V,
            'gating_m': m,
            'gating_h': h, 
            'gating_n': n,
            'pressure': P_t,
            'capacitance': C_m_t,
            'capacitive_current': I_cap_t,
            'membrane_thickness': d_t,
            'success': True
        }
        
        logger.info("Simulation completed successfully")
        return results
    
    def plot_results(self, results, save_path=None):
        """Plot the key simulation results as specified in the PDF"""
        
        if not results['success']:
            logger.error("Cannot plot - simulation failed")
            return None
        
        fig, axes = plt.subplots(4, 1, figsize=(12, 14))
        
        t_ms = results['time'] * 1000  # Convert to ms
        
        # Plot 1: Pressure profile
        axes[0].plot(t_ms, results['pressure']/1000, 'b-', linewidth=2, label='P_US(t)')
        axes[0].set_ylabel('Pressure (kPa)')
        axes[0].set_title('Input Pressure Profile')
        axes[0].grid(True, alpha=0.3)
        axes[0].legend()
        
        # Plot 2: Capacitive current (the key mechanism)
        axes[1].plot(t_ms, results['capacitive_current']*1e6, 'r-', linewidth=2, label='I_cap(t)')
        axes[1].set_ylabel('Capacitive Current (μA/m²)')
        axes[1].set_title('Capacitive Current: V(t) × dC_m/dt')
        axes[1].grid(True, alpha=0.3)
        axes[1].legend()
        
        # Plot 3: Membrane potential response
        axes[2].plot(t_ms, results['voltage']*1000, 'g-', linewidth=2, label='V(t)')
        axes[2].set_ylabel('Membrane Potential (mV)')
        axes[2].set_title('Membrane Potential Response')
        axes[2].grid(True, alpha=0.3)
        axes[2].legend()
        
        # Plot 4: Membrane capacitance
        axes[3].plot(t_ms, results['capacitance']*1e6, 'purple', linewidth=2, label='C_m(t)')
        axes[3].set_ylabel('Capacitance (μF/cm²)')
        axes[3].set_xlabel('Time (ms)')
        axes[3].set_title('Time-Varying Membrane Capacitance')
        axes[3].grid(True, alpha=0.3)
        axes[3].legend()
        
        plt.tight_layout()
        
        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Results figure saved to {save_path}")
        
        return fig
    
    def pressure_response_curve(self, pressure_range=None, save_path=None):
        """Generate pressure-response curve by varying P_max"""
        if pressure_range is None:
            pressure_range = np.linspace(10e3, 200e3, 20)  # 10-200 kPa
        
        responses = []
        original_P_max = self.P_max
        
        logger.info("Generating pressure-response curve...")
        
        for P_max in pressure_range:
            self.P_max = P_max
            
            # Run simulation
            results = self.simulate(t_span=(0, 150e-3))
            
            if results['success']:
                # Check for action potential (threshold crossing)
                V_max = np.max(results['voltage'])
                spike_threshold = -55e-3  # -55 mV threshold
                fired = V_max > spike_threshold
                responses.append(1 if fired else 0)
            else:
                responses.append(0)
        
        # Restore original pressure
        self.P_max = original_P_max
        
        # Plot pressure-response curve
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(pressure_range/1000, responses, 'bo-', linewidth=2, markersize=8)
        ax.set_xlabel('Peak Pressure (kPa)')
        ax.set_ylabel('Action Potential Fired (0/1)')
        ax.set_title('Pressure-Response Curve')
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-0.1, 1.1)
        
        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Pressure-response curve saved to {save_path}")
        
        return pressure_range, responses, fig


def main():
    """Main function demonstrating the ultrasonic membrane simulation"""
    
    print("=" * 60)
    print("ULTRASONIC MEMBRANE CAPACITANCE SIMULATION")  
    print("Based on URENIMOD_simulation_RadiationForce_v1.pdf")
    print("=" * 60)
    
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
    
    print(f"Key Parameters:")
    print(f"  Peak pressure: {P_max/1000:.1f} kPa")
    print(f"  Pulse center: {t_centre*1000:.1f} ms")
    print(f"  Pulse width: {sigma_t*1000:.1f} ms")
    print(f"  Axon diameter: {diameter*1e6:.1f} μm")
    print(f"  Young's modulus: {E/1e6:.1f} MPa")
    print("")
    
    # Create simulator with defined parameters
    sim = UltrasonicMembraneSimulator(
        P_max=P_max, t_centre=t_centre, sigma_t=sigma_t,
        length=length, diameter=diameter, d0=d0, E=E,
        C_m0=C_m0, g_Na=g_Na, g_K=g_K, g_L=g_L,
        E_Na=E_Na, E_K=E_K, E_L=E_L, V_rest=V_rest
    )
    
    # Run simulation
    results = sim.simulate(t_span=(0, 100e-3))  # 100 ms simulation
    
    if results['success']:
        # Create figures directory
        os.makedirs('figures', exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Plot main results
        main_fig_path = os.path.join('figures', f'membrane_simulation_{timestamp}.png')
        fig1 = sim.plot_results(results, save_path=main_fig_path)
        
        # Generate pressure-response curve
        curve_fig_path = os.path.join('figures', f'pressure_response_{timestamp}.png')
        pressure_range, responses, fig2 = sim.pressure_response_curve(save_path=curve_fig_path)
        
        # Show summary
        V_max = np.max(results['voltage']) * 1000  # Convert to mV
        fired = V_max > -20  # -20 mV threshold
        
        print(f"\nSimulation Results:")
        print(f"  Peak pressure: {sim.P_max/1000:.1f} kPa")
        print(f"  Maximum voltage: {V_max:.1f} mV")
        print(f"  Action potential fired: {'Yes' if fired else 'No'}")
        print(f"  Max capacitive current: {np.max(np.abs(results['capacitive_current']))*1e6:.2f} μA/m²")
        
        plt.show()
        
        return sim, results
    
    else:
        print("Simulation failed!")
        return None, None


if __name__ == "__main__":
    main() 