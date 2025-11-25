from scipy.optimize import minimize
import pandas as pd
import numpy as np
import math


#=============================== LAMMPS cossine-periodic to gromacs harmonic ===============================
def get_lammps_cos_per_params(C, B, n):                    
    theta = np.linspace(0, 2 * np.pi, 1000)
    theta_angle = theta * 180 / np.pi

    if B == 1:
        periodic = 3   
    else:
        periodic = 2

    theta0 = 180 / n * (periodic)

    E_list = []
    for i in range(len(theta)):
        E = 2 / (n**2) * C * (1 - (B * (-1**n) * (np.cos(n * theta[i]))))
        E_list.append(E)

    # Define the range for fitting
    fit_range_start = theta0 - 100/n
    fit_range_end = theta0 + 100/n

    # Filter data within the fitting range
    theta_fit = theta_angle[(theta_angle >= fit_range_start) & (theta_angle <= fit_range_end)]
    E_fit = np.array(E_list)[(theta_angle >= fit_range_start) & (theta_angle <= fit_range_end)]

    # Define the objective function to minimize (sum of squared differences)
    def objective_function(keff):
        E_harm_fit = keff * (theta_fit - theta0)**2
        return np.sum((E_fit - E_harm_fit)**2)

    # Find the optimal keff using optimization
    initial_keff = C / 4000  # Starting point for optimization
    result = minimize(objective_function, initial_keff)
    optimal_keff = result.x[0]

    return optimal_keff, theta0


#=============================== LAMMPS fourier to gromacs gromos style ===============================
def get_fourier_gromos_params(K, C0, C1, C2):
    theta = np.linspace(0, 2 * np.pi, 1000)
    theta_angle = theta * 180 / np.pi

    E_list = []
    for i in range(len(theta)):
        E = K * (C0 + C1 * (np.cos(theta[i])) + C2 * (np.cos(2 * theta[i])))
        E_list.append(E)

    # For the Fourier potential, the equilibrium angle is where the derivative is zero.
    # Finding the exact equilibrium angle can be complex, so we'll use the minimum of the calculated energies as a proxy for theta0.
    min_energy_index = np.argmin(E_list)
    theta0_approx = theta_angle[min_energy_index]

    # Define a fitting range around the approximate equilibrium angle
    # Adjust the range based on the expected width of the potential well
    # For simplicity here, we'll assume fitting over the whole range for now, similar to the previous function.
    theta_fit = theta_angle
    E_fit = np.array(E_list)

    # Convert theta_fit back to radians for the new equation
    theta_fit_rad = theta_fit * np.pi / 180

    # Define the objective function to minimize (sum of squared differences)
    def objective_function(params):
        k, theta_eq_rad = params
        # Ensure theta_eq_rad is within a valid range for np.cos
        theta_eq_rad = np.arctan2(np.sin(theta_eq_rad), np.cos(theta_eq_rad))

        E_new_eq_fit = k/2 * (np.cos(theta_fit_rad) - np.cos(theta_eq_rad))**2
        return np.sum((E_fit - E_new_eq_fit)**2)

    # Find the optimal k and theta_eq using optimization
    # Initial guess for k and theta_eq (in radians)
    initial_k = K / 1000  # Starting point for optimization, adjust as needed
    initial_theta_eq_rad = theta0_approx * np.pi / 180 # Starting point for optimization

    initial_params = [initial_k, initial_theta_eq_rad]

    # Use bounds to keep theta_eq_rad within a reasonable range (e.g., -2*pi to 2*pi) to help the optimizer
    bounds = [(None, None), (-2 * np.pi, 2 * np.pi)]


    result = minimize(objective_function, initial_params, bounds=bounds)
    optimal_k, optimal_theta_eq_rad = result.x
    optimal_theta_eq_angle = (optimal_theta_eq_rad * 180 / np.pi) % 360 # Ensure angle is within 0-360

    return optimal_k, optimal_theta_eq_angle

#=============================== LAMMPS fourier to gromacs harmonic ===============================
def get_fourier_harmonic_params(K_fourier, C0, C1, C2):
    """
    Converts LAMMPS Fourier (kcal/mol) → GROMACS harmonic (kJ/mol/rad^2)
    by fitting the harmonic potential near theta = 0 (radians).
    """

    # --- 1. Angles for fitting (in RADIANS) ---
    theta = np.linspace(-np.pi, np.pi, 1000)
    theta_deg = theta * 180 / np.pi  # only for defining the fit interval

    # --- 2. Compute LAMMPS Fourier energy in kcal/mol ---
    E_fourier = K_fourier * (C0 + C1*np.cos(theta) + C2*np.cos(2*theta))

    # --- 3. Convert energies to kJ/mol BEFORE fitting ---
    E_fourier = E_fourier * 4.184  # kcal → kJ

    # --- 4. Select fitting range: -60° to +60° ---
    mask = (theta_deg >= -60) & (theta_deg <= 60)
    theta_fit = theta[mask]      # radians
    E_fit = E_fourier[mask]      # kJ/mol

    # --- 5. Harmonic reference angle ---
    theta0 = 0.0  # radians

    # --- 6. Objective function: least squares ---
    def objective(k):
        E_harm = 0.5 * k * (theta_fit - theta0)**2
        return np.sum((E_fit - E_harm)**2)

    # --- 7. Initial guess: correct scale (kJ/mol/rad²) ---
    initial_k = K_fourier * 4.184

    # --- 8. Optimize k ---
    result = minimize(objective, x0=[initial_k], bounds=[(0, None)])
    k_opt = float(result.x[0])  # already in kJ/mol/rad²

    return k_opt, theta0
