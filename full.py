import math
import pandas as pd
import scipy.constants as const
import scipy.integrate as integrate
from scipy.interpolate import interp1d
import numpy as np
from scipy.integrate import solve_ivp

# Constants
h = const.h  # Planck constant
c = const.c  # Speed of light
k = const.k  # Boltzmann constant
T = 5778  # Temperature of the Sun's surface in Kelvin
lam_min = 300e-9  # Minimum wavelength in meters
lam_max = 2500e-9  # Maximum wavelength in meters

# Planck distribution function
def planck_distribution(lam):
    return (2 * h * c**2 / lam**5) / (np.exp(h * c / (lam * k * T)) - 1)

# Integrate Planck's distribution over the wavelength range
total_radiance, _ = integrate.quad(planck_distribution, lam_min, lam_max)

# Average photon energy for a range of wavelengths
average_wavelength = (lam_min + lam_max) / 2
average_photon_energy = h * c / average_wavelength

# Solar angle calculation functions
def solar_declination(d, h):
    fracyear = 2 * math.pi / 365 * (d - 1 + (h - 12) / 24)
    delta = math.degrees(0.006918 - 0.399912 * math.cos(fracyear) + 0.070257 * math.sin(fracyear) - 0.006758 * math.cos(2 * fracyear) + 0.000907 * math.sin(2 * fracyear) - 0.002697 * math.cos(3 * fracyear) + 0.00148 * math.sin(3 * fracyear))
    return delta

def calc_elevation_angle(d, ho, longitude=97.74, latitude=30.266666):
    h = ho + 6
    fracyear = 2 * math.pi / 365 * (d - 1 + (h - 12) / 24)
    eqtime = 229.17 * (0.000075 + 0.001868 * math.cos(fracyear) - 0.032077 * math.sin(fracyear) - 0.014615 * math.cos(2 * fracyear) - 0.040849 * math.sin(2 * fracyear))
    time_offset = eqtime - 4 * longitude
    tst = (h * 60 + time_offset) % 1440
    ha = math.radians(tst / 4 - 180) if tst < 720 else math.radians(tst / 4 - 180 - 360)
    delta = math.radians(solar_declination(d, h))
    zenith = math.degrees(math.acos(math.sin(math.radians(latitude)) * math.sin(delta) + math.cos(math.radians(latitude)) * math.cos(delta) * math.cos(ha)))
    psi = 90 - zenith
    return psi

def angle_of_incidence(solar_azimuth, solar_elevation, panel_azimuth, panel_tilt):
    """
    Calculate the angle of incidence of sunlight on the solar panel.
    solar_azimuth: Azimuth angle of the sun (degrees)
    solar_elevation: Elevation angle of the sun (degrees)
    panel_azimuth: Azimuth angle of the solar panel (degrees)
    panel_tilt: Tilt angle of the solar panel from the horizontal (degrees)
    """
    # Convert angles to radians
    solar_azimuth_rad = math.radians(solar_azimuth)
    solar_elevation_rad = math.radians(solar_elevation)
    panel_azimuth_rad = math.radians(panel_azimuth)
    panel_tilt_rad = math.radians(panel_tilt)

    # Calculate the angle of incidence
    cos_incidence = (math.sin(solar_elevation_rad) * math.cos(panel_tilt_rad) +
                     math.cos(solar_elevation_rad) * math.sin(panel_tilt_rad) *
                     math.cos(solar_azimuth_rad - panel_azimuth_rad))
    incidence_angle = math.degrees(math.acos(min(max(cos_incidence, -1), 1)))
    return incidence_angle

def air_mass_coefficient(elevation_angle):
    """
    Calculate the air mass coefficient using a simplified version of the Kasten-Young formula.
    This coefficient represents the path length that sunlight travels through the atmosphere.
    """
    if elevation_angle <= 0:
        return 0  # No sunlight during night or when the sun is below the horizon
    zenith_angle = 90 - elevation_angle
    zenith_angle_rad = math.radians(zenith_angle)
    return 1 / (math.cos(zenith_angle_rad) + 0.50572 * (96.07995 - zenith_angle)**(-1.6364))

def adjusted_sunlight_intensity(solar_elevation_angle):
    """
    Adjust the sunlight intensity based on the solar elevation angle and atmospheric effects.
    """
    am_coefficient = air_mass_coefficient(solar_elevation_angle)
    intensity_at_top_of_atmosphere = 1361  # Solar constant in W/m²
    concentration_factor = 2500
    if am_coefficient == 0:
        return 0
    return intensity_at_top_of_atmosphere * concentration_factor * math.cos(math.radians(90 - solar_elevation_angle)) / am_coefficient

def corrected_sunlight_intensity(incidence_angle, base_intensity):
    """
    Correct the base sunlight intensity based on the angle of incidence.
    incidence_angle: Angle of incidence (degrees)
    base_intensity: Intensity of sunlight before correction (W/m²)
    """
    # Calculate the intensity correction factor
    if incidence_angle > 90:
        return 0  # No sunlight is received if the incidence angle is more than 90 degrees
    correction_factor = math.cos(math.radians(incidence_angle))
    return base_intensity * correction_factor


# Calculation of the number of photons
def calculate_num_photons(total_energy, average_photon_energy, area):
    energy_per_second = total_energy * area
    return energy_per_second / average_photon_energy

# Energy-related calculations
file_path = '/Users/rohanshankar/Downloads/OpticalPropertiesOfSilicon.xlsx'  # Adjust the file path as needed
df = pd.read_excel(file_path)

# Extract wavelengths and absorption coefficients into arrays
# Assuming the Excel file has columns 'Wavelength' (in nanometers) and 'AbsorptionCoefficient' (in cm^-1)
wavelengths_nm = df['wavelength(nm)'].to_numpy()  # Wavelengths in nanometers
absorption_coefficients_cm = df['a(/cm)'].to_numpy()  # Absorption coefficients in cm^-1

# Convert wavelengths to meters and absorption coefficients to m^-1
wavelengths_m = wavelengths_nm * 1e-9  # Convert nanometers to meters
absorption_coefficients_m = absorption_coefficients_cm * 100  # Convert cm^-1 to m^-1

# Create an interpolation function for absorption coefficients
absorption_interpolator = interp1d(wavelengths_m, absorption_coefficients_m, kind='cubic', fill_value="extrapolate")

def silicon_absorption_coefficient(wavelength):
    """
    Calculate the absorption coefficient of silicon for a given wavelength using interpolation.
    wavelength: Wavelength of the incoming light (in meters)
    """
    return absorption_interpolator(wavelength)

# Update the photon calculation to include wavelength-dependent absorption
def calculate_photon_absorption(intensity, lam_min, lam_max):
    """
    Calculate the number of photons absorbed by silicon over the solar spectrum.
    intensity: Intensity of sunlight (W/m²)
    lam_min, lam_max: Wavelength range (in meters)
    """
    def integrand(lam):
        energy_per_photon = h * c / lam
        num_photons = intensity / energy_per_photon
        absorption_efficiency = silicon_absorption_coefficient(lam)
        return num_photons * absorption_efficiency

    absorbed_photons, _ = integrate.quad(integrand, lam_min, lam_max)
    return absorbed_photons

panel_azimuth = 180  # Example value, facing south
panel_tilt = 30  # Example value, 30 degrees tilt from the horizontal

# Solar simulation
results = []
for day in range(1, 8):  # Example: one week
    for hour in range(24):
        elevation_angle = calc_elevation_angle(day, hour)
        intensity = adjusted_sunlight_intensity(elevation_angle)  # Base sunlight intensity
        absorbed_photons = calculate_photon_absorption(intensity, lam_min, lam_max)
        exciton_generation_rate = absorbed_photons  # Assuming one exciton per absorbed photon

        # Store the results
        results.append([day, hour, elevation_angle, intensity, absorbed_photons, exciton_generation_rate])

# Convert results to a DataFrame for analysis and output
df_results = pd.DataFrame(results, columns=["Day", "Hour", "Elevation Angle", "Intensity", "Absorbed Photons", "Exciton Generation Rate"])

# Constants for CMC simulation
consts = {
    'N': 7,
    'ω': [1.13, 0.96, 0.89, 1.09, 0.98, 1.04, 0.86],
    'γ': [0.09780352, 0.2154919, 0.2159708, 0.17640289, 0.15794933, 0.09051265, 0.20707075],
    'γφ': [0.333, 0.292, 0.364, 0.228, 0.317, 0.251, 0.375],
    'Γ': [0.064, 0.045, 0.058, 0.077, 0.051, 0.068, 0.044],
    'J': [[0.45, 0.03, 0.38, 0.47, 0.29, 0.08, 0.41],
          [0.10, 0.40, 0.48, 0.11, 0.30, 0.17, 0.15],
          [0.01, 0.43, 0.19, 0.36, 0.05, 0.14, 0.22],
          [0.35, 0.33, 0.27, 0.25, 0.49, 0.44, 0.42],
          [0.39, 0.26, 0.18, 0.50, 0.07, 0.20, 0.16],
          [0.48, 0.37, 0.45, 0.25, 0.31, 0.17, 0.35],
          [0.02, 0.19, 0.34, 0.28, 0.46, 0.40, 0.21]]
}

# Pauli matrices for CMC simulation
pauli = {
    '+': np.array([[0, 1], [0, 0]]),
    '-': np.array([[0, 0], [1, 0]]),
    'z': np.array([[1, 0], [0, -1]])
}

# CMC simulation functions
def sigma(i, op):
    return np.kron(np.eye(2**i), np.kron(pauli[op], np.eye(2**(consts['N']-i-1))))

def H():
    N = consts['N']
    ω = consts['ω']
    J = consts['J']
    H_val = np.zeros((2**N, 2**N))
    for i in range(N):
        H_val += ω[i] * sigma(i, 'z')
        for j in range(i+1, N):
            H_val += J[i][j] * (sigma(i, '+') @ sigma(j, '-') + sigma(i, '-') @ sigma(j, '+'))
    return H_val

def anticomm(A, B, C):
    return np.dot(A, np.dot(B, C)) + np.dot(B, np.dot(A, C))

def decay(rho, i):
    return consts['γ'][i] * (sigma(i, '-') @ rho @ sigma(i, '+') - 0.5 * anticomm(sigma(i, '+'), sigma(i, '-'), rho))

def dephase(rho, i):
    return consts['γφ'][i] * (sigma(i, 'z') @ rho @ sigma(i, 'z') - rho)

def emit(rho, i):
    return consts['Γ'][i] * (sigma(i, '-') @ rho @ sigma(i, '+'))

def rhs(t, rho_vec):
    rho = rho_vec.reshape((2**consts['N'], 2**consts['N']))
    dρdt = -1j * np.dot(H(), rho) + 1j * np.dot(rho, H())
    for i in range(consts['N']):
        dρdt += decay(rho, i) + dephase(rho, i) + emit(rho, i)
    return dρdt.flatten()

def simulate_CMC(G_initial):
    rho0 = np.eye(2**consts['N']) / 2**consts['N']
    times = np.linspace(0, 10, 100)
    result = solve_ivp(rhs, (0, 10), rho0.flatten(), t_eval=times)
    P_excited = np.kron([[0, 0], [0, 1]], np.eye(2**(consts['N']-1)))
    G_final = G_initial * np.trace(np.dot(P_excited, result.y[:, -1].reshape((2**consts['N'], 2**consts['N']))))
    return G_final

# Apply CMC simulation to solar simulation results
df_results['Final_Exciton_Rate_After_CMC'] = df_results['Exciton Generation Rate'].apply(simulate_CMC)

conversion_factor = .689 #describes how efficiency of energy conversion from excetonic to electric (0 - 1)
energy_used = 0 #describes what % of the energy is used for bodily functions like respiration. 1 means that all absorbed energy is used for functioning. 0 means 0%
# convert the excitonic rate to an actual energy value

def exciton_rate_to_energy(excitonic_rate):
    power = excitonic_rate * 15.01e-3 * 1.60218e-19       #power = excitonic rate * coefficient * conversion to watts
    return power
    
# Converts excitonic energy to electric energy
# Excitonic energy is valuable but not directly usable by most appliances. We need to convert it into electrical
# energy. The efficiency of this conversion process is represented by the factor η, with values between 0 (0%
# efficient) and 1 (100% efficient)
# Storage Mechanism:
# Once converted, the energy is stored in specialized storage cells, similar to batteries. The amount of energy Es
# stored at any given time is:

def energyConversion(excitonic_energy):
    convertedEnergy = excitonic_energy*conversion_factor
    usableEnergy = convertedEnergy - (convertedEnergy * energy_used)
    return usableEnergy

df_results['Excitonic Energy'] = df_results['Final_Exciton_Rate_After_CMC'].apply(exciton_rate_to_energy)
df_results['Usable Electric Energy'] = df_results['Excitonic Energy'].apply(energyConversion)

# Save the results
df_results.to_csv('solar_energy_CMC_results.csv', index=False)