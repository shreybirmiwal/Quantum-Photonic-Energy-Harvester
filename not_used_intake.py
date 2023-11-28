##########################INTAKE#######################################

# Planck distribution function
def planck_distribution(lam):
    return (2 * h * c**2 / lam**5) / (np.exp(h * c / (lam * k * T)) - 1)

# Integrate Planck's distribution over the wavelength range
total_radiance, _ = integrate.quad(planck_distribution, lam_min, lam_max)


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










##########################QPC#######################################

# Calculation of the number of photons
def calculate_num_photons(total_energy, average_photon_energy, area):
    energy_per_second = total_energy * area
    return energy_per_second / average_photon_energy

panel_azimuth = 180  # Example value, facing south
panel_tilt = 30  # Example value, 30 degrees tilt from the horizontal\

# Average photon energy for a range of wavelengths
average_wavelength = (lam_min + lam_max) / 2
average_photon_energy = h * c / average_wavelength