import numpy as np

def calculate_entropy(self, hu_values):
    """Calculate entropy of HU values"""
    if len(hu_values) == 0:
        return 0

    # Create histogram
    hist, bins = np.histogram(hu_values, bins=50)
    hist = hist[hist > 0]  # Remove zero counts

    # Calculate probabilities
    probabilities = hist / np.sum(hist)

    # Calculate entropy
    entropy = -np.sum(probabilities * np.log2(probabilities))
    return entropy

def find_max_cross_sectional_area(self, mask_3d):
    """Find the slice with maximum cross-sectional area"""
    max_area = 0
    max_slice = 0

    for i in range(mask_3d.shape[0]):
        area = np.sum(mask_3d[i] > 0)
        if area > max_area:
            max_area = area
            max_slice = i

    return max_slice, max_area

def calculate_max_diameter_2d(self, mask_2d, spacing):
    """Calculate maximum diameter in 2D slice"""
    if np.sum(mask_2d) == 0:
        return 0

    # Find coordinates of mask pixels
    coords = np.where(mask_2d > 0)
    if len(coords[0]) == 0:
        return 0

    # Convert to physical coordinates
    y_coords = coords[0] * spacing[1]  # spacing[1] is y-spacing
    x_coords = coords[1] * spacing[0]  # spacing[0] is x-spacing

    # Calculate all pairwise distances
    max_distance = 0
    for i in range(len(y_coords)):
        for j in range(i + 1, len(y_coords)):
            distance = np.sqrt((y_coords[i] - y_coords[j]) ** 2 + (x_coords[i] - x_coords[j]) ** 2)
            max_distance = max(max_distance, distance)

    return max_distance

def calculate_max_diameter_3d(self, mask_3d, spacing):
    """Calculate maximum diameter in 3D"""
    if np.sum(mask_3d) == 0:
        return 0

    # Find coordinates of mask voxels
    coords = np.where(mask_3d > 0)
    if len(coords[0]) == 0:
        return 0

    # Convert to physical coordinates
    z_coords = coords[0] * spacing[2]  # spacing[2] is z-spacing
    y_coords = coords[1] * spacing[1]  # spacing[1] is y-spacing
    x_coords = coords[2] * spacing[0]  # spacing[0] is x-spacing

    # Sample points for efficiency (if too many points)
    n_points = len(z_coords)
    if n_points > 1000:  # Sample for efficiency
        indices = np.random.choice(n_points, 1000, replace=False)
        z_coords = z_coords[indices]
        y_coords = y_coords[indices]
        x_coords = x_coords[indices]

    # Calculate all pairwise distances
    max_distance = 0
    for i in range(len(z_coords)):
        for j in range(i + 1, len(z_coords)):
            distance = np.sqrt((z_coords[i] - z_coords[j]) ** 2 +
                               (y_coords[i] - y_coords[j]) ** 2 +
                               (x_coords[i] - x_coords[j]) ** 2)
            max_distance = max(max_distance, distance)

    return max_distance

def calculate_mass(self, hu_values, voxel_count):
    """Calculate mass: Mass = Volume × [(Attenuation + 1,000) × 0.001]"""
    if len(hu_values) == 0:
        return 0, 0, 0

    # Calculate mean attenuation
    mean_attenuation = np.mean(hu_values)

    # Calculate volume
    volume = voxel_count * self.voxel_volume

    # Calculate mass
    mass = volume * ((mean_attenuation + 1000) * 0.001)

    return mass, mean_attenuation, volume


def calculate_analysis(self):
    if not self.validate_inputs():
        return

    try:
        # Get parameters
        tggo = float(self.tggo_var.get())
        ts = float(self.ts_var.get())
        current_staging = self.staging_var.get()

        # Get masks
        nodule_mask_bool = self.mask_data > 0
        solid_mask = (self.lung_data >= ts) & nodule_mask_bool
        ggo_mask = (self.lung_data >= tggo) & (self.lung_data < ts) & nodule_mask_bool

        # Find maximum cross-sectional area slice
        max_slice, max_area = self.find_max_cross_sectional_area(nodule_mask_bool)

        # Get HU values
        nodule_hu_values = self.lung_data[nodule_mask_bool]
        solid_hu_values = self.lung_data[solid_mask]
        ggo_hu_values = self.lung_data[ggo_mask]

        # Calculate basic metrics
        whole_voxel_count = np.sum(nodule_mask_bool)
        solid_voxel_count = np.sum(solid_mask)
        ggo_voxel_count = np.sum(ggo_mask)

        # Calculate masses and volumes
        mass_whole, mean_hu_whole, volume_whole = self.calculate_mass(nodule_hu_values, whole_voxel_count)
        mass_solid, mean_hu_solid, volume_solid = self.calculate_mass(solid_hu_values, solid_voxel_count)
        mass_ggo, mean_hu_ggo, volume_ggo = self.calculate_mass(ggo_hu_values, ggo_voxel_count)

        # Calculate percentages
        sm_percent = (mass_solid / mass_whole * 100) if mass_whole > 0 else 0
        gm_percent = (mass_ggo / mass_whole * 100) if mass_whole > 0 else 0
        sv_percent = (volume_solid / volume_whole * 100) if volume_whole > 0 else 0
        gv_percent = (volume_ggo / volume_whole * 100) if volume_whole > 0 else 0

        # Calculate size metrics on axial image (max cross-sectional slice)
        nodule_2d = nodule_mask_bool[max_slice]
        solid_2d = solid_mask[max_slice]
        ggo_2d = ggo_mask[max_slice]

        # Calculate 2D diameters
        mtsa = self.calculate_max_diameter_2d(nodule_2d, self.voxel_spacing)
        mssa = self.calculate_max_diameter_2d(solid_2d, self.voxel_spacing)
        mgsa = self.calculate_max_diameter_2d(ggo_2d, self.voxel_spacing)

        # Calculate 3D diameters
        mtsv = self.calculate_max_diameter_3d(nodule_mask_bool, self.voxel_spacing)
        mssv = self.calculate_max_diameter_3d(solid_mask, self.voxel_spacing)
        mgsv = self.calculate_max_diameter_3d(ggo_mask, self.voxel_spacing)

        # Calculate ratios
        mssa_percent = (mssa / mtsa * 100) if mtsa > 0 else 0

        # Calculate entropy
        entropy_hu = self.calculate_entropy(nodule_hu_values)

        # Display results
        self.display_comprehensive_results(
            tggo, ts, current_staging,
            mtsa, mssa, mgsa, mtsv, mssv, mgsv,
            volume_whole, volume_solid, volume_ggo,
            mass_whole, mass_solid, mass_ggo,
            mean_hu_whole, entropy_hu,
            sm_percent, gm_percent, sv_percent, gv_percent, mssa_percent
        )

        # Update visualization
        self.update_visualization(max_slice, nodule_mask_bool, solid_mask, ggo_mask)

        # Staging recommendation
        new_staging, adjustment = self.adjust_staging(current_staging, sm_percent)
        self.show_staging_recommendation(current_staging, new_staging, adjustment, sm_percent)

    except Exception as e:
        messagebox.showerror("Error", f"Error during calculation: {str(e)}")

def display_comprehensive_results(tggo, ts, current_staging,
                                  mtsa, mssa, mgsa, mtsv, mssv, mgsv,
                                  tv, sv, gv, tm, sm, gm, mhu, ehu,
                                  sm_percent, gm_percent, sv_percent, gv_percent, mssa_percent):

            results = f"""
    === COMPREHENSIVE LUNG NODULE ANALYSIS ===
    
    Parameters:
      Ground Glass Threshold (TGGO): {tggo} HU
      Solid Threshold (TS): {ts} HU
      Current cT Staging: {current_staging}
    
    Size Measurements on Axial Image:
      Maximal Total Size on Axial (MTSA): {mtsa:.2f} mm
      Maximal Solid Size on Axial (MSSA): {mssa:.2f} mm
      MSSA-to-MTSA Ratio (MSSA%): {mssa_percent:.2f}%
      Maximal GGO Size on Axial (MGSA): {mgsa:.2f} mm
    
    Size Measurements by Volume:
      Maximal Total Size by Volume (MTSV): {mtsv:.2f} mm
      Maximal Solid Size by Volume (MSSV): {mssv:.2f} mm
      Maximal GGO Size by Volume (MGSV): {mgsv:.2f} mm
    
    Volume Analysis:
      Total Volume (TV): {tv:.2f} mm³
      Volume of Solid Component (SV): {sv:.2f} mm³
      Volume of GGO Component (GV): {gv:.2f} mm³
      SV-to-TV Ratio (SV%): {sv_percent:.2f}%
      GV-to-TV Ratio (GV%): {gv_percent:.2f}%
    
    Mass Analysis:
      Total Mass (TM): {tm:.2f} mg
      Solid Mass (SM): {sm:.2f} mg
      GGO Mass (GM): {gm:.2f} mg
      SM-to-TM Ratio (SM%): {sm_percent:.0f}%
      GM-to-TM Ratio (GM%): {gm_percent:.0f}%
    
    Intensity Analysis:
      Mean HU (MHU): {mhu:.2f} HU
      Entropy of HU (EHU): {ehu:.4f}
    
    Key Clinical Metrics:
      Solid Mass Percentage (SM%): {sm_percent:.0f}%
    """

