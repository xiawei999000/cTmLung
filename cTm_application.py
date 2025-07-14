import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import numpy as np
# import os
import sys
# from scipy import ndimage
# import math

# Check and import required packages
try:
    import SimpleITK as sitk
except ImportError:
    print("Error: SimpleITK is not installed. Please install it using: pip install SimpleITK")
    sys.exit(1)

try:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib.figure import Figure
except ImportError:
    print("Error: matplotlib is not installed. Please install it using: pip install matplotlib")
    sys.exit(1)


class LungNoduleAnalyzer:
    def __init__(self, root):
        self.root = root
        self.root.title("Modified clinical T staging for early-stage LUAD")
        self.root.geometry("1600x900")

        # Data storage
        self.lung_image = None
        self.nodule_mask = None
        self.lung_data = None
        self.mask_data = None
        self.voxel_volume = None
        self.voxel_spacing = None

        # Staging options
        self.staging_options = ["cT1a", "cT1b", "cT1c", "cT2a"]

        self.create_widgets()

    def create_widgets(self):
        # Main frame
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Create left and right frames
        left_frame = ttk.Frame(main_frame)
        left_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), padx=(0, 10))

        right_frame = ttk.Frame(main_frame)
        right_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Create analysis section (left side)
        self.create_analysis_section(left_frame)

        # Create visualization section (right side)
        self.create_visualization_section(right_frame)

        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(0, weight=1)

    def create_analysis_section(self, parent):
        # File loading section
        file_frame = ttk.LabelFrame(parent, text="File Loading", padding="10")
        file_frame.grid(row=0, column=0, sticky=(tk.W, tk.E), pady=(0, 10))

        # Lung image loading
        ttk.Label(file_frame, text="Lung CT volume :").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.lung_file_var = tk.StringVar()
        ttk.Entry(file_frame, textvariable=self.lung_file_var, width=40).grid(row=0, column=1, padx=(5, 5))
        ttk.Button(file_frame, text="Browse", command=self.load_lung_image).grid(row=0, column=2)

        # Nodule mask loading
        ttk.Label(file_frame, text="Nodule Mask volume :").grid(row=1, column=0, sticky=tk.W, pady=2)
        self.mask_file_var = tk.StringVar()
        ttk.Entry(file_frame, textvariable=self.mask_file_var, width=40).grid(row=1, column=1, padx=(5, 5))
        ttk.Button(file_frame, text="Browse", command=self.load_mask_image).grid(row=1, column=2)

        # Help text
        help_text = "Instructions: Load CT image and mask, set parameters, then click 'Start Quantitative Analysis'\nNote: volume data with nii, nii.gz, or nrrd format recommended."
        ttk.Label(file_frame, text=help_text, font=("Arial", 8), foreground="gray").grid(
            row=2, column=0, columnspan=4, pady=(5, 0))

        # Parameters section
        param_frame = ttk.LabelFrame(parent, text="Parameters", padding="10")
        param_frame.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=(0, 10))

        # cT staging selection
        ttk.Label(param_frame, text="Current cT Staging (9th edition):").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.staging_var = tk.StringVar(value="") # 'cT1a'
        staging_combo = ttk.Combobox(param_frame, textvariable=self.staging_var, values=self.staging_options,
                                     state="readonly", width=10)
        staging_combo.grid(row=0, column=1, padx=(5, 5))

        # Threshold settings
        ttk.Label(param_frame, text="Ground Glass Threshold (TGGO):").grid(row=1, column=0, sticky=tk.W, pady=2)
        self.tggo_var = tk.StringVar(value="-800")
        ttk.Entry(param_frame, textvariable=self.tggo_var, width=10).grid(row=1, column=1, padx=(5, 20))

        ttk.Label(param_frame, text="Solid Threshold (TS):").grid(row=1, column=2, sticky=tk.W, pady=2)
        self.ts_var = tk.StringVar(value="-250")
        ttk.Entry(param_frame, textvariable=self.ts_var, width=10).grid(row=1, column=3, padx=(5, 5))

        # Help text
        help_text = "Note: TGGO < TS. Solid: [TS, +∞), Ground glass: [TGGO, TS)"
        ttk.Label(param_frame, text=help_text, font=("Arial", 8), foreground="gray").grid(
            row=2, column=0, columnspan=4, pady=(5, 0))

        # Calculate button
        calc_button = ttk.Button(parent, text="Start Quantitative Analysis", command=self.calculate_analysis)
        calc_button.grid(row=2, column=0, pady=10)

        # Results display section
        result_frame = ttk.LabelFrame(parent, text="Analysis Results", padding="10")
        result_frame.grid(row=3, column=0, sticky=(tk.W, tk.E), pady=(0, 10))

        # Results text box
        self.result_text = tk.Text(result_frame, height=25, width=70, font=("Courier", 8))
        scrollbar = ttk.Scrollbar(result_frame, orient="vertical", command=self.result_text.yview)
        self.result_text.configure(yscrollcommand=scrollbar.set)
        self.result_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))

        # Staging recommendation section
        staging_frame = ttk.LabelFrame(parent, text="Staging Recommendation", padding="10")
        staging_frame.grid(row=4, column=0, sticky=(tk.W, tk.E), pady=(0, 10))

        self.staging_label = ttk.Label(staging_frame, text="", font=("Arial", 11, "bold"))
        self.staging_label.grid(row=0, column=0, pady=5)

        # Help text
        help_text = "Note: a modified clinical T category (cTm) was obtained by downstaging the patients with SM%<45% and upstaging the patients with SM%>75%."
        ttk.Label(staging_frame, text=help_text, font=("Arial", 8), foreground="gray").grid(
            row=2, column=0, columnspan=4, pady=(5, 0))

        # Configure grid weights for analysis section
        parent.columnconfigure(0, weight=1)
        result_frame.columnconfigure(0, weight=1)
        result_frame.rowconfigure(0, weight=1)

    def create_visualization_section(self, parent):
        # Visualization frame
        viz_frame = ttk.LabelFrame(parent, text="Visualization", padding="10")
        viz_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Create matplotlib figure
        self.fig = Figure(figsize=(8, 8), dpi=100)
        self.ax1 = self.fig.add_subplot(211)  # Top subplot
        self.ax2 = self.fig.add_subplot(212)  # Bottom subplot

        # Create canvas
        self.canvas = FigureCanvasTkAgg(self.fig, viz_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Add toolbar for matplotlib
        toolbar_frame = ttk.Frame(viz_frame)
        toolbar_frame.grid(row=1, column=0, sticky=(tk.W, tk.E))

        # # Instructions
        # instructions = ttk.Label(viz_frame,
        #                          text="Instructions: Load CT image and mask, set parameters, then click 'Start Quantitative Analysis'",
        #                          font=("Arial", 9), foreground="gray")
        # instructions.grid(row=2, column=0, pady=(5, 0))

        # Configure grid weights for visualization section
        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(0, weight=1)
        viz_frame.columnconfigure(0, weight=1)
        viz_frame.rowconfigure(0, weight=1)

        # Initialize with empty plots
        self.ax1.text(0.5, 0.5, 'Original CT Image\n(Load data)',
                      ha='center', va='center', transform=self.ax1.transAxes, fontsize=12)
        self.ax1.set_xticks([])
        self.ax1.set_yticks([])

        self.ax2.text(0.5, 0.5, 'Component Overlay\n(Red: Solid, Green: GGO)',
                      ha='center', va='center', transform=self.ax2.transAxes, fontsize=12)
        self.ax2.set_xticks([])
        self.ax2.set_yticks([])

        self.canvas.draw()

    def load_lung_image(self):
        file_path = filedialog.askopenfilename(
            title="Select Lung Image File",
            filetypes=[("All files", "*.*"), ("NIfTI files", "*.nii.gz"), ("NIfTI files", "*.nii")]
        )
        if file_path:
            self.lung_file_var.set(file_path)
            try:
                self.lung_image = sitk.ReadImage(file_path)
                self.lung_data = sitk.GetArrayFromImage(self.lung_image)
                self.voxel_spacing = self.lung_image.GetSpacing()
                self.voxel_volume = np.prod(self.voxel_spacing)
                messagebox.showinfo("Success", "Lung image loaded successfully!")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load lung image: {str(e)}")
                self.lung_image = None
                self.lung_data = None

    def load_mask_image(self):
        file_path = filedialog.askopenfilename(
            title="Select Nodule Mask File",
            filetypes=[("All files", "*.*"), ("NIfTI files", "*.nii.gz"), ("NIfTI files", "*.nii")]
        )
        if file_path:
            self.mask_file_var.set(file_path)
            try:
                self.nodule_mask = sitk.ReadImage(file_path)
                self.mask_data = sitk.GetArrayFromImage(self.nodule_mask)
                messagebox.showinfo("Success", "Nodule mask loaded successfully!")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load nodule mask: {str(e)}")
                self.nodule_mask = None
                self.mask_data = None

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

    def adjust_staging(self, current_staging, sm_percent):
        """Adjust staging based on SM%"""
        current_index = self.staging_options.index(current_staging)

        if sm_percent < 45:  # Down-staging
            new_index = max(0, current_index - 1)
            adjustment = "DOWN-STAGING"
        elif sm_percent > 75:  # Up-staging
            new_index = min(len(self.staging_options) - 1, current_index + 1)
            adjustment = "UP-STAGING"
        else:  # No change
            new_index = current_index
            adjustment = "NO CHANGE"

        return self.staging_options[new_index], adjustment

    def validate_inputs(self):
        """Validate inputs"""
        if self.lung_data is None or self.mask_data is None:
            messagebox.showerror("Error", "Please load both lung image and nodule mask files!")
            return False

        if self.lung_data.shape != self.mask_data.shape:
            messagebox.showerror("Error", "Lung image and mask must have the same dimensions!")
            return False

        # 添加这个检查
        if not self.staging_var.get():
            messagebox.showerror("Error", "Please select a cT staging!")
            return False

        try:
            tggo = float(self.tggo_var.get())
            ts = float(self.ts_var.get())

            if tggo >= ts:
                messagebox.showerror("Error", "Ground glass threshold must be less than solid threshold!")
                return False

            return True
        except ValueError:
            messagebox.showerror("Error", "Please enter valid numerical values for thresholds!")
            return False

    def update_visualization(self, max_slice, nodule_mask_bool, solid_mask, ggo_mask):
        """Update visualization with original and overlay images"""
        # Clear previous plots
        self.ax1.clear()
        self.ax2.clear()

        # Get the slice data
        lung_slice = self.lung_data[max_slice]
        nodule_slice = nodule_mask_bool[max_slice]
        solid_slice = solid_mask[max_slice]
        ggo_slice = ggo_mask[max_slice]

        # Display original image
        self.ax1.imshow(lung_slice, cmap='gray', vmin=-1000, vmax=400)
        self.ax1.contour(nodule_slice, colors='yellow', linewidths=2)
        self.ax1.set_title(f'Original CT Image with Nodule Contour (Slice {max_slice})', fontsize=12)
        self.ax1.axis('off')

        # Display overlay image
        self.ax2.imshow(lung_slice, cmap='gray', vmin=-1000, vmax=400)

        # Create colored overlays
        solid_overlay = np.zeros((*solid_slice.shape, 4))
        ggo_overlay = np.zeros((*ggo_slice.shape, 4))

        # Red for solid (semi-transparent)
        solid_overlay[solid_slice > 0] = [1, 0, 0, 0.6]

        # Green for GGO (semi-transparent)
        ggo_overlay[ggo_slice > 0] = [0, 1, 0, 0.6]

        self.ax2.imshow(solid_overlay)
        self.ax2.imshow(ggo_overlay)
        self.ax2.set_title('CT Image with Component Overlay (Red: Solid, Green: GGO)', fontsize=12)
        self.ax2.axis('off')

        # Add statistics text
        solid_pixels = np.sum(solid_slice > 0)
        ggo_pixels = np.sum(ggo_slice > 0)
        total_pixels = np.sum(nodule_slice > 0)

        stats_text = f'Solid pixels: {solid_pixels}\nGGO pixels: {ggo_pixels}\nTotal pixels: {total_pixels}'
        self.ax2.text(0.02, 0.98, stats_text, transform=self.ax2.transAxes,
                      verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        # Adjust layout
        self.fig.tight_layout()

        # Refresh canvas
        self.canvas.draw()

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

    def display_comprehensive_results(self, tggo, ts, current_staging,
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

        # Clear and display results
        self.result_text.delete(1.0, tk.END)
        self.result_text.insert(1.0, results)

    def show_staging_recommendation(self, current_staging, new_staging, adjustment, sm_percent):
        if adjustment == "DOWN-STAGING":
            recommendation = f"SM% = {sm_percent:.0f}% < 45%\n{adjustment} Recommended: {current_staging} → {new_staging}"
            color = "blue"
        elif adjustment == "UP-STAGING":
            recommendation = f"SM% = {sm_percent:.0f}% > 75%\n{adjustment} Recommended: {current_staging} → {new_staging}"
            color = "red"
        else:
            recommendation = f"SM% = {sm_percent:.0f}% (45% ≤ SM% ≤ 75%)\n{adjustment}: Maintain {current_staging}"
            color = "green"

        self.staging_label.config(text=recommendation, foreground=color)


def check_dependencies():
    """Check if all required dependencies are available"""
    required_packages = ['SimpleITK', 'numpy', 'matplotlib', 'scipy', 'tkinter']
    missing_packages = []

    try:
        import SimpleITK
    except ImportError:
        missing_packages.append('SimpleITK')

    try:
        import numpy
    except ImportError:
        missing_packages.append('numpy')

    try:
        import matplotlib
    except ImportError:
        missing_packages.append('matplotlib')

    try:
        import scipy
    except ImportError:
        missing_packages.append('scipy')

    try:
        import tkinter
    except ImportError:
        missing_packages.append('tkinter')

    if missing_packages:
        print("Missing required packages:")
        for package in missing_packages:
            print(f"  - {package}")
        print("\nPlease install missing packages using:")
        print("pip install SimpleITK numpy matplotlib scipy")
        if 'tkinter' in missing_packages:
            print("Note: tkinter usually comes with Python.")
        return False

    return True


def main():
    if not check_dependencies():
        input("Press Enter to exit...")
        return

    try:
        root = tk.Tk()
        app = LungNoduleAnalyzer(root)
        root.mainloop()
    except Exception as e:
        print(f"Error starting application: {str(e)}")
        input("Press Enter to exit...")


if __name__ == "__main__":
    main()
