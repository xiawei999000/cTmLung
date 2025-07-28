import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import numpy as np
import sys

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
        self.root.geometry("1200x900")

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

        # Create analysis section
        self.create_analysis_section(main_frame)

        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(0, weight=1)
        main_frame.rowconfigure(0, weight=0)  # File loading
        main_frame.rowconfigure(1, weight=0)  # Parameters
        main_frame.rowconfigure(2, weight=0)  # Button
        main_frame.rowconfigure(3, weight=1)  # Visualization gets more space
        main_frame.rowconfigure(4, weight=0)  # Staging recommendation

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
        self.staging_var = tk.StringVar(value="")  # 'cT1a'
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

        # Visualization section (replacing Analysis Results)
        viz_frame = ttk.Frame(parent)
        viz_frame.grid(row=3, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), pady=(0, 10))

        # Create matplotlib figure
        self.fig = Figure(figsize=(10, 6), dpi=100)
        self.ax1 = self.fig.add_subplot(121)  # Left subplot
        self.ax2 = self.fig.add_subplot(122)  # Right subplot

        # Create canvas
        self.canvas = FigureCanvasTkAgg(self.fig, viz_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Initialize with empty plots
        self.ax1.text(0.5, 0.5, 'Original CT Image\nwith Nodule Contour',
                      ha='center', va='center', transform=self.ax1.transAxes, fontsize=12)
        self.ax1.set_xticks([])
        self.ax1.set_yticks([])

        self.ax2.text(0.5, 0.5, 'CT Image with Component Overlay\n(Red: Solid, Green: GGO)',
                      ha='center', va='center', transform=self.ax2.transAxes, fontsize=12)
        self.ax2.set_xticks([])
        self.ax2.set_yticks([])

        self.canvas.draw()

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
        viz_frame.columnconfigure(0, weight=1)
        viz_frame.rowconfigure(0, weight=1)

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
        self.ax1.set_title(f'Original CT with Nodule Contour (Slice {max_slice})', fontsize=12)
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

            # Calculate basic metrics
            whole_voxel_count = np.sum(nodule_mask_bool)
            solid_voxel_count = np.sum(solid_mask)

            # Calculate masses and volumes
            mass_whole, mean_hu_whole, volume_whole = self.calculate_mass(nodule_hu_values, whole_voxel_count)
            mass_solid, mean_hu_solid, volume_solid = self.calculate_mass(solid_hu_values, solid_voxel_count)

            # Calculate percentages
            sm_percent = (mass_solid / mass_whole * 100) if mass_whole > 0 else 0

            # Update visualization
            self.update_visualization(max_slice, nodule_mask_bool, solid_mask, ggo_mask)

            # Staging recommendation
            new_staging, adjustment = self.adjust_staging(current_staging, sm_percent)
            self.show_staging_recommendation(current_staging, new_staging, adjustment, sm_percent)

        except Exception as e:
            messagebox.showerror("Error", f"Error during calculation: {str(e)}")

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
