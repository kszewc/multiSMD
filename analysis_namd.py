import numpy as np
from MDAnalysis import Universe
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import os
import sys
import matplotlib.pyplot as plt
import argparse
from typing import List, Tuple, Optional

def distance(a: np.ndarray, b: np.ndarray) -> float:
    """Calculate Euclidean distance between two points."""
    return np.linalg.norm(a - b)

def extract_smd_forces(log_file: str, fix_point: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract force vs time and force vs distance data from SMD log file.
    
    Returns:
        force_time: array with [time, force]
        force_dist: array with [distance, force]
    """
    force_time_data = []
    force_dist_data = []
    
    with open(log_file, 'r') as f:
        for line in f:
            if line.startswith('SMD  ') and len(line.split()) >= 8:
                parts = line.split()
                t = float(parts[1])  # timestep
                r = np.array([float(parts[2]), float(parts[3]), float(parts[4])])  # COM position
                f_vec = np.array([float(parts[5]), float(parts[6]), float(parts[7])])  # Pulling force
                
                force_magnitude = np.linalg.norm(f_vec)
                dist = distance(r, fix_point)
                
                force_time_data.append([t, force_magnitude])
                force_dist_data.append([dist, force_magnitude])
    
    return np.array(force_time_data), np.array(force_dist_data)

def calculate_hydrogen_bonds(psf_file: str, dcd_file: str, sel_const: str, sel_pull: str) -> np.ndarray:
    """Calculate hydrogen bonds between selections over trajectory."""
    u = Universe(psf_file, dcd_file)
    
    # Hydrogen bonds in both directions
    hbonds_forward = HBA(
        universe=u,
        hydrogens_sel="protein and name H*",
        donors_sel=sel_const,
        acceptors_sel=sel_pull
    )
    hbonds_forward.run()
    
    hbonds_reverse = HBA(
        universe=u,
        hydrogens_sel="protein and name H*",
        donors_sel=sel_pull,
        acceptors_sel=sel_const
    )
    hbonds_reverse.run()
    
    # Combine both directions
    total_hbonds = hbonds_forward.count_by_time() + hbonds_reverse.count_by_time()
    
    # Create time array (assuming 0.05 ns per frame)
    time_points = np.arange(len(total_hbonds)) * 0.05
    return np.column_stack((time_points, total_hbonds))

def detect_smd_directories(base_dir: str) -> List[str]:
    """Find all SMD simulation directories."""
    return [d for d in os.listdir(base_dir) if d.startswith('SMD_theta') and os.path.isdir(os.path.join(base_dir, d))]

def detect_repeats(smd_dir: str) -> int:
    """Detect number of repeats in an SMD directory."""
    repeats = 0
    while os.path.exists(os.path.join(smd_dir, f'mdrun{repeats}.log')):
        repeats += 1
    return repeats

def process_smd_simulations(base_dir: str, sel_const: str, sel_pull: str, 
                           extract_forces: bool = True, calculate_HB: bool = True,
                           red_factor: int = 5) -> None:
    """Main function to process all SMD simulations."""
    
    # Find reference structure for COM calculation
    u_ref = Universe(os.path.join(base_dir, 'SMD_constraints.pdb'))
    fix_point = u_ref.select_atoms(sel_const).center_of_mass()
    
    # Find PSF file
    psf_files = [f for f in os.listdir(base_dir) if f.endswith('.psf')]
    if not psf_files:
        raise FileNotFoundError("No PSF file found in base directory")
    psf_file = os.path.join(base_dir, psf_files[0])
    
    # Detect SMD directories and process them
    smd_dirs = detect_smd_directories(base_dir)
    bunch_vectors = []
    
    for smd_dir in smd_dirs:
        full_dir = os.path.join(base_dir, smd_dir)
        n_repeats = detect_repeats(full_dir)
        
        print(f"Processing {smd_dir} with {n_repeats} repeats")
        
        for n in range(n_repeats):
            log_file = os.path.join(full_dir, f'mdrun{n}.log')
            dcd_file = os.path.join(full_dir, f'md{n}.dcd')
            
            # Extract forces
            if extract_forces and os.path.exists(log_file):
                force_time, force_dist = extract_smd_forces(log_file, fix_point)
                
                np.savetxt(os.path.join(full_dir, f'smd_force_time{n}.dat'), 
                          force_time, header='time [fs] vs force of pulling')
                np.savetxt(os.path.join(full_dir, f'smd_force_dist{n}.dat'), 
                          force_dist, header='distance between COMs vs force of pulling')
            
            # Calculate hydrogen bonds
            if calculate_HB and os.path.exists(dcd_file):
                print('Calculating hydrogen bonds...')
                hb_data = calculate_hydrogen_bonds(psf_file, dcd_file, sel_const, sel_pull)
                np.savetxt(os.path.join(full_dir, f'smd_hb_time{n}.dat'), 
                          hb_data, header='time [ns] vs number of hydrogen bonds between domains')
        
        # Extract pulling vector from directory name
        parts = smd_dir.split('_')
        theta = float(parts[2])
        phi = float(parts[4])
        
        x = np.cos(np.deg2rad(phi)) * np.sin(np.deg2rad(theta))
        y = np.sin(np.deg2rad(phi)) * np.sin(np.deg2rad(theta))
        z = np.cos(np.deg2rad(theta))
        
        bunch_vectors.append([x, y, z])
    
    # Save pulling vectors
    np.savetxt(os.path.join(base_dir, 'bunch_of_vectors.dat'), bunch_vectors)

def create_summary_plots(base_dir: str, red_factor: int = 5) -> None:
    """Create summary plots for all SMD simulations."""
    smd_dirs = detect_smd_directories(base_dir)
    
    # Default sorting for common pulling directions
    default_order = [
        'SMD_theta_0_phi_0', 'SMD_theta_45_phi_0', 'SMD_theta_45_phi_90',
        'SMD_theta_45_phi_180', 'SMD_theta_45_phi_270', 'SMD_theta_90_phi_0',
        'SMD_theta_90_phi_90', 'SMD_theta_90_phi_180', 'SMD_theta_90_phi_270'
    ]
    
    # Use default order if available, otherwise sort alphabetically
    sorted_dirs = [d for d in default_order if d in smd_dirs]
    if not sorted_dirs:
        sorted_dirs = sorted(smd_dirs)
    
    n_dirs = len(sorted_dirs)
    fig, ax = plt.subplots(n_dirs, 3, sharex='col', figsize=(10, 3 * n_dirs))
    
    if n_dirs == 1:
        ax = ax.reshape(1, -1)
    
    # Find global maxima for consistent scaling
    max_force = 0
    max_HB = 0
    
    for dir_name in sorted_dirs:
        n_repeats = detect_repeats(os.path.join(base_dir, dir_name))
        for n in range(n_repeats):
            force_file = os.path.join(base_dir, dir_name, f'smd_force_time{n}.dat')
            hb_file = os.path.join(base_dir, dir_name, f'smd_hb_time{n}.dat')
            
            if os.path.exists(force_file):
                force_data = np.loadtxt(force_file)
                max_force = max(max_force, np.max(force_data[:, 1]) if force_data.size > 0 else 0)
            
            if os.path.exists(hb_file):
                hb_data = np.loadtxt(hb_file)
                max_HB = max(max_HB, np.max(hb_data[:, 1]) if hb_data.size > 0 else 0)
    
    # Add some padding to maxima
    max_force *= 1.1
    max_HB *= 1.1
    
    # Create plots
    for i, dir_name in enumerate(sorted_dirs):
        n_repeats = detect_repeats(os.path.join(base_dir, dir_name))
        
        # Collect data from all repeats
        all_force_time = []
        all_force_dist = []
        all_hb_time = []
        
        for n in range(n_repeats):
            force_time_file = os.path.join(base_dir, dir_name, f'smd_force_time{n}.dat')
            force_dist_file = os.path.join(base_dir, dir_name, f'smd_force_dist{n}.dat')
            hb_time_file = os.path.join(base_dir, dir_name, f'smd_hb_time{n}.dat')
            
            if os.path.exists(force_time_file):
                data = np.loadtxt(force_time_file)
                all_force_time.append(data)
            
            if os.path.exists(force_dist_file):
                data = np.loadtxt(force_dist_file)
                all_force_dist.append(data)
            
            if os.path.exists(hb_time_file):
                data = np.loadtxt(hb_time_file)
                all_hb_time.append(data)
        
        # Plot force vs time
        if all_force_time:
            times = all_force_time[0][::red_factor, 0] / 500000  # Convert fs to ns
            forces = np.array([data[::red_factor, 1] for data in all_force_time])
            
            mean_force = np.mean(forces, axis=0)
            std_force = np.std(forces, axis=0)
            
            ax[i, 0].fill_between(times, mean_force + std_force, mean_force - std_force, 
                                 alpha=0.3, label='±1 SD')
            ax[i, 0].plot(times, mean_force, linewidth=2, label='Mean')
            ax[i, 0].set_ylabel('Force [pN]')
            ax[i, 0].set_ylim(0, max_force)
        
        # Plot force vs distance
        if all_force_dist:
            distances = np.array([data[::red_factor, 0] for data in all_force_dist])
            forces = np.array([data[::red_factor, 1] for data in all_force_dist])
            times = all_force_time[0][::red_factor, 0] if all_force_time else None
            
            mean_force = np.mean(forces, axis=0)
            
            if times is not None:
                scatter = ax[i, 1].scatter(np.mean(distances, axis=0), mean_force, 
                                          c=times/500000, s=10, cmap='viridis')
            else:
                ax[i, 1].plot(np.mean(distances, axis=0), mean_force, linewidth=2)
            
            ax[i, 1].set_ylabel('Force [pN]')
            ax[i, 1].set_ylim(0, max_force)
        
        # Plot hydrogen bonds vs time
        if all_hb_time:
            times = all_hb_time[0][:, 0]
            hb_counts = np.array([data[:, 1] for data in all_hb_time])
            
            mean_hb = np.mean(hb_counts, axis=0)
            std_hb = np.std(hb_counts, axis=0)
            
            ax[i, 2].fill_between(times, mean_hb + std_hb, mean_hb - std_hb, alpha=0.3)
            ax[i, 2].plot(times, mean_hb, linewidth=2)
            ax[i, 2].set_ylabel('HB number')
            ax[i, 2].set_ylim(0, max_HB)
        
        # Set title with direction info
        ax[i, 0].set_title(dir_name, fontsize=10)
    
    # Set common labels
    for i in range(n_dirs):
        ax[i, 0].set_ylabel('Force [pN]')
        ax[i, 1].set_ylabel('Force [pN]')
        ax[i, 2].set_ylabel('HB number')
    
    ax[-1, 0].set_xlabel('Time [ns]')
    ax[-1, 1].set_xlabel(r'Distance [$\AA$]')
    ax[-1, 2].set_xlabel('Time [ns]')
    
    # Add colorbar for force vs distance plot
    #if any(all_force_dist for _ in sorted_dirs):
    #    cbar = plt.colorbar(scatter, ax=ax[:, 1].ravel().tolist())
    #    cbar.set_label('Time [ns]')
    fig.subplots_adjust()
    fig.set_tight_layout(True)
    plt.savefig(os.path.join(base_dir, 'all_smd_results.png'), dpi=600, bbox_inches='tight')
    plt.show()

def main():
    """Main function with argument parsing."""
    parser = argparse.ArgumentParser(description='Analyze SMD simulation data')
    parser.add_argument('directory', help='SMD output directory')
    parser.add_argument('sel_const', help='Selection for constrained part')
    parser.add_argument('sel_pull', help='Selection for pulled part')
    parser.add_argument('--repeats', type=int, help='Number of repeats (auto-detected if not specified)')
    parser.add_argument('--no-forces', action='store_true', help='Skip force extraction')
    parser.add_argument('--no-hb', action='store_true', help='Skip hydrogen bond calculation')
    parser.add_argument('--red-factor', type=int, default=5, help='Data reduction factor for plotting')
    parser.add_argument('--plot-only', action='store_true', help='Only create plots, skip data processing')
    
    args = parser.parse_args()
    
    if not args.plot_only:
        process_smd_simulations(
            args.directory,
            args.sel_const,
            args.sel_pull,
            extract_forces=not args.no_forces,
            calculate_HB=not args.no_hb,
            red_factor=args.red_factor
        )
    
    create_summary_plots(args.directory, args.red_factor)

if __name__ == '__main__':
    main()
