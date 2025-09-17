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

def read_gromacs_xvg(filename: str, skip_rows: int = 17) -> np.ndarray:
    """
    Read GROMACS XVG file, handling comments and metadata.
    
    Args:
        filename: Path to XVG file
        skip_rows: Number of rows to skip (default for GROMACS)
    
    Returns:
        Array with data from file
    """
    data_lines = []
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith(('#', '@')):
                data_lines.append(line.strip())
    
    # Parse data
    data = []
    for line in data_lines:
        if line:
            values = [float(x) for x in line.split()]
            data.append(values)
    
    return np.array(data)

def calculate_hydrogen_bonds_gromacs(pdb_file: str, xtc_file: str, 
                                   sel_const: str, sel_pull: str) -> np.ndarray:
    """
    Calculate hydrogen bonds between selections for GROMACS trajectory.
    
    Args:
        pdb_file: PDB structure file
        xtc_file: XTC trajectory file
        sel_const: Selection for constrained part
        sel_pull: Selection for pulled part
    
    Returns:
        Array with [time, hbond_count]
    """
    u = Universe(pdb_file, xtc_file)
    
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
    
    # Create time array (assuming 0.2 ns per frame for GROMACS)
    time_points = np.arange(len(total_hbonds)) * 0.2
    return np.column_stack((time_points, total_hbonds))

def detect_smd_directories(base_dir: str) -> List[str]:
    """Find all SMD simulation directories."""
    return [d for d in os.listdir(base_dir) if d.startswith('SMD_theta') and os.path.isdir(os.path.join(base_dir, d))]

def detect_gromacs_repeats(smd_dir: str) -> int:
    """Detect number of repeats in a GROMACS SMD directory."""
    repeats = 1
    while os.path.exists(os.path.join(smd_dir, f'smd{repeats}_pullf.xvg')):
        repeats += 1
    return repeats - 1  # Subtract 1 because we start from 1

def process_gromacs_smd_simulations(base_dir: str, pdb_file: str, sel_const: str, 
                                   sel_pull: str, calculate_HB: bool = True) -> None:
    """
    Process GROMACS SMD simulations.
    
    Args:
        base_dir: Base directory containing SMD simulations
        pdb_file: PDB structure file
        sel_const: Selection for constrained part
        sel_pull: Selection for pulled part
        calculate_HB: Whether to calculate hydrogen bonds
    """
    # Get reference structure for COM calculation
    u_ref = Universe(os.path.join(base_dir, pdb_file))
    fix_point = u_ref.select_atoms(sel_const).center_of_mass()
    
    # Detect SMD directories
    smd_dirs = detect_smd_directories(base_dir)
    bunch_vectors = []
    
    for smd_dir in smd_dirs:
        full_dir = os.path.join(base_dir, smd_dir)
        n_repeats = detect_gromacs_repeats(full_dir)
        
        print(f"Processing {smd_dir} with {n_repeats} repeats")
        
        for n in range(1, n_repeats + 1):
            pullf_file = os.path.join(full_dir, f'smd{n}_pullf.xvg')
            pullx_file = os.path.join(full_dir, f'smd{n}_pullx.xvg')
            xtc_file = os.path.join(full_dir, f'smd{n}.xtc')
            
            # Calculate hydrogen bonds if requested
            if calculate_HB and os.path.exists(xtc_file):
                print(f'Calculating hydrogen bonds for repeat {n}...')
                hb_data = calculate_hydrogen_bonds_gromacs(
                    os.path.join(base_dir, pdb_file), xtc_file, sel_const, sel_pull
                )
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

def create_gromacs_summary_plots(base_dir: str, red_factor: int = 5) -> None:
    """
    Create summary plots for GROMACS SMD simulations in NAMD format.
    
    Args:
        base_dir: Base directory containing SMD simulations
        red_factor: Data reduction factor for plotting
    """
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
    
    # Create 3-column layout like in NAMD version
    fig, ax = plt.subplots(n_dirs, 3, sharex='col', figsize=(10, 3 * n_dirs))
    
    if n_dirs == 1:
        ax = ax.reshape(1, -1)
    
    # Find global maxima for consistent scaling
    max_force = 0
    max_HB = 0
    
    for dir_name in sorted_dirs:
        n_repeats = detect_gromacs_repeats(os.path.join(base_dir, dir_name))
        for n in range(1, n_repeats + 1):
            force_file = os.path.join(base_dir, dir_name, f'smd{n}_pullf.xvg')
            hb_file = os.path.join(base_dir, dir_name, f'smd_hb_time{n}.dat')
            
            if os.path.exists(force_file):
                force_data = read_gromacs_xvg(force_file)
                if force_data.size > 0:
                    max_force = max(max_force, np.max(force_data[:, 1]))
            
            if os.path.exists(hb_file):
                hb_data = np.loadtxt(hb_file)
                if hb_data.size > 0:
                    max_HB = max(max_HB, np.max(hb_data[:, 1]))
    
    # Add padding to maxima
    max_force = max_force * 1.1 if max_force > 0 else 500
    max_HB = max_HB * 1.1 if max_HB > 0 else 40
    
    # Create plots in NAMD format
    for i, dir_name in enumerate(sorted_dirs):
        n_repeats = detect_gromacs_repeats(os.path.join(base_dir, dir_name))
        
        # Collect data from all repeats
        all_force_time = []
        all_force_dist = []
        all_hb_time = []
        
        for n in range(1, n_repeats + 1):
            force_file = os.path.join(base_dir, dir_name, f'smd{n}_pullf.xvg')
            dist_file = os.path.join(base_dir, dir_name, f'smd{n}_pullx.xvg')
            hb_file = os.path.join(base_dir, dir_name, f'smd_hb_time{n}.dat')
            
            if os.path.exists(force_file):
                force_data = read_gromacs_xvg(force_file)
                all_force_time.append(force_data)
            
            if os.path.exists(dist_file):
                dist_data = read_gromacs_xvg(dist_file)
                all_force_dist.append(dist_data)
            
            if os.path.exists(hb_file):
                hb_data = np.loadtxt(hb_file)
                all_hb_time.append(hb_data)
        
        # Plot force vs time (first column)
        if all_force_time:
            # Find minimum length to avoid dimension mismatches
            min_length = min(len(data) for data in all_force_time)
            times = all_force_time[0][:min_length:red_factor, 0] / 1000  # Convert ps to ns
            forces = np.array([data[:min_length:red_factor, 1] for data in all_force_time])
            
            mean_force = np.mean(forces, axis=0)
            std_force = np.std(forces, axis=0)
            
            ax[i, 0].fill_between(times, mean_force + std_force, mean_force - std_force, 
                                 alpha=0.5, linewidth=1)
            ax[i, 0].plot(times, mean_force, linewidth=2)
            ax[i, 0].set_ylim(0, max_force)
            ax[i, 0].set_xlim(0, 10)
            ax[i, 0].set_ylabel('Force [pN]')
        
        # Plot force vs distance (second column)
        if all_force_dist and all_force_time:
            min_length = min(len(all_force_dist[0]), len(all_force_time[0]))
            distances = np.array([data[:min_length:red_factor, 1] for data in all_force_dist])
            forces = np.array([data[:min_length:red_factor, 1] for data in all_force_time])
            times = all_force_time[0][:min_length:red_factor, 0] / 1000  # For coloring
            
            mean_distance = np.mean(distances, axis=0)
            mean_force = np.mean(forces, axis=0)
            
            scatter = ax[i, 1].scatter(mean_distance, mean_force, s=2, c=times, cmap='viridis')
            ax[i, 1].set_ylim(0, max_force)
            ax[i, 1].set_ylabel('Force [pN]')
        
        # Plot hydrogen bonds vs time (third column)
        if all_hb_time:
            min_length = min(len(data) for data in all_hb_time)
            times = all_hb_time[0][:min_length, 0]
            hb_counts = np.array([data[:min_length, 1] for data in all_hb_time])
            
            mean_hb = np.mean(hb_counts, axis=0)
            std_hb = np.std(hb_counts, axis=0)
            
            ax[i, 2].fill_between(times, mean_hb + std_hb, mean_hb - std_hb, 
                                 alpha=0.5, linewidth=1)
            ax[i, 2].plot(times, mean_hb, linewidth=2)
            ax[i, 2].set_ylim(0, max_HB)
            ax[i, 2].set_xlim(0, 10)
            ax[i, 2].set_ylabel('HB number')
        
        # Set title with direction info
        parts = dir_name.split('_')
        theta = parts[2]
        phi = parts[4]
        ax[i, 0].set_title(f'θ={theta}, φ={phi}', fontsize=10)
    
    # Set common labels
    ax[-1, 0].set_xlabel('Time [ns]')
    ax[-1, 1].set_xlabel(r'Distance [$\AA$]')
    ax[-1, 2].set_xlabel('Time [ns]')
    
    # Add colorbar for force vs distance plot
    if any(all_force_dist for _ in sorted_dirs):
        cbar = plt.colorbar(scatter, ax=ax[:, 1].ravel().tolist())
        cbar.set_label('Time [ns]')
    
    fig.tight_layout()
    plt.savefig(os.path.join(base_dir, 'all_smd_results.png'), dpi=600, bbox_inches='tight')
    plt.show()

def main():
    """Main function with argument parsing for GROMACS SMD analysis."""
    parser = argparse.ArgumentParser(description='Analyze GROMACS SMD simulation data')
    parser.add_argument('directory', help='SMD output directory')
    parser.add_argument('pdb_file', help='PDB structure file')
    parser.add_argument('sel_const', help='Selection for constrained part')
    parser.add_argument('sel_pull', help='Selection for pulled part')
    parser.add_argument('--no-hb', action='store_true', help='Skip hydrogen bond calculation')
    parser.add_argument('--red-factor', type=int, default=5, help='Data reduction factor for plotting')
    parser.add_argument('--plot-only', action='store_true', help='Only create plots, skip data processing')
    
    args = parser.parse_args()
    
    if not args.plot_only:
        process_gromacs_smd_simulations(
            args.directory,
            args.pdb_file,
            args.sel_const,
            args.sel_pull,
            calculate_HB=not args.no_hb
        )
    
    create_gromacs_summary_plots(args.directory, args.red_factor)

if __name__ == '__main__':
    main()
