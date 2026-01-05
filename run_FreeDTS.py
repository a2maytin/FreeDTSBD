#!/usr/bin/env python3
"""
Run a series of FreeDTS simulations with varying parameters.

Main focus: Varying the number of threads (nt_list).
Also supports varying area, volume, curvature, and other parameters.
"""

import os
import subprocess
from pathlib import Path


def generate_input_files(template_input, input_folder, traj_tag, nt_list,
                         area_traj=None, gammaV_list=None, 
                         curvature_list=None, scurvature_list=None,
                         sim_steps=None):
    """
    Generate input.dts files with modified parameters.
    Creates one input file per thread count in nt_list.
    
    Parameters:
    -----------
    template_input : str or Path
        Path to template input.dts file
    input_folder : str or Path
        Folder where input files will be created
    traj_tag : str
        Tag to include in filenames
    nt_list : list
        List of thread counts (determines number of input files)
    area_traj : list, optional
        List of area target values (one per thread count)
    gammaV_list : list, optional
        List of volume target values (one per thread count)
    curvature_list : list, optional
        List of global curvature values (one per thread count)
    scurvature_list : list, optional
        List of spontaneous curvature values (one per thread count)
    sim_steps : int, optional
        Number of simulation steps (if None, uses template value)
    """
    print("Generating input files...")
    input_folder = Path(input_folder)
    input_folder.mkdir(parents=True, exist_ok=True)
    template_input = Path(template_input)
    
    # Read template
    with open(template_input, "r") as f:
        template_lines = f.readlines()
    
    # Generate one input file per thread count
    for i, threads in enumerate(nt_list):
        # Get parameter values for this simulation (use first value if list is shorter)
        area = area_traj[i] if area_traj and i < len(area_traj) else None
        gammaV = gammaV_list[i] if gammaV_list and i < len(gammaV_list) else None
        curvature = curvature_list[i] if curvature_list and i < len(curvature_list) else None
        scurvature = scurvature_list[i] if scurvature_list and i < len(scurvature_list) else None
        
        # Create filename using thread count (matching old script pattern)
        filename_append = threads
        new_input_path = input_folder / f"input_{traj_tag}_{filename_append}.dts"
        
        with open(new_input_path, "w") as f:
            for line in template_lines:
                modified = False
                
                # Modify VolumeCoupling
                if line.startswith("VolumeCoupling") and gammaV is not None:
                    # Format: VolumeCoupling = SecondOrder <param1> <param2> <target>
                    # Try to preserve existing param1 and param2 if present
                    parts = line.split()
                    if len(parts) >= 5 and parts[2] == "SecondOrder":
                        param1 = parts[3] if len(parts) > 3 else "0"
                        param2 = parts[4] if len(parts) > 4 else "10000"
                        f.write(f"VolumeCoupling = SecondOrder {param1} {param2} {gammaV}\n")
                    else:
                        f.write(f"VolumeCoupling = SecondOrder 0 10000 {gammaV}\n")
                    modified = True
                
                # Modify TotalAreaCoupling
                elif line.startswith("TotalAreaCoupling") and area is not None:
                    # Format: TotalAreaCoupling = HarmonicPotential <stiffness> <target>
                    parts = line.split()
                    if len(parts) >= 4 and parts[2] == "HarmonicPotential":
                        stiffness = parts[3] if len(parts) > 3 else "1000"
                        f.write(f"TotalAreaCoupling = HarmonicPotential {stiffness} {area}\n")
                    else:
                        f.write(f"TotalAreaCoupling = HarmonicPotential 1000 {area}\n")
                    modified = True
                
                # Modify GlobalCurvatureCoupling
                elif line.startswith("GlobalCurvatureCoupling") and curvature is not None:
                    # Format: GlobalCurvatureCoupling = HarmonicPotential <stiffness> <target>
                    parts = line.split()
                    if len(parts) >= 4 and parts[2] == "HarmonicPotential":
                        stiffness = parts[3] if len(parts) > 3 else "150"
                        f.write(f"GlobalCurvatureCoupling = HarmonicPotential {stiffness} {curvature}\n")
                    else:
                        f.write(f"GlobalCurvatureCoupling = HarmonicPotential 150 {curvature}\n")
                    modified = True
                
                # Modify Kappa (spontaneous curvature)
                elif line.startswith("Kappa") and scurvature is not None:
                    # Format: Kappa = <kappa> <param2> <spontaneous_curvature>
                    parts = line.split()
                    if len(parts) >= 4:
                        kappa = parts[2] if len(parts) > 2 else "60"
                        param2 = parts[3] if len(parts) > 3 else "0"
                        f.write(f"Kappa = {kappa} {param2} {scurvature}\n")
                    else:
                        f.write(f"Kappa = 60 0 {scurvature}\n")
                    modified = True
                
                # Modify Set_Steps
                elif line.startswith("Set_Steps") and sim_steps is not None:
                    f.write(f"Set_Steps = 1 {sim_steps}\n")
                    modified = True
                
                # Keep original line if not modified
                if not modified:
                    f.write(line)
        
        print(f"Created {new_input_path}")


def run_simulations(template_input, output_folder, top_file, topq_file,
                    dts_executable, traj_tag, nt_list, 
                    area_traj=None, gammaV_list=None,
                    curvature_list=None, scurvature_list=None,
                    sim_steps=None, replicates=1):
    """
    Run FreeDTS simulations with varying thread counts.
    Creates one folder per run containing input.dts and all outputs.
    
    Parameters:
    -----------
    template_input : str or Path
        Path to template input.dts file
    output_folder : str or Path
        Base folder for simulation run folders
    top_file : str or Path
        Path to topology.top file
    topq_file : str or Path, optional
        Path to topology.q file (if needed)
    dts_executable : str or Path
        Path to DTS executable
    traj_tag : str
        Tag used in filenames
    nt_list : list
        List of thread counts (one per simulation)
    area_traj : list, optional
        List of area values (one per thread count)
    gammaV_list : list, optional
        List of volume values (one per thread count)
    curvature_list : list, optional
        List of global curvature values (one per thread count)
    scurvature_list : list, optional
        List of spontaneous curvature values (one per thread count)
    sim_steps : int, optional
        Number of simulation steps
    replicates : int
        Number of replicates to run for each configuration
    """
    print("Running simulations...")
    # Convert all paths to absolute paths before changing directories
    template_input = Path(template_input).resolve()
    output_folder = Path(output_folder).resolve()
    top_file = Path(top_file).resolve() if top_file else None
    topq_file = Path(topq_file).resolve() if topq_file else None
    dts_executable = Path(dts_executable).resolve()
    
    # Check that required files exist
    if not template_input.exists():
        print(f"ERROR: Template input file not found: {template_input}")
        return []
    
    if not dts_executable.exists():
        print(f"ERROR: DTS executable not found: {dts_executable}")
        return []
    
    if top_file and not top_file.exists():
        print(f"WARNING: Topology file not found: {top_file}")
        print("  Continuing anyway, but simulation may fail if topology is required.")
    
    output_folder.mkdir(parents=True, exist_ok=True)
    
    original_dir = os.getcwd()
    os.chdir(output_folder)
    
    # Read template
    with open(template_input, "r") as f:
        template_lines = f.readlines()
    
    num_sims = len(nt_list)
    
    try:
        for j in range(1, replicates + 1):
            for idx, threads in enumerate(nt_list):
                # Create run directory (one folder per run)
                filename_append = threads
                if replicates > 1:
                    run_dir = Path(f"run_{traj_tag}_{filename_append}_{j}")
                else:
                    run_dir = Path(f"run_{traj_tag}_{filename_append}")
                
                run_dir.mkdir(exist_ok=True)
                os.chdir(run_dir)
                
                # Get parameter values for this simulation
                area = area_traj[idx] if area_traj and idx < len(area_traj) else None
                gammaV = gammaV_list[idx] if gammaV_list and idx < len(gammaV_list) else None
                curvature = curvature_list[idx] if curvature_list and idx < len(curvature_list) else None
                scurvature = scurvature_list[idx] if scurvature_list and idx < len(scurvature_list) else None
                
                # Generate input.dts file directly in this run folder
                local_input_file = Path("input.dts")
                with open(local_input_file, "w") as f:
                    for line in template_lines:
                        modified = False
                        
                        # Modify VolumeCoupling
                        if line.startswith("VolumeCoupling") and gammaV is not None:
                            parts = line.split()
                            # Handle both "VolumeCoupling = No" and "VolumeCoupling = SecondOrder ..."
                            if len(parts) >= 5 and parts[2] == "SecondOrder":
                                param1 = parts[3] if len(parts) > 3 else "0"
                                param2 = parts[4] if len(parts) > 4 else "10000"
                                f.write(f"VolumeCoupling = SecondOrder {param1} {param2} {gammaV}\n")
                            else:
                                # Replace "No" or any other value with SecondOrder
                                f.write(f"VolumeCoupling = SecondOrder 0 10000 {gammaV}\n")
                            modified = True
                        
                        # Modify TotalAreaCoupling
                        elif line.startswith("TotalAreaCoupling") and area is not None:
                            parts = line.split()
                            # Handle both "TotalAreaCoupling = No" and "TotalAreaCoupling = HarmonicPotential ..."
                            if len(parts) >= 4 and parts[2] == "HarmonicPotential":
                                stiffness = parts[3] if len(parts) > 3 else "1000"
                                f.write(f"TotalAreaCoupling = HarmonicPotential {stiffness} {area}\n")
                            else:
                                # Replace "No" or any other value with HarmonicPotential
                                f.write(f"TotalAreaCoupling = HarmonicPotential 1000 {area}\n")
                            modified = True
                        
                        # Modify GlobalCurvatureCoupling
                        elif line.startswith("GlobalCurvatureCoupling") and curvature is not None:
                            parts = line.split()
                            # Handle both "GlobalCurvatureCoupling = No" and "GlobalCurvatureCoupling = HarmonicPotential ..."
                            if len(parts) >= 4 and parts[2] == "HarmonicPotential":
                                stiffness = parts[3] if len(parts) > 3 else "150"
                                f.write(f"GlobalCurvatureCoupling = HarmonicPotential {stiffness} {curvature}\n")
                            else:
                                # Replace "No" or any other value with HarmonicPotential
                                f.write(f"GlobalCurvatureCoupling = HarmonicPotential 150 {curvature}\n")
                            modified = True
                        
                        # Modify Kappa (spontaneous curvature)
                        elif line.startswith("Kappa") and scurvature is not None:
                            parts = line.split()
                            if len(parts) >= 4:
                                kappa = parts[2] if len(parts) > 2 else "60"
                                param2 = parts[3] if len(parts) > 3 else "0"
                                f.write(f"Kappa = {kappa} {param2} {scurvature}\n")
                            else:
                                f.write(f"Kappa = 60 0 {scurvature}\n")
                            modified = True
                        
                        # Modify Set_Steps
                        elif line.startswith("Set_Steps") and sim_steps is not None:
                            f.write(f"Set_Steps = 1 {sim_steps}\n")
                            modified = True
                        
                        # Keep original line if not modified
                        if not modified:
                            f.write(line)
                
                print(f"Created input.dts in {run_dir}")
                
                # Copy topology files to run directory
                files_copied = []
                if top_file and top_file.exists():
                    # Copy as topology.top (required by DTS)
                    result = subprocess.run(["cp", str(top_file), "topology.top"], check=False, capture_output=True)
                    if result.returncode == 0:
                        files_copied.append("topology.top")
                    # Also copy as top.top for compatibility
                    subprocess.run(["cp", str(top_file), "top.top"], check=False)
                else:
                    print(f"WARNING: Topology file not found: {top_file}")
                    print("  Simulation may fail without topology file.")
                
                # Copy optional topology-related files
                if topq_file and topq_file.exists():
                    result = subprocess.run(["cp", str(topq_file), "topol.q"], check=False, capture_output=True)
                    if result.returncode == 0:
                        files_copied.append("topol.q")
                
                # Check for dna.q in the same directory as top_file
                if top_file and top_file.exists():
                    top_dir = top_file.parent
                    dna_q = top_dir / "dna.q"
                    if dna_q.exists():
                        result = subprocess.run(["cp", str(dna_q), "dna.q"], check=False, capture_output=True)
                        if result.returncode == 0:
                            files_copied.append("dna.q")
                
                if files_copied:
                    print(f"  Copied files to run directory: {', '.join(files_copied)}")
                else:
                    print(f"  WARNING: No additional files were copied to run directory")
                
                # Run simulation
                print(f"\n{'='*60}")
                print(f"Running simulation: threads={threads}, replicate={j}/{replicates}")
                print(f"Run folder: {run_dir}")
                print(f"Current directory: {os.getcwd()}")
                print(f"{'='*60}")
                
                # Build command - use relative paths since we're in the run directory
                # This matches how the user runs it manually
                cmd_str = f"{dts_executable} -in input.dts -top topology.top -nt {threads}"
                
                print(f"Command: {cmd_str}")
                print(f"Working directory: {os.getcwd()}")
                print(f"Input file exists: {Path('input.dts').exists()}")
                print(f"Topology file exists: {Path('topology.top').exists()}")
                
                # Run using shell=True to match manual execution
                result = subprocess.run(cmd_str, shell=True, check=False)
                
                if result.returncode == 0:
                    print(f"✓ Simulation completed successfully")
                else:
                    print(f"✗ Simulation exited with code {result.returncode}")
                    # Check for log files that might contain error information
                    log_file = Path("dts.log")
                    if log_file.exists():
                        print(f"  Checking dts.log for errors...")
                        try:
                            with open(log_file, 'r') as f:
                                lines = f.readlines()
                                # Print last 20 lines of log file
                                print(f"  Last 20 lines of dts.log:")
                                for line in lines[-20:]:
                                    print(f"    {line.rstrip()}")
                        except Exception as e:
                            print(f"  Could not read log file: {e}")
                    else:
                        print(f"  No dts.log file found (simulation may have crashed before creating log)")
                
                os.chdir("..")
    
    finally:
        os.chdir(original_dir)
    
    print("\nFinished running FreeDTS.")


if __name__ == "__main__":
    # ========================================================================
    # USER CONFIGURATION: Modify these parameters for your simulations
    # ========================================================================
    
    # Paths
    dts_executable = Path("../version_2.1/DTS")  # Path to DTS executable
    input_template = Path("./input_template.dts")  # Template input file
    top_file = Path("./topology.top")  # Topology file
    # topq_file is optional - script will look for topol.q and dna.q in same dir as top_file
    topq_file = None  # Or specify explicitly: Path("toy_system/topol.q")
    output_folder = Path(".")  # Where to store simulation run folders
    
    # Simulation parameters
    sim_steps = 1000000  # Number of simulation steps
    replicates = 1  # Number of replicates per configuration
    traj_tag = "debug"  # Tag for filenames
    
    # Parameter lists (main focus: nt_list for varying threads)
    nt_list = [4, 15, 30]  # Number of threads - THIS IS THE MAIN PARAMETER TO VARY
    
    # Optional: Other parameters (set to None or same values if not varying)
    area_traj = [0.34, 0.34, 0.34]  # Area targets
    gammaV_list = [1.0, 1.0, 1.0]  # Volume targets
    curvature_list = [0.0, 0.0, 0.0]  # Global curvature c_G0
    scurvature_list = [0.0, 0.0, 0.0]  # Spontaneous curvature c_0
    
    # ========================================================================
    # Generate input files and run simulations
    # ========================================================================
    
    print("Starting program...")
    print(f"Number of thread configurations: {len(nt_list)}")
    print(f"Thread counts: {nt_list}")
    print(f"Replicates per configuration: {replicates}")
    
    # Run simulations (input files are generated in each run folder)
    run_simulations(
        template_input=input_template,
        output_folder=output_folder,
        top_file=top_file,
        topq_file=topq_file,
        dts_executable=dts_executable,
        traj_tag=traj_tag,
        nt_list=nt_list,
        area_traj=area_traj,
        gammaV_list=gammaV_list,
        curvature_list=curvature_list,
        scurvature_list=scurvature_list,
        sim_steps=sim_steps,
        replicates=replicates
    )
    
    print("\nAll simulations completed!")
