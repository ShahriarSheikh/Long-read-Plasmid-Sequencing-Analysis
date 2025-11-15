#!/usr/bin/env python3
"""
Plasmid Clone Validation Workflow Runner
Interactive Python script for running wf-clone-validation workflow
VERSION WITH PROGRESS SPINNER
"""

import os
import sys
import subprocess
import shutil
import tempfile
import glob
import threading
import time
from pathlib import Path
from typing import List, Dict, Tuple

# Global variable for assembly tool selection
ASSEMBLY_TOOL = "flye"

class ProgressSpinner:
    """A simple spinning progress indicator"""
    def __init__(self, message: str = "Working"):
        self.message = message
        self.spinning = False
        self.spinner_chars = "|/-\\"
        self.spinner_thread = None
        
    def start(self):
        """Start the spinner"""
        self.spinning = True
        self.spinner_thread = threading.Thread(target=self._spin)
        self.spinner_thread.daemon = True
        self.spinner_thread.start()
        
    def stop(self):
        """Stop the spinner"""
        self.spinning = False
        if self.spinner_thread:
            self.spinner_thread.join(timeout=1)
        # Clear the spinner line
        sys.stdout.write('\r' + ' ' * (len(self.message) + 10) + '\r')
        sys.stdout.flush()
        
    def _spin(self):
        """Internal spinning method"""
        i = 0
        while self.spinning:
            char = self.spinner_chars[i % len(self.spinner_chars)]
            elapsed_time = int(time.time() - self.start_time) if hasattr(self, 'start_time') else 0
            minutes, seconds = divmod(elapsed_time, 60)
            time_str = f"{minutes:02d}:{seconds:02d}"
            
            sys.stdout.write(f'\r{Colors.BLUE}{char} {self.message} ({time_str}){Colors.NC}')
            sys.stdout.flush()
            time.sleep(0.1)
            i += 1
            
    def __enter__(self):
        """Context manager entry"""
        self.start_time = time.time()
        self.start()
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit"""
        self.stop()

class Colors:
    """ANSI color codes for terminal output"""
    RED = '\033[0;31m'
    GREEN = '\033[0;32m'
    YELLOW = '\033[1;33m'
    BLUE = '\033[0;34m'
    BOLD = '\033[1m'
    NC = '\033[0m'  # No Color

def print_colored(color: str, message: str) -> None:
    """Print colored message to terminal"""
    print(f"{color}{message}{Colors.NC}")

def convert_windows_path(path: str) -> str:
    """Convert Windows path to WSL/Linux path format"""
    path = path.strip().strip('"').strip("'")
    
    # Check if it's a Windows path (contains drive letter)
    if len(path) >= 2 and path[1] == ':' and path[0].isalpha():
        drive = path[0].lower()
        rest = path[2:].replace('\\', '/')
        if rest.startswith('/'):
            rest = rest[1:]
        return f"/mnt/{drive}/{rest}"
    
    return path

def validate_directory(directory: str, description: str) -> bool:
    """Validate that directory exists"""
    if not os.path.isdir(directory):
        print_colored(Colors.RED, f"Error: {description} directory does not exist: {directory}")
        return False
    return True

def get_user_input(prompt: str, default: str = None) -> str:
    """Get user input with optional default value"""
    if default:
        user_input = input(f"{prompt} [{default}]: ").strip()
        return user_input if user_input else default
    else:
        return input(f"{prompt}: ").strip()

def get_barcodes(input_dir: str) -> List[str]:
    """Get list of barcodes (subdirectories with fastq files)"""
    barcodes = []
    
    try:
        for item in os.listdir(input_dir):
            item_path = os.path.join(input_dir, item)
            if os.path.isdir(item_path):
                # Check if directory contains fastq files
                fastq_files = glob.glob(os.path.join(item_path, "*.fastq*"))
                if fastq_files:
                    barcodes.append(item)
    except PermissionError:
        print_colored(Colors.RED, f"Permission denied accessing directory: {input_dir}")
        return []
    
    return sorted(barcodes)

def create_complete_sample_sheet(barcodes: List[str], approx_size: int) -> str:
    """Create complete sample sheet for all barcodes"""
    temp_dir = tempfile.gettempdir()
    sample_sheet_path = os.path.join(temp_dir, "complete_sample_sheet.csv")
    
    with open(sample_sheet_path, 'w') as f:
        f.write("alias,barcode,type,approx_size\n")
        for barcode in barcodes:
            # Create alias that doesn't start with 'barcode'
            if barcode.startswith('barcode'):
                alias = f"sample_{barcode[7:]}"  # Remove 'barcode' prefix and add 'sample_'
            else:
                alias = f"sample_{barcode}"
            f.write(f"{alias},{barcode},test_sample,{approx_size}\n")
    
    return sample_sheet_path

def create_sample_sheet(barcode: str, approx_size: int) -> str:
    """Create temporary sample sheet for a single barcode"""
    temp_dir = tempfile.gettempdir()
    sample_sheet_path = os.path.join(temp_dir, f"sample_sheet_{barcode}.csv")
    
    # Create alias that doesn't start with 'barcode'
    if barcode.startswith('barcode'):
        alias = f"sample_{barcode[7:]}"  # Remove 'barcode' prefix and add 'sample_'
    else:
        alias = f"sample_{barcode}"
    
    with open(sample_sheet_path, 'w') as f:
        f.write("alias,barcode,type,approx_size\n")
        f.write(f"{alias},{barcode},test_sample,{approx_size}\n")
    
    return sample_sheet_path

def run_quick_analysis_all_barcodes(input_dir: str, output_dir: str, barcodes: List[str]) -> bool:
    """Run quick analysis on all barcodes simultaneously"""
    print_colored(Colors.BLUE, f"Running quick analysis on ALL barcodes: {', '.join(barcodes)}")
    print_colored(Colors.YELLOW, f"Parameters: size=7000, quality=9, coverage=2, assembler={ASSEMBLY_TOOL}")
    
    if ASSEMBLY_TOOL == "canu":
        print_colored(Colors.YELLOW, "Note: Canu is much slower than Flye - this may take 30+ minutes")
    
    # Debug: Check if all barcode directories exist and contain fastq files
    print_colored(Colors.BLUE, "Checking barcode directories:")
    for barcode in barcodes:
        barcode_path = os.path.join(input_dir, barcode)
        if os.path.exists(barcode_path):
            fastq_files = glob.glob(os.path.join(barcode_path, "*.fastq*"))
            print_colored(Colors.GREEN, f"  {barcode}: {len(fastq_files)} fastq files found")
        else:
            print_colored(Colors.RED, f"  {barcode}: directory not found!")
    
    # Create complete sample sheet for all barcodes
    sample_sheet_path = create_complete_sample_sheet(barcodes, 7000)
    
    # Debug: show the sample sheet content
    print_colored(Colors.BLUE, f"Created sample sheet: {sample_sheet_path}")
    try:
        with open(sample_sheet_path, 'r') as f:
            content = f.read()
            print_colored(Colors.YELLOW, "Sample sheet content:")
            for line_num, line in enumerate(content.strip().split('\n'), 1):
                print(f"  {line_num}: {line}")
    except Exception as e:
        print_colored(Colors.RED, f"Error reading sample sheet: {e}")
    
    try:
        # Build the nextflow command for all barcodes
        cmd = [
            "nextflow", "run", "epi2me-labs/wf-clone-validation",
            "--fastq", input_dir,
            "--sample_sheet", sample_sheet_path,
            "--out_dir", output_dir,
            "--approx_size", "7000",
            "--min_quality", "9",
            "--assm_coverage", "2",
            "--assembly_tool", ASSEMBLY_TOOL,
            "--threads", "4",
            "-resume"
        ]
        
        print_colored(Colors.BLUE, f"Running command: {' '.join(cmd)}")
        print_colored(Colors.YELLOW, "Starting workflow - watch for progress...")
        print()
        
        # Run the workflow with progress spinner
        spinner_message = f"Running {ASSEMBLY_TOOL} assembly on {len(barcodes)} barcodes"
        with ProgressSpinner(spinner_message):
            result = subprocess.run(cmd, capture_output=False, text=True)
        
        # Don't clean up sample sheet immediately for debugging
        print_colored(Colors.BLUE, f"Workflow finished with return code: {result.returncode}")
        
        if result.returncode == 0:
            print_colored(Colors.GREEN, "✓ Quick analysis completed successfully")
            # Clean up sample sheet only on success
            if os.path.exists(sample_sheet_path):
                os.remove(sample_sheet_path)
            return True
        else:
            print_colored(Colors.RED, "✗ Quick analysis failed")
            print_colored(Colors.YELLOW, f"Sample sheet left at: {sample_sheet_path} for debugging")
            return False
            
    except Exception as e:
        print_colored(Colors.RED, f"Error running quick analysis: {str(e)}")
        print_colored(Colors.YELLOW, f"Sample sheet left at: {sample_sheet_path} for debugging")
        return False

def run_workflow(input_dir: str, output_dir: str, barcode: str, 
                approx_size: int, min_quality: int, coverage: int) -> bool:
    """Run the wf-clone-validation workflow for a single barcode"""
    
    print_colored(Colors.BLUE, f"Running workflow for {barcode}...")
    print_colored(Colors.YELLOW, f"Parameters: size={approx_size}, quality={min_quality}, coverage={coverage}, assembler={ASSEMBLY_TOOL}")
    
    # Set up output directory for this barcode
    barcode_output_dir = os.path.join(output_dir, barcode)
    os.makedirs(barcode_output_dir, exist_ok=True)
    
    # Create sample sheet for this barcode
    sample_sheet_path = create_sample_sheet(barcode, approx_size)
    
    try:
        # Build the nextflow command
        cmd = [
            "nextflow", "run", "epi2me-labs/wf-clone-validation",
            "--fastq", input_dir,
            "--sample_sheet", sample_sheet_path,
            "--out_dir", barcode_output_dir,
            "--approx_size", str(approx_size),
            "--min_quality", str(min_quality),
            "--assm_coverage", str(coverage),
            "--assembly_tool", ASSEMBLY_TOOL,
            "--threads", "4",
            "-resume"
        ]
        
        print_colored(Colors.BLUE, f"Running command: {' '.join(cmd)}")
        
        # Run the workflow with progress spinner
        spinner_message = f"Processing {barcode} with {ASSEMBLY_TOOL}"
        with ProgressSpinner(spinner_message):
            result = subprocess.run(cmd, capture_output=False, text=True)
        
        # Clean up sample sheet
        if os.path.exists(sample_sheet_path):
            os.remove(sample_sheet_path)
        
        if result.returncode == 0:
            print_colored(Colors.GREEN, f"✓ Workflow completed successfully for {barcode}")
            return True
        else:
            print_colored(Colors.RED, f"✗ Workflow failed for {barcode}")
            return False
            
    except Exception as e:
        print_colored(Colors.RED, f"Error running workflow for {barcode}: {str(e)}")
        # Clean up sample sheet
        if os.path.exists(sample_sheet_path):
            os.remove(sample_sheet_path)
        return False

def get_barcode_parameters(barcodes: List[str]) -> Dict[str, Tuple[int, int]]:
    """Get parameters for each barcode from user input"""
    parameters = {}
    
    print_colored(Colors.BLUE, "Please provide parameters for each barcode based on the quick run results:")
    print()
    
    for barcode in barcodes:
        print_colored(Colors.YELLOW, f"Parameters for {barcode}:")
        
        size = int(get_user_input("  Plasmid size (bp)", "7000"))
        quality = int(get_user_input("  Minimum quality score", "9"))
        
        parameters[barcode] = (size, quality)
        print()
    
    return parameters

def get_single_barcode_parameters(barcode: str) -> Tuple[int, int, int]:
    """Get parameters for a single barcode (for reprocessing)"""
    print_colored(Colors.YELLOW, f"Parameters for {barcode}:")
    
    size = int(get_user_input("  Plasmid size (bp)", "7000"))
    quality = int(get_user_input("  Minimum quality score", "9"))
    coverage = int(get_user_input("  Coverage depth", "60"))
    
    return size, quality, coverage

def check_dependencies() -> bool:
    """Check if required dependencies are installed"""
    dependencies = ["nextflow"]
    missing = []
    
    for dep in dependencies:
        if shutil.which(dep) is None:
            missing.append(dep)
    
    if missing:
        print_colored(Colors.RED, f"Missing dependencies: {', '.join(missing)}")
        print_colored(Colors.YELLOW, "Please install the missing dependencies and try again.")
        return False
    
    return True

def main():
    """Main script execution"""
    print_colored(Colors.GREEN, "=== Plasmid Clone Validation Workflow Runner ===")
    print()
    
    # Check dependencies
    if not check_dependencies():
        sys.exit(1)
    
    # Get input directory
    input_path = get_user_input("Enter the path to the parent folder containing barcode subfolders")
    input_dir = convert_windows_path(input_path)
    
    if not validate_directory(input_dir, "Input"):
        sys.exit(1)
    
    # Get output directory
    output_path = get_user_input("Enter the path to the results folder")
    output_dir = convert_windows_path(output_path)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Get barcodes
    print_colored(Colors.BLUE, "Scanning for barcodes...")
    barcodes = get_barcodes(input_dir)
    
    if not barcodes:
        print_colored(Colors.RED, f"No barcodes with FASTQ files found in {input_dir}")
        sys.exit(1)
    
    print_colored(Colors.GREEN, f"Found barcodes: {', '.join(barcodes)}")
    print()
    
    # Phase 1: Quick run on ALL barcodes (ALWAYS use Flye for speed)
    print_colored(Colors.BLUE, "=== PHASE 1: Quick Analysis ===")
    print_colored(Colors.YELLOW, "Running quick analysis on ALL barcodes with Flye (fast) to get size estimates...")
    print_colored(Colors.YELLOW, "Parameters: size=7000, quality=9, coverage=2, assembler=flye")
    
    # Set Flye for quick analysis
    global ASSEMBLY_TOOL
    ASSEMBLY_TOOL = "flye"
    
    # Create quick results directory
    quick_results_dir = os.path.join(output_dir, "quick_results")
    os.makedirs(quick_results_dir, exist_ok=True)
    
    # Run quick analysis on ALL barcodes with Flye
    if run_quick_analysis_all_barcodes(input_dir, quick_results_dir, barcodes):
        print_colored(Colors.GREEN, "Quick analysis completed for all barcodes!")
        print_colored(Colors.YELLOW, f"Please examine the results in: {quick_results_dir}")
        print_colored(Colors.YELLOW, "Look at the wf-clone-validation-report.html and read length distribution plots for each barcode")
        print_colored(Colors.YELLOW, "Note the optimal plasmid size and quality scores for each barcode")
        print()
        input("Press Enter when you've examined the results and are ready to proceed...")
    else:
        print_colored(Colors.RED, "Quick analysis failed. Please check the error messages above.")
        sys.exit(1)
    
    # Phase 2: Get parameters and run full analysis
    print_colored(Colors.BLUE, "=== PHASE 2: Full Analysis ===")
    print()
    
    # NOW ask for assembly tool preference for the full analysis
    print_colored(Colors.BLUE, "Choose assembly tool for the full analysis:")
    print("1. Flye (default, faster)")
    print("2. Canu (slower but more consistent)")
    assembly_choice = get_user_input("Select assembly tool (1/2)", "1")
    
    if assembly_choice == "2":
        ASSEMBLY_TOOL = "canu"
        print_colored(Colors.GREEN, "Using Canu assembler for full analysis")
        print_colored(Colors.YELLOW, "Note: Canu will be much slower but more consistent")
    else:
        ASSEMBLY_TOOL = "flye"
        print_colored(Colors.GREEN, "Using Flye assembler for full analysis")
    
    print()
    
    # Get parameters for each barcode
    barcode_params = get_barcode_parameters(barcodes)
    
    # Get global coverage parameter
    global_coverage = int(get_user_input("Enter the coverage depth to apply to all barcodes", "60"))
    
    # Confirmation before proceeding
    print()
    print_colored(Colors.YELLOW, f"Ready to start full analysis with {ASSEMBLY_TOOL}. This will:")
    print_colored(Colors.YELLOW, "- Erase temporary quick results")
    print_colored(Colors.YELLOW, "- Process each barcode with specified parameters")
    print_colored(Colors.YELLOW, "- Save results in separate subfolders")
    print()
    input("Press Enter to continue (Ctrl+C to abort)...")
    
    # Clean up quick results
    print_colored(Colors.BLUE, "Cleaning up temporary results...")
    shutil.rmtree(quick_results_dir, ignore_errors=True)
    
    # Process each barcode individually
    for barcode in barcodes:
        size, quality = barcode_params[barcode]
        
        if run_workflow(input_dir, output_dir, barcode, size, quality, global_coverage):
            print_colored(Colors.GREEN, f"✓ Completed processing {barcode}")
        else:
            print_colored(Colors.RED, f"✗ Failed to process {barcode}")
        print()
    
    print_colored(Colors.GREEN, "=== All barcodes processed! ===")
    
    # Reprocessing loop
    while True:
        print()
        reprocess = get_user_input("Do you need to reprocess any barcode? (y/N)", "n").lower()
        
        if reprocess not in ['y', 'yes']:
            break
        
        # Show available barcodes
        print_colored(Colors.BLUE, f"Available barcodes: {', '.join(barcodes)}")
        target_barcode = get_user_input("Which barcode do you want to reprocess?")
        
        # Validate barcode exists
        if target_barcode not in barcodes:
            print_colored(Colors.RED, f"Invalid barcode: {target_barcode}")
            continue
        
        # Get parameters for reprocessing
        size, quality, coverage = get_single_barcode_parameters(target_barcode)
        
        # Remove existing results for this barcode
        target_output = os.path.join(output_dir, target_barcode)
        if os.path.exists(target_output):
            print_colored(Colors.YELLOW, f"Removing existing results for {target_barcode}...")
            shutil.rmtree(target_output, ignore_errors=True)
        
        # Reprocess
        if run_workflow(input_dir, output_dir, target_barcode, size, quality, coverage):
            print_colored(Colors.GREEN, f"✓ Successfully reprocessed {target_barcode}")
        else:
            print_colored(Colors.RED, f"✗ Failed to reprocess {target_barcode}")
    
    print_colored(Colors.GREEN, "=== Workflow complete! ===")
    print_colored(Colors.BLUE, f"Results are saved in: {output_dir}")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print_colored(Colors.YELLOW, "\nScript interrupted by user. Exiting...")
        sys.exit(0)
    except Exception as e:
        print_colored(Colors.RED, f"An unexpected error occurred: {str(e)}")
        sys.exit(1)
