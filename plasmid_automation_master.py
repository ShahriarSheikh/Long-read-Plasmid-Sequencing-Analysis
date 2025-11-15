#!/usr/bin/env python3
"""
Plasmid Analysis Automation Pipeline - Standalone Version
Performs Steps 2-5 of the analysis workflow (post plasmid_runner)
Run plasmid_runner.py separately first, then use this script
"""

import os
import sys
import subprocess
import shutil
import glob
import re
import tempfile
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO

# Color codes for terminal output
class Colors:
    RED = '\033[0;31m'
    GREEN = '\033[0;32m'
    YELLOW = '\033[1;33m'
    BLUE = '\033[0;34m'
    CYAN = '\033[0;36m'
    BOLD = '\033[1m'
    NC = '\033[0m'

def print_colored(color: str, message: str) -> None:
    """Print colored message to terminal"""
    print(f"{color}{message}{Colors.NC}")

def print_header(title: str) -> None:
    """Print section header"""
    print("\n" + "=" * 70)
    print_colored(Colors.CYAN + Colors.BOLD, f"  {title}")
    print("=" * 70)

def handle_error(error_message: str, barcode: str = None) -> bool:
    """
    Handle an error by asking user whether to continue or stop
    Returns: True to continue, False to stop
    """
    print_colored(Colors.RED, f"✗ ERROR: {error_message}")
    
    if barcode:
        print_colored(Colors.YELLOW, f"\nProblem with barcode: {barcode}")
    
    print()
    print_colored(Colors.YELLOW, "What would you like to do?")
    print("  1. Skip this barcode and continue with others")
    print("  2. Stop the automation")
    
    while True:
        choice = input("Enter your choice (1/2): ").strip()
        
        if choice == "1":
            print_colored(Colors.BLUE, f"⏭️  Skipping {barcode if barcode else 'this step'} and continuing...")
            return True
        elif choice == "2":
            print_colored(Colors.RED, "🛑 Stopping automation as requested.")
            return False
        else:
            print_colored(Colors.RED, "Invalid choice. Please enter 1 or 2.")

def convert_windows_path(path: str) -> str:
    """Convert Windows path to WSL/Linux path format"""
    path = path.strip().strip('"').strip("'")
    
    # Handle WSL localhost paths - multiple variations
    # \\wsl.localhost\Ubuntu\path or //wsl.localhost/Ubuntu/path
    if '\\wsl.localhost\\Ubuntu\\' in path or '//wsl.localhost/Ubuntu/' in path or '\\\\wsl.localhost\\Ubuntu\\' in path:
        # Remove all WSL localhost prefixes
        path = path.replace('\\\\wsl.localhost\\Ubuntu\\', '/')
        path = path.replace('\\wsl.localhost\\Ubuntu\\', '/')
        path = path.replace('//wsl.localhost/Ubuntu/', '/')
        path = path.replace('\\wsl.localhost\\Ubuntu', '')
        path = path.replace('//wsl.localhost/Ubuntu', '')
        # Convert remaining backslashes to forward slashes
        path = path.replace('\\', '/')
        # Ensure it starts with /
        if not path.startswith('/'):
            path = '/' + path
        return path
    
    # Check if it's a Windows path (contains drive letter)
    if len(path) >= 2 and path[1] == ':' and path[0].isalpha():
        drive = path[0].lower()
        rest = path[2:].replace('\\', '/')
        if rest.startswith('/'):
            rest = rest[1:]
        return f"/mnt/{drive}/{rest}"
    
    # Default: just convert backslashes to forward slashes
    return path.replace('\\', '/')

def get_path_from_user(prompt: str, description: str, must_exist: bool = True) -> str:
    """Get and validate a path from user"""
    print()
    print_colored(Colors.CYAN, f"=== {description} ===")
    
    while True:
        path = input(f"{prompt}: ").strip()
        if not path:
            print_colored(Colors.RED, "Path cannot be empty. Please try again.")
            continue
        
        linux_path = convert_windows_path(path)
        
        if must_exist and not os.path.exists(linux_path):
            print_colored(Colors.RED, f"✗ Path not found: {linux_path}")
            retry = input("Try again? (y/n): ").strip().lower()
            if retry != 'y':
                sys.exit(1)
            continue
        
        print_colored(Colors.GREEN, f"✓ Path accepted: {linux_path}")
        return linux_path

def get_barcodes_from_output(output_dir: str) -> List[str]:
    """Get list of barcode folders from output directory"""
    barcodes = []
    for item in os.listdir(output_dir):
        item_path = os.path.join(output_dir, item)
        if os.path.isdir(item_path) and item not in ['quick_results', 'work', '.nextflow']:
            barcodes.append(item)
    return sorted(barcodes)

def get_search_date() -> Optional[datetime]:
    """Get the date to use for searching downsampled FASTQ files"""
    print()
    print_colored(Colors.CYAN, "=== Downsampled FASTQ Search Date ===")
    print_colored(Colors.YELLOW, "Enter the date when plasmid_runner was executed.")
    print_colored(Colors.YELLOW, "This helps filter downsampled FASTQ files by modification date.")
    print()
    
    date_formats = [
        ('%Y-%m-%d', '2025-10-17'),
        ('%m/%d/%Y', '10/17/2025'),
        ('%d/%m/%Y', '17/10/2025'),
        ('%Y%m%d', '20251017'),
    ]
    
    print("Accepted formats:")
    for fmt, example in date_formats:
        print(f"  • {example}")
    print()
    
    while True:
        date_input = input("Enter search date (or press Enter to skip date filtering): ").strip()
        
        if date_input == "":
            print_colored(Colors.YELLOW, "⚠ No date specified - will search all FASTQ files without date filtering")
            return None
        
        # Try parsing with different formats
        for fmt, _ in date_formats:
            try:
                parsed_date = datetime.strptime(date_input, fmt)
                print_colored(Colors.GREEN, f"✓ Date accepted: {parsed_date.strftime('%Y-%m-%d')}")
                return parsed_date
            except ValueError:
                continue
        
        print_colored(Colors.RED, "✗ Invalid date format. Please try again.")

def get_plasmid_sizes(barcodes: List[str]) -> Dict[str, int]:
    """Get plasmid sizes for each barcode from user"""
    print()
    print_colored(Colors.CYAN, "=== Plasmid Sizes ===")
    print_colored(Colors.YELLOW, "Please provide the plasmid size for each barcode:")
    print()
    
    plasmid_sizes = {}
    for barcode in barcodes:
        while True:
            try:
                size_input = input(f"  {barcode} plasmid size (bp) [default: 7000]: ").strip()
                if size_input == "":
                    size = 7000
                else:
                    size = int(size_input)
                
                if size <= 0:
                    print_colored(Colors.RED, "Plasmid size must be positive")
                    continue
                
                plasmid_sizes[barcode] = size
                break
            except ValueError:
                print_colored(Colors.RED, "Please enter a valid number")
    
    print()
    print_colored(Colors.GREEN, "✓ Plasmid sizes recorded:")
    for barcode, size in plasmid_sizes.items():
        print(f"    {barcode}: {size} bp")
    
    return plasmid_sizes

def find_file_in_folder(folder: str, pattern: str) -> Optional[str]:
    """Find a file matching pattern in folder"""
    matches = glob.glob(os.path.join(folder, pattern))
    if matches:
        return matches[0]
    return None

def extract_barcode_and_date(fasta_path: str) -> Tuple[str, datetime]:
    """Extract barcode name and creation date from FASTA file"""
    filename = os.path.basename(fasta_path)
    barcode_name = filename.replace('.final.fasta', '')
    
    file_stat = os.stat(fasta_path)
    file_date = datetime.fromtimestamp(file_stat.st_mtime).date()
    
    return barcode_name, file_date

def find_downsampled_fastq(barcode_name: str, search_date: Optional[datetime], work_dir: str) -> Optional[str]:
    """Find the correct downsampled FASTQ file"""
    patterns = [
        f"**/{barcode_name}.downsampled*.fastq",
        f"**/{barcode_name}.downsampled*.fq",
        f"**/*{barcode_name}*downsampled*.fastq",
        f"**/*{barcode_name}*downsampled*.fq",
    ]
    
    found_files = []
    for pattern in patterns:
        full_pattern = os.path.join(work_dir, pattern)
        matches = glob.glob(full_pattern, recursive=True)
        for match in matches:
            if os.path.isfile(match):
                found_files.append(match)
    
    found_files = list(set(found_files))
    
    if not found_files:
        return None
    
    # Filter by date if provided
    if search_date:
        date_filtered = []
        target_date = search_date.date()
        
        for file_path in found_files:
            file_stat = os.stat(file_path)
            file_mod_date = datetime.fromtimestamp(file_stat.st_mtime).date()
            if file_mod_date == target_date:
                date_filtered.append(file_path)
        
        if date_filtered:
            found_files = date_filtered
        else:
            # If no exact date match, show warning but continue with all files
            print_colored(Colors.YELLOW, f"    ⚠ No files found matching date {target_date}, using all found files")
    
    # Sort by size (largest first), then by time (most recent first)
    found_files.sort(key=lambda x: (os.path.getsize(x), os.path.getmtime(x)), reverse=True)
    
    return found_files[0]

def run_step2_pysam_analysis(output_dir: str, work_dir: str, barcodes: List[str], search_date: Optional[datetime]) -> None:
    """Run detailed pysam analysis for each barcode"""
    print_header("STEP 2: Detailed Pysam-Based Consensus Analysis")
    
    if search_date:
        print_colored(Colors.BLUE, f"Using search date: {search_date.strftime('%Y-%m-%d')}")
    else:
        print_colored(Colors.YELLOW, "No date filter applied - searching all FASTQ files")
    
    for idx, barcode in enumerate(barcodes, 1):
        print()
        print_colored(Colors.CYAN, f"[{idx}/{len(barcodes)}] Processing barcode: {barcode}")
        
        barcode_folder = os.path.join(output_dir, barcode)
        
        fasta_file = find_file_in_folder(barcode_folder, "*.final.fasta")
        if not fasta_file:
            if not handle_error(f"No *.final.fasta file found in {barcode_folder}", barcode):
                sys.exit(1)
            continue
        
        print_colored(Colors.GREEN, f"✓ Found FASTA: {os.path.basename(fasta_file)}")
        
        barcode_name, file_date = extract_barcode_and_date(fasta_file)
        print_colored(Colors.BLUE, f"  Barcode name: {barcode_name}")
        print_colored(Colors.BLUE, f"  FASTA file date: {file_date}")
        
        print_colored(Colors.BLUE, "  Searching for downsampled FASTQ file...")
        fastq_file = find_downsampled_fastq(barcode_name, search_date, work_dir)
        
        if not fastq_file:
            if not handle_error(f"No matching downsampled FASTQ file found for {barcode}", barcode):
                sys.exit(1)
            continue
        
        print_colored(Colors.GREEN, f"✓ Found FASTQ: {os.path.basename(fastq_file)}")
        print_colored(Colors.BLUE, f"  Size: {os.path.getsize(fastq_file) / (1024*1024):.1f} MB")
        
        dest_fastq = os.path.join(barcode_folder, os.path.basename(fastq_file))
        print_colored(Colors.BLUE, "  Copying FASTQ to barcode folder...")
        shutil.copy2(fastq_file, dest_fastq)
        print_colored(Colors.GREEN, f"✓ FASTQ copied")
        
        print_colored(Colors.BLUE, "  Running minimap2 alignment...")
        bam_file = run_minimap2_alignment(fasta_file, dest_fastq, barcode_folder)
        
        if not bam_file:
            if not handle_error(f"Alignment failed for {barcode}", barcode):
                sys.exit(1)
            continue
        
        print_colored(Colors.GREEN, f"✓ Alignment completed")
        
        print_colored(Colors.BLUE, "  Running pysam analysis...")
        csv_file = analyze_with_pysam(bam_file, fasta_file, barcode_folder)
        
        if not csv_file:
            if not handle_error(f"Pysam analysis failed for {barcode}", barcode):
                sys.exit(1)
            continue
        
        print_colored(Colors.GREEN, f"✓ Analysis completed: Consensus Sequence Confidence Report.csv")

def run_minimap2_alignment(fasta_file: str, fastq_file: str, output_dir: str) -> Optional[str]:
    """Run minimap2 alignment"""
    output_sam = os.path.join(output_dir, "alignment.sam")
    output_bam = os.path.join(output_dir, "alignment.bam")
    output_sorted = os.path.join(output_dir, "alignment_sorted.bam")
    
    try:
        minimap2_cmd = ['minimap2', '-ax', 'map-ont', '-y', fasta_file, fastq_file]
        with open(output_sam, 'w') as sam_out:
            subprocess.run(minimap2_cmd, stdout=sam_out, stderr=subprocess.PIPE, check=True)
        
        subprocess.run(['samtools', 'view', '-bS', output_sam, '-o', output_bam], check=True)
        subprocess.run(['samtools', 'sort', output_bam, '-o', output_sorted], check=True)
        subprocess.run(['samtools', 'index', output_sorted], check=True)
        
        os.remove(output_sam)
        os.remove(output_bam)
        
        return output_sorted
    except Exception as e:
        print_colored(Colors.RED, f"Alignment error: {e}")
        return None

def analyze_with_pysam(bam_file: str, fasta_file: str, output_dir: str) -> Optional[str]:
    """Perform pysam analysis"""
    try:
        import pysam
        
        ref_seq = ""
        with open(fasta_file, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    ref_seq += line.strip().upper()
        
        results = []
        
        with pysam.AlignmentFile(bam_file, "rb") as samfile:
            for pileupcolumn in samfile.pileup(stepper='all'):
                if pileupcolumn.pos >= len(ref_seq):
                    continue
                
                pos = pileupcolumn.pos + 1
                ref_base = ref_seq[pileupcolumn.pos]
                
                pos_stats = {
                    'pos': pos,
                    'ref': ref_base,
                    'reads_all': pileupcolumn.nsegments,
                    'matches': 0,
                    'mismatches': 0,
                    'deletions': 0,
                    'insertions': 0,
                    'A': 0, 'C': 0, 'T': 0, 'G': 0,
                    'homopolymer': 0,
                    'methylation': 0
                }
                
                for pileupread in pileupcolumn.pileups:
                    if pileupread.is_del:
                        pos_stats['deletions'] += 1
                    elif not pileupread.is_refskip and pileupread.query_position is not None:
                        read_base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                        if read_base in 'ATGC':
                            pos_stats[read_base] += 1
                            if read_base == ref_base:
                                pos_stats['matches'] += 1
                            else:
                                pos_stats['mismatches'] += 1
                    
                    if pileupread.indel > 0:
                        pos_stats['insertions'] += 1
                
                pos_stats['homopolymer'] = calculate_homopolymer_length(ref_seq, pileupcolumn.pos)
                
                results.append(pos_stats)
        
        df = pd.DataFrame(results)
        csv_file = os.path.join(output_dir, "Consensus Sequence Confidence Report.csv")
        df.to_csv(csv_file, index=False)
        
        return csv_file
    except Exception as e:
        print_colored(Colors.RED, f"Pysam analysis error: {e}")
        return None

def calculate_homopolymer_length(sequence: str, position: int) -> int:
    """Calculate homopolymer length at position"""
    if position >= len(sequence):
        return 0
    
    base = sequence[position]
    if base not in 'ATGC':
        return 0
    
    left_count = 0
    for i in range(position - 1, -1, -1):
        if sequence[i] == base:
            left_count += 1
        else:
            break
    
    right_count = 0
    for i in range(position + 1, len(sequence)):
        if sequence[i] == base:
            right_count += 1
        else:
            break
    
    total_length = left_count + 1 + right_count
    return total_length if total_length >= 2 else 0

def run_step3_fasta_to_ab1(output_dir: str, barcodes: List[str]) -> None:
    """Convert FASTA files to AB1 format"""
    print_header("STEP 3: FASTA to AB1 Conversion")
    
    for idx, barcode in enumerate(barcodes, 1):
        print()
        print_colored(Colors.CYAN, f"[{idx}/{len(barcodes)}] Processing barcode: {barcode}")
        
        barcode_folder = os.path.join(output_dir, barcode)
        
        fasta_file = find_file_in_folder(barcode_folder, "*.final.fasta")
        if not fasta_file:
            if not handle_error(f"No *.final.fasta file found in {barcode_folder}", barcode):
                sys.exit(1)
            continue
        
        print_colored(Colors.BLUE, f"  Converting {os.path.basename(fasta_file)} to AB1...")
        
        seq_count = 0
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_count += 1
        
        if seq_count > 1:
            if not handle_error(f"Multiple sequences found in FASTA file", barcode):
                sys.exit(1)
            continue
        
        filename = os.path.basename(fasta_file).replace('.final.fasta', '')
        output_ab1 = os.path.join(barcode_folder, f"{filename}.ab1")
        
        try:
            for record in SeqIO.parse(fasta_file, "fasta"):
                create_ab1_file(record.seq, output_ab1, record.id)
                if len(record.seq) > 800:
                    print_colored(Colors.YELLOW, f"  ⚠ Sequence truncated to 800bp for AB1 compatibility")
        except Exception as e:
            if not handle_error(f"AB1 conversion failed: {e}", barcode):
                sys.exit(1)
            continue
        
        print_colored(Colors.GREEN, f"✓ AB1 file created: {os.path.basename(output_ab1)}")

def create_ab1_file(sequence, output_path: str, sample_name: str) -> None:
    """Create AB1 file with realistic chromatogram"""
    import struct
    
    seq_str = str(sequence).upper()[:800]
    seq_length = len(seq_str)
    
    scan_rate = 24
    total_scans = seq_length * scan_rate
    
    trace_a = np.zeros(total_scans, dtype=np.int16)
    trace_c = np.zeros(total_scans, dtype=np.int16)
    trace_g = np.zeros(total_scans, dtype=np.int16)
    trace_t = np.zeros(total_scans, dtype=np.int16)
    
    np.random.seed(hash(seq_str) % 2**32)
    
    for i, base in enumerate(seq_str):
        center_scan = i * scan_rate + scan_rate // 2
        if center_scan < total_scans:
            peak_height = np.random.randint(15000, 25000)
            peak_width = 8
            
            for offset in range(-peak_width, peak_width + 1):
                scan_pos = center_scan + offset
                if 0 <= scan_pos < total_scans:
                    intensity = int(peak_height * np.exp(-0.5 * (offset / 3.0) ** 2))
                    
                    if base == 'A':
                        trace_a[scan_pos] = min(trace_a[scan_pos] + intensity, 32000)
                    elif base == 'C':
                        trace_c[scan_pos] = min(trace_c[scan_pos] + intensity, 32000)
                    elif base == 'G':
                        trace_g[scan_pos] = min(trace_g[scan_pos] + intensity, 32000)
                    elif base == 'T':
                        trace_t[scan_pos] = min(trace_t[scan_pos] + intensity, 32000)
    
    baseline = 200
    noise_level = 50
    for i in range(total_scans):
        noise = np.random.randint(0, noise_level + 1)
        trace_a[i] = min(trace_a[i] + baseline + noise, 32000)
        trace_c[i] = min(trace_c[i] + baseline + noise, 32000)
        trace_g[i] = min(trace_g[i] + baseline + noise, 32000)
        trace_t[i] = min(trace_t[i] + baseline + noise, 32000)
    
    peak_locations = [i * scan_rate + scan_rate // 2 for i in range(seq_length)]
    
    entries = []
    data_offset = 128
    
    def add_entry(name, number, data_type, element_size, data):
        nonlocal data_offset, entries
        while data_offset % 4 != 0:
            data_offset += 1
        
        data_bytes = b''
        num_elements = 0
        
        if isinstance(data, str):
            data_bytes = data.encode('ascii') + b'\x00'
            num_elements = len(data_bytes)
        elif isinstance(data, (list, np.ndarray)):
            data_bytes = np.array(data, dtype='>i2').tobytes()
            num_elements = len(data)
        
        entry = {
            'name': name[:4].ljust(4, '\x00'), 'number': number,
            'data_type': data_type, 'element_size': element_size,
            'num_elements': num_elements, 'data_size': len(data_bytes),
            'data_offset': data_offset, 'data': data_bytes
        }
        entries.append(entry)
        data_offset += len(data_bytes)
        return entry
    
    add_entry('PBAS', 1, 2, 1, seq_str)
    add_entry('PLOC', 1, 4, 2, peak_locations)
    add_entry('DATA', 9, 4, 2, trace_a)
    add_entry('DATA', 10, 4, 2, trace_c)
    add_entry('DATA', 11, 4, 2, trace_g)
    add_entry('DATA', 12, 4, 2, trace_t)
    add_entry('SMPL', 1, 19, 1, sample_name)
    add_entry('MCHN', 1, 19, 1, 'Python Automation')
    add_entry('FWO_', 1, 2, 1, 'ACGT')
    
    with open(output_path, 'wb') as f:
        directory_pos = data_offset
        while directory_pos % 4 != 0:
            directory_pos += 1
        
        f.write(b'ABIF')
        f.write(struct.pack('>H', 101))
        f.write(b'TDIR')
        f.write(struct.pack('>I', 1))
        f.write(struct.pack('>H', 1027))
        f.write(struct.pack('>H', 28))
        f.write(struct.pack('>I', len(entries)))
        f.write(struct.pack('>I', len(entries) * 28))
        f.write(struct.pack('>I', directory_pos))
        f.write(struct.pack('>I', 0))
        f.write(b'\x00' * (128 - f.tell()))
        
        for entry in entries:
            while f.tell() < entry['data_offset']:
                f.write(b'\x00')
            entry['actual_offset'] = f.tell()
            f.write(entry['data'])
        
        while f.tell() < directory_pos:
            f.write(b'\x00')
        
        for entry in entries:
            name_bytes = entry['name'].encode('ascii')
            f.write(name_bytes)
            f.write(struct.pack('>I', entry['number']))
            f.write(struct.pack('>H', entry['data_type']))
            f.write(struct.pack('>H', entry['element_size']))
            f.write(struct.pack('>I', entry['num_elements']))
            f.write(struct.pack('>I', entry['data_size']))
            f.write(struct.pack('>I', entry['actual_offset']))
            f.write(struct.pack('>I', 0))

def run_step4_coverage_plots(output_dir: str, barcodes: List[str]) -> None:
    """Generate coverage plots for each barcode"""
    print_header("STEP 4: Coverage Plot Generation")
    
    for idx, barcode in enumerate(barcodes, 1):
        print()
        print_colored(Colors.CYAN, f"[{idx}/{len(barcodes)}] Processing barcode: {barcode}")
        
        barcode_folder = os.path.join(output_dir, barcode)
        
        csv_file = os.path.join(barcode_folder, "Consensus Sequence Confidence Report.csv")
        if not os.path.exists(csv_file):
            if not handle_error(f"CSV file not found: {csv_file}", barcode):
                sys.exit(1)
            continue
        
        print_colored(Colors.BLUE, f"  Generating coverage plot from CSV...")
        
        try:
            df = pd.read_csv(csv_file)
            
            plt.figure(figsize=(12, 6))
            plt.fill_between(df['pos'], 0, df['reads_all'], color='black', alpha=1.0)
            plt.xlabel('Position', fontsize=12)
            plt.ylabel('Coverage (reads_all)', fontsize=12)
            plt.ylim(bottom=0)
            plt.xlim(left=0)
            plt.tight_layout()
            
            output_path = os.path.join(barcode_folder, "Read Coverage Plot.png")
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            print_colored(Colors.GREEN, f"✓ Coverage plot saved: Read Coverage Plot.png")
            
        except Exception as e:
            if not handle_error(f"Plot generation failed: {e}", barcode):
                sys.exit(1)
            continue

def run_step5_clone_validation(input_dir: str, output_dir: str, barcodes: List[str], 
                                plasmid_sizes: Dict[str, int], ecoli_ref: str) -> None:
    """Run clone validation for each barcode"""
    print_header("STEP 5: Clone Validation Reports")
    
    for idx, barcode in enumerate(barcodes, 1):
        print()
        print_colored(Colors.CYAN, f"[{idx}/{len(barcodes)}] Processing barcode: {barcode}")
        
        barcode_input = os.path.join(input_dir, barcode)
        barcode_output = os.path.join(output_dir, barcode)
        
        print_colored(Colors.BLUE, f"  Reading FASTQ files from {barcode}...")
        
        fastq_files = glob.glob(os.path.join(barcode_input, "*.fastq*"))
        
        if not fastq_files:
            if not handle_error(f"No FASTQ files found in {barcode_input}", barcode):
                sys.exit(1)
            continue
        
        all_sequences = []
        for fastq_file in fastq_files:
            if fastq_file.endswith('.gz'):
                import gzip
                handle = gzip.open(fastq_file, 'rt')
            else:
                handle = open(fastq_file, 'r')
            
            for record in SeqIO.parse(handle, "fastq"):
                all_sequences.append(str(record.seq))
            handle.close()
        
        print_colored(Colors.BLUE, f"  Analyzing {len(all_sequences)} sequences...")
        
        plasmid_size = plasmid_sizes.get(barcode, 7000)
        
        mole_pct, mass_pct = calculate_nmer_distribution(all_sequences, plasmid_size)
        
        contamination = check_contamination_minimap2(all_sequences, ecoli_ref)
        
        report_path = os.path.join(barcode_output, "CLONE VALIDATION REPORT.txt")
        
        with open(report_path, 'w') as f:
            f.write("CLONE VALIDATION REPORT\n")
            f.write(f"Barcode: {barcode}\n")
            f.write("=" * 50 + "\n")
            f.write(f"Total sequences: {len(all_sequences):,}\n")
            f.write(f"Total bases: {sum(len(s) for s in all_sequences):,}\n")
            f.write(f"Plasmid size used: {plasmid_size} bp\n\n")
            f.write("N-mer Distribution:\n")
            f.write(f"{'':>10} {'1-mer (%)':>10} {'2-mer (%)':>10} {'3-mer (%)':>10}\n")
            f.write(f"{'moles':>10} {mole_pct[1]:>9.1f} {mole_pct[2]:>9.1f} {mole_pct[3]:>9.1f}\n")
            f.write(f"{'mass':>10} {mass_pct[1]:>9.1f} {mass_pct[2]:>9.1f} {mass_pct[3]:>9.1f}\n\n")
            f.write(f"***** E. coli genomic contamination: {contamination:.1f}%\n")
            f.write("=" * 50 + "\n")
        
        print_colored(Colors.GREEN, f"✓ Clone validation report saved: CLONE VALIDATION REPORT.txt")

def calculate_nmer_distribution(sequences: List[str], plasmid_size: int) -> Tuple[Dict, Dict]:
    """Calculate n-mer distribution"""
    total_sequences = len(sequences)
    total_bases = sum(len(seq) for seq in sequences)
    
    if total_sequences == 0 or plasmid_size == 0:
        return {1:0, 2:0, 3:0}, {1:0, 2:0, 3:0}
    
    monomer_min = plasmid_size * 0.8
    monomer_max = plasmid_size * 1.2
    dimer_min = plasmid_size * 1.8
    dimer_max = plasmid_size * 2.2
    multimer_min = plasmid_size * 2.5
    
    monomers = []
    dimers = []
    multimers = []
    
    for seq in sequences:
        seq_len = len(seq)
        if monomer_min <= seq_len <= monomer_max:
            monomers.append(seq)
        elif dimer_min <= seq_len <= dimer_max:
            dimers.append(seq)
        elif seq_len >= multimer_min:
            multimers.append(seq)
        else:
            if seq_len < monomer_min:
                monomers.append(seq)
            elif seq_len < dimer_min:
                if abs(seq_len - plasmid_size) < abs(seq_len - (2 * plasmid_size)):
                    monomers.append(seq)
                else:
                    dimers.append(seq)
            else:
                multimers.append(seq)
    
    mole_percentages = {
        1: (len(monomers) / total_sequences) * 100,
        2: (len(dimers) / total_sequences) * 100,
        3: (len(multimers) / total_sequences) * 100
    }
    
    monomer_bases = sum(len(seq) for seq in monomers)
    dimer_bases = sum(len(seq) for seq in dimers)
    multimer_bases = sum(len(seq) for seq in multimers)
    
    mass_percentages = {
        1: (monomer_bases / total_bases) * 100,
        2: (dimer_bases / total_bases) * 100,
        3: (multimer_bases / total_bases) * 100
    }
    
    return mole_percentages, mass_percentages

def check_contamination_minimap2(sequences: List[str], ecoli_ref: str) -> float:
    """Check contamination using minimap2"""
    try:
        sample_size = min(5000, len(sequences))
        sample_sequences = sequences[::len(sequences)//sample_size] if len(sequences) > sample_size else sequences
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_fasta:
            temp_fasta_path = temp_fasta.name
            for i, seq in enumerate(sample_sequences):
                temp_fasta.write(f">seq_{i}\n{seq}\n")
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.paf', delete=False) as temp_paf:
            temp_paf_path = temp_paf.name
        
        cmd = ['minimap2', '-x', 'map-ont', '--secondary=no', '-t', '4', ecoli_ref, temp_fasta_path]
        with open(temp_paf_path, 'w') as paf_file:
            subprocess.run(cmd, stdout=paf_file, stderr=subprocess.PIPE)
        
        contaminated_seqs = set()
        with open(temp_paf_path, 'r') as paf:
            for line in paf:
                if line.strip():
                    fields = line.strip().split('\t')
                    seq_id = fields[0]
                    query_len = int(fields[1])
                    alignment_len = int(fields[10])
                    if alignment_len > query_len * 0.3:
                        contaminated_seqs.add(seq_id)
        
        contamination_percent = (len(contaminated_seqs) / len(sample_sequences)) * 100
        
        os.unlink(temp_fasta_path)
        os.unlink(temp_paf_path)
        
        return contamination_percent
    except Exception:
        return 0.0

def main():
    """Main automation pipeline"""
    print_colored(Colors.GREEN + Colors.BOLD, "\n" + "="*70)
    print_colored(Colors.GREEN + Colors.BOLD, "  PLASMID ANALYSIS AUTOMATION - STANDALONE VERSION")
    print_colored(Colors.GREEN + Colors.BOLD, "  (Run plasmid_runner.py separately first)")
    print_colored(Colors.GREEN + Colors.BOLD, "="*70)
    
    print()
    print_colored(Colors.YELLOW, "This script performs Steps 2-5 of the analysis:")
    print("  • Step 2: Pysam-based consensus analysis")
    print("  • Step 3: FASTA to AB1 conversion")
    print("  • Step 4: Coverage plot generation")
    print("  • Step 5: Clone validation reports")
    print()
    print_colored(Colors.CYAN, "You will be asked to provide all necessary paths at the beginning.")
    print()
    
    # Get all required paths
    input_dir = get_path_from_user(
        "Enter INPUT folder (barcode parent folder with original FASTQ files)",
        "Input Directory",
        must_exist=True
    )
    
    output_dir = get_path_from_user(
        "Enter OUTPUT folder (where plasmid_runner saved results)",
        "Output Directory",
        must_exist=True
    )
    
    work_dir = get_path_from_user(
        "Enter WORK directory (for finding downsampled FASTQ files)",
        "Work Directory",
        must_exist=True
    )
    
    ecoli_ref = get_path_from_user(
        "Enter E. coli reference genome path",
        "E. coli Reference",
        must_exist=True
    )
    
    # Detect barcodes
    print()
    print_colored(Colors.CYAN, "=== Detecting Barcodes ===")
    barcodes = get_barcodes_from_output(output_dir)
    
    if not barcodes:
        print_colored(Colors.RED, f"✗ No barcode folders found in {output_dir}")
        sys.exit(1)
    
    print_colored(Colors.GREEN, f"✓ Found {len(barcodes)} barcodes:")
    for bc in barcodes:
        print(f"    • {bc}")
    
    # Get plasmid sizes
    plasmid_sizes = get_plasmid_sizes(barcodes)
    
    # Get search date for FASTQ files
    search_date = get_search_date()
    
    # Confirmation
    print()
    print_header("READY TO START")
    print_colored(Colors.YELLOW, "Summary:")
    print(f"  Input directory:  {input_dir}")
    print(f"  Output directory: {output_dir}")
    print(f"  Work directory:   {work_dir}")
    print(f"  E. coli ref:      {ecoli_ref}")
    print(f"  Barcodes:         {len(barcodes)}")
    print(f"  Search date:      {search_date.strftime('%Y-%m-%d') if search_date else 'No date filter'}")
    print()
    
    response = input("Start automation? (y/n): ").strip().lower()
    if response != 'y':
        print_colored(Colors.YELLOW, "Automation cancelled.")
        sys.exit(0)
    
    # Run all steps
    run_step2_pysam_analysis(output_dir, work_dir, barcodes, search_date)
    print_colored(Colors.GREEN, "\n✓ STEP 2 COMPLETED")
    
    run_step3_fasta_to_ab1(output_dir, barcodes)
    print_colored(Colors.GREEN, "\n✓ STEP 3 COMPLETED")
    
    run_step4_coverage_plots(output_dir, barcodes)
    print_colored(Colors.GREEN, "\n✓ STEP 4 COMPLETED")
    
    run_step5_clone_validation(input_dir, output_dir, barcodes, plasmid_sizes, ecoli_ref)
    print_colored(Colors.GREEN, "\n✓ STEP 5 COMPLETED")
    
    # Final summary
    print_header("AUTOMATION COMPLETE!")
    print_colored(Colors.GREEN + Colors.BOLD, f"✓ Successfully processed {len(barcodes)} barcodes")
    print_colored(Colors.CYAN, f"\nResults saved in: {output_dir}")
    print()
    print_colored(Colors.YELLOW, "Generated files per barcode:")
    print("  • Consensus Sequence Confidence Report.csv")
    print("  • *.ab1 (chromatogram file)")
    print("  • Read Coverage Plot.png")
    print("  • CLONE VALIDATION REPORT.txt")
    print("  • alignment_sorted.bam (alignment file)")
    print()
    print_colored(Colors.BLUE, "Note: Some files may be missing for barcodes that were skipped due to errors.")
    print()

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print_colored(Colors.YELLOW, "\n\nAutomation interrupted by user. Exiting...")
        sys.exit(0)
    except Exception as e:
        print_colored(Colors.RED, f"\n\nFATAL ERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)