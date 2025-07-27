"""
NOEmd : a computer program to automate the comparison between experimental NMR data and molecular dynamics simulations

This script processes molecular dynamics simulation data against experimental restraints,
from NMR (NOE) data. It performs the following key functions:

1.  Parses a STAR/.mr file for distance restraints and a PSF file for system topology.
2.  Maps the atom names from the restraint file (e.g., using Markley nomenclature)
    to the corresponding atom IDs in the PSF file.
3.  Handles degenerate restraints (e.g., methyl groups) by identifying all
    possible atom pairs.
4.  Generates an input file for CARMA to calculate distances for these pairs
    across a DCD trajectory.
5.  Invokes CARMA to perform the distance calculations.
6.  Processes the output from CARMA, resolving degeneracies by selecting the
    minimum distance for each restraint group in each frame.
7.  Calculates <r^-3> and <r^-6> time averages for each restraint.
8.  Generates a final report comparing the experimental distances to the
    simulation-derived averages and calculates violations.
"""

import re
import sys
import os
import subprocess
import time
import math
import numpy as np
from collections import defaultdict
import mmap
import datetime
import concurrent.futures
from multiprocessing import cpu_count

# Color codes
class Colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

# Amino acid specific mappings from Markley nomenclature
markley_map = {
    "ALA": {"QB": ["HB1", "HB2", "HB3"]},
    "VAL": {
        "QG1": ["HG11", "HG12", "HG13"],
        "QG2": ["HG21", "HG22", "HG23"]
    },
    "LEU": {
        "QD1": ["HD11", "HD12", "HD13"],
        "QD2": ["HD21", "HD22", "HD23"]
    },
    "ILE": {
        "QG1": ["HG12", "HG13"],
        "QG2": ["HG21", "HG22", "HG23"],
        "QD1": ["HD11", "HD12", "HD13"]
    },
    "THR": {
        "QG1": ["HG1"],
        "HG2": "HG2",
        "QG2": ["HG21", "HG22", "HG23"],
        "HG3": "HG3"
    },
    "SER": {
        "QB": ["HB2", "HB3"],
        "HG": "HG"
    },
    "CYS": {
        "QB": ["HB2", "HB3"],
        "HG": "HG"
    },
    "MET": {
        "QB": ["HB2", "HB3"],
        "QG": ["HG2", "HG3"],
        "QE": ["HE1", "HE2", "HE3"]
    },
    "PHE": {
        "QB": ["HB2", "HB3"],
        "QD": ["HD1", "HD2"],
        "QE": ["HE1", "HE2"],
        "HZ": "HZ"
    },
    "TYR": {
        "QB": ["HB2", "HB3"],
        "QD": ["HD1", "HD2"],
        "QE": ["HE1", "HE2"],
        "HH": "HH"
    },
    "TRP": {
        "QB": ["HB2", "HB3"],
        "QD1": "HD1",
        "QE3": "HE3",
        "QZ2": "HZ2",
        "QZ3": "HZ3",
        "QH2": "HH2"
    },
    "GLN": {
        "QB": ["HB2", "HB3"],
        "QG": ["HG2", "HG3"],
        "QE": ["HE21", "HE22"]
    },
    "GLU": {
        "QB": ["HB2", "HB3"],
        "QG": ["HG2", "HG3"],
        "QD": ["HD21", "HD22"]
    },
    "ASN": {
        "QB": ["HB2", "HB3"],
        "QD": ["HD21", "HD22"]
    },
    "ASP": {
        "QB": ["HB2", "HB3"],
        "QD": ["HD21", "HD22"]
    },
    "LYS": {
        "QB": ["HB2", "HB3"],
        "QG": ["HG2", "HG3"],
        "QD": ["HD2", "HD3"],
        "QE": ["HE2", "HE3"],
        "QZ": ["HZ1", "HZ2", "HZ3"]
    },
    "ARG": {
        "QB": ["HB2", "HB3"],
        "QG": ["HG2", "HG3"],
        "QD": ["HD2", "HD3"],
        "QH": ["HH11", "HH12", "HH21", "HH22"]
    },
    "HIS": {
        "QB": ["HB2", "HB3"],
        "QD": ["HD1"],
        "QE": ["HE1"]
    },
    "PRO": {
        "QB": ["HB2", "HB3"],
        "QG": ["HG2", "HG3"],
        "QD": ["HD2", "HD3"]
    }
}

# Generate reverse mappings
reverse_markley_map = {}
for res_name, mappings in markley_map.items():
    res_reverse_map = {}
    for k, v in mappings.items():
        if isinstance(v, list):
            for item in v:
                res_reverse_map[item] = k
        else:
            res_reverse_map[v] = k
    reverse_markley_map[res_name] = res_reverse_map

def parse_mr_file(mr_file_path):
    """Returns atom pairs AND their experimental distances."""
    entries = []
    entry_pattern = re.compile(
        r"(\d+)\s+(\w+)\s+(\w+)\s+(\d+)\s+(\w+)\s+(\w+)\s+([\d.]+)"
    )
    with open(mr_file_path, 'r') as mr_file:
        for line in mr_file:
            if line.startswith('*') or not line.strip():
                continue
            match = entry_pattern.search(line)
            if match:
                atom1 = (match.group(1), match.group(2), match.group(3))
                atom2 = (match.group(4), match.group(5), match.group(6))
                distance = float(match.group(7))
                entries.append((atom1, atom2, distance))
    return entries

def find_all_atom_entries(psf_file_path, residue_number, residue_name, atom_type):
    """Optimized PSF file parser with precompiled regex"""
    entries = []
    atom_types_to_try = {atom_type}
    
    if residue_name in markley_map and atom_type in markley_map[residue_name]:
        mapped = markley_map[residue_name][atom_type]
        if isinstance(mapped, list):
            atom_types_to_try.update(mapped)
        else:
            atom_types_to_try.add(mapped)
    
    if residue_name in reverse_markley_map and atom_type in reverse_markley_map[residue_name]:
        atom_types_to_try.add(reverse_markley_map[residue_name][atom_type])
    
    pattern = re.compile(rf"^\s*(\d+)\s+\w\s+{residue_number}\s+{residue_name}\s+({'|'.join(atom_types_to_try)})\s")
    
    with open(psf_file_path, 'r') as psf_file:
        for line in psf_file:
            match = pattern.match(line)
            if match:
                parts = line.strip().split()
                parts[4] = atom_type  # Keep original atom name
                entries.append(" ".join(parts))
    
    return entries if entries else [f"???? ? {residue_number:>4} {residue_name:>3} {atom_type:>4}"]

def find_atom_ids_for_pairs(mr_file_path, psf_file_path):
    """Returns formatted pairs with their experimental distances."""
    entries = parse_mr_file(mr_file_path)
    formatted_output = []
    seen_pairs = set()

    # Build cache for all unique atoms
    unique_atoms = set()
    for atom1, atom2, _ in entries:
        unique_atoms.add(atom1)
        unique_atoms.add(atom2)
    
    atom_entries_cache = {}
    for atom in unique_atoms:
        atom_entries = find_all_atom_entries(psf_file_path, *atom)
        if not any("????" in entry for entry in atom_entries):
            atom_entries_cache[atom] = atom_entries

    # Process pairs with distances
    for atom1, atom2, exp_distance in entries:
        if atom1 not in atom_entries_cache or atom2 not in atom_entries_cache:
            continue
            
        atom1_entries = atom_entries_cache[atom1]
        atom2_entries = atom_entries_cache[atom2]

        for a1_entry in atom1_entries:
            for a2_entry in atom2_entries:
                atom1_id = a1_entry.split()[0]
                atom2_id = a2_entry.split()[0]
                
                pair_id = tuple(sorted((atom1_id, atom2_id)))
                if pair_id in seen_pairs:
                    continue
                seen_pairs.add(pair_id)
                
                a1_parts = a1_entry.split()
                a2_parts = a2_entry.split()
                
                is_degenerate = (len(atom1_entries) > 1 or len(atom2_entries) > 1)
                
                formatted_line = (
                    f"{atom1_id:>4} {atom2_id:>4} {exp_distance:6.2f} *** "
                    f"{a1_parts[0]:>4} {a1_parts[1]:<1} {a1_parts[2]:>4} {a1_parts[3]:>3} {atom1[2]:>4} <==> "
                    f"{a2_parts[0]:>4} {a2_parts[1]:<1} {a2_parts[2]:>4} {a2_parts[3]:>3} {atom2[2]:>4}"
                    f"{' (Degenerate)' if is_degenerate else ''}"
                )
                formatted_output.append(formatted_line)

    return formatted_output, len(formatted_output) == 0

def resolve_and_filter_degeneracies(formatted_atom_pairs, carma_file='carma.distances'):
    """Resolve degeneracies with exact spacing: 6 spaces for index, 4 spaces between values"""
    try:
        # First pass: group columns by their restraint pattern
        column_groups = defaultdict(list)
        
        # Create unique identifiers for each restraint group
        for col_idx, pair_desc in enumerate(formatted_atom_pairs):
            parts = pair_desc.split("***")[1].split("<==>")
            left_psf = parts[0].strip().split()
            right_psf_part = parts[1].strip()
            right_psf = right_psf_part.split(" (Degenerate)")[0].strip().split()

            # Ensure both left and right parts have enough elements to prevent errors
            if len(left_psf) < 5 or len(right_psf) < 5:
                continue
            
            # Correct key generation to match other functions
            k1 = (left_psf[2], left_psf[3], left_psf[4])  # ResNum, ResName, AtomName_MR
            k2 = (right_psf[2], right_psf[3], right_psf[4])
            group_key = tuple(sorted((k1, k2)))
            
            column_groups[group_key].append(col_idx)
        
        # Get total number of frames for progress reporting
        total_frames = 0
        with open(carma_file, 'r') as f:
            for line in f:
                if line.strip():
                    total_frames += 1
        
        if total_frames == 0:
            print(f"{Colors.FAIL}ERROR: Empty distance file{Colors.ENDC}")
            return None
        
        print(f"{Colors.OKBLUE}Processing {total_frames:,} frames...{Colors.ENDC}")
        
        # Process file with exact spacing
        temp_file = carma_file + ".tmp"
        processed_frames = 0
        start_time = time.time()
        
        with open(carma_file, 'r') as infile, open(temp_file, 'w') as outfile:
            for line in infile:
                if not line.strip():
                    continue
                
                parts = line.strip().split()
                frame_idx = parts[0]
                distances = parts[1:]
                
                new_distances = []
                # For each group, find the smallest distance
                for columns in column_groups.values():
                    if not columns:
                        continue
                        
                    # Get all distances for this group in current frame
                    group_distances = []
                    for col in columns:
                        if col < len(distances):
                            try:
                                group_distances.append(float(distances[col]))
                            except (ValueError, IndexError):
                                continue
                    
                    if group_distances:
                        # Keep the smallest distance in this group
                        min_distance = min(group_distances)
                        new_distances.append(f"{min_distance:.4f}")
                
                # Format with exact spacing:
                # 6 leading spaces for frame index
                # 4 spaces between each value
                new_line = f"      {frame_idx}"  # 6 spaces before index
                for dist in new_distances:
                    new_line += f"    {dist}"  # 4 spaces before each value
                
                outfile.write(new_line + "\n")
                
                processed_frames += 1
                
                # Update progress
                if processed_frames % 100 == 0 or processed_frames == total_frames:
                    elapsed = time.time() - start_time
                    remaining = (elapsed / processed_frames) * (total_frames - processed_frames) if processed_frames > 0 else 0
                    print(
                        f"\r{Colors.OKBLUE}Progress: {processed_frames}/{total_frames} "
                        f"({processed_frames/total_frames*100:.1f}%) | "
                        f"Elapsed: {elapsed:.1f}s | "
                        f"Remaining: {remaining:.1f}s{Colors.ENDC}",
                        end="", flush=True
                    )
        
        print()  
        
        # Replace original file
        os.replace(temp_file, carma_file)
        
        # Print summary
        num_groups = len(column_groups)
        print(f"\n{Colors.OKGREEN}Filtered {carma_file} - kept {num_groups} distance groups (minimum distance per group){Colors.ENDC}")
        
        # Return the actual keys
        return list(column_groups.keys())
        
    except Exception as e:
        print(f"\n{Colors.FAIL}ERROR in resolve_and_filter_degeneracies: {str(e)}{Colors.ENDC}")
        return None

def process_noe_data(formatted_atom_pairs, carma_file='carma.distances', filtered_column_keys=None):
    """Process NOE data from the (potentially filtered) carma.distances file.
    formatted_atom_pairs: The original list from find_atom_ids_for_pairs. Used to map MR restraints.
    filtered_column_keys: The list of unique restraint keys that define the columns in the *filtered* carma.distances file.
    """
    if not os.path.exists(carma_file):
        print(f"{Colors.FAIL}ERROR: {carma_file} not found for NOE processing.{Colors.ENDC}")
        return None
    if not formatted_atom_pairs:
        print(f"{Colors.FAIL}ERROR: No formatted_atom_pairs for NOE processing.{Colors.ENDC}")
        return None
    if filtered_column_keys is None: # If carma.distances wasn't filtered (e.g. resolve_and_filter_degeneracies was skipped or failed)
        print(f"{Colors.WARNING}Assuming carma.distances columns map directly to formatted_atom_pairs.{Colors.ENDC}")
        # Fallback: create effective filtered_column_keys and mapping if not provided
        # This path assumes carma.distances columns correspond to the first occurrence of each unique restraint type
        # as defined by the grouping logic in resolve_and_filter_degeneracies.
        temp_col_groups = defaultdict(list)
        temp_filtered_keys = []
        for idx, pair_desc in enumerate(formatted_atom_pairs):
            if "ERROR PARSING" in pair_desc: continue
            parts = pair_desc.split("***")[1].split("<==>")
            left_psf = parts[0].strip().split()
            right_psf = parts[1].strip().split(" (Degenerate)")[0].strip().split()
            if len(left_psf) < 5 or len(right_psf) < 5: continue
            k1 = (left_psf[2], left_psf[3], left_psf[4])
            k2 = (right_psf[2], right_psf[3], right_psf[4])
            gkey = tuple(sorted((k1, k2)))
            if gkey not in temp_col_groups:
                temp_filtered_keys.append(gkey)
            temp_col_groups[gkey].append(idx)
        filtered_column_keys = temp_filtered_keys


    num_expected_cols = len(filtered_column_keys)
    sum_r3 = np.zeros(num_expected_cols)
    sum_r6 = np.zeros(num_expected_cols)
    n_frames = 0

    print(f"{Colors.OKBLUE}Processing NOE data from {carma_file} with {num_expected_cols} expected restraint columns.{Colors.ENDC}")

    try:
        with open(carma_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line_strip = line.strip()
                if not line_strip or line_strip.startswith('#'):
                    continue 
                
                parts = line_strip.split()
                # Frame_idx is parts[0]
                dist_values_from_file = parts[1:]

                if len(dist_values_from_file) != num_expected_cols:
                    print(f"{Colors.WARNING}Line {line_num} in {carma_file}: Expected {num_expected_cols} distance values, found {len(dist_values_from_file)}. Skipping line.{Colors.ENDC}")
                    continue
                
                for i in range(num_expected_cols):
                    try:
                        d = float(dist_values_from_file[i])
                        # Add to sums (add epsilon to avoid division by zero if d is exactly 0)
                        sum_r3[i] += 1.0 / ((d + 1e-9) ** 3) 
                        sum_r6[i] += 1.0 / ((d + 1e-9) ** 6)
                    except ValueError:
                        print(f"{Colors.WARNING}Line {line_num}, Col {i+1}: Error parsing distance '{dist_values_from_file[i]}'. Skipping value.{Colors.ENDC}")
                        
                        continue 
                n_frames += 1
                if n_frames % 1000 == 0:
                    print(f"\r{Colors.OKCYAN}Processed {n_frames} frames for NOE averaging...{Colors.ENDC}", end="", flush=True)
        
        print(f"\n{Colors.OKGREEN}Loaded NOE data from {n_frames:,} frames.{Colors.ENDC}")
        if n_frames == 0:
            print(f"{Colors.FAIL}No frames processed from {carma_file}. Cannot calculate NOE averages.{Colors.ENDC}")
            return None

    except Exception as e:
        print(f"{Colors.FAIL}ERROR reading or processing {carma_file}: {e}{Colors.ENDC}")
        import traceback
        traceback.print_exc()
        return None

    # Create mapping from original formatted_atom_pairs index to the index in noe_data/filtered_column_keys
    # This is crucial for run_noe_averaging to find the correct averaged data.
    column_mapping_original_to_filtered = {}
    for original_idx, pair_desc in enumerate(formatted_atom_pairs):
        if "ERROR PARSING" in pair_desc: continue
        parts = pair_desc.split("***")[1].split("<==>")
        left_psf = parts[0].strip().split()
        right_psf = parts[1].strip().split(" (Degenerate)")[0].strip().split()
        if len(left_psf) < 5 or len(right_psf) < 5: continue
        
        atom1_key_tuple = (left_psf[2], left_psf[3], left_psf[4]) # ResNum, ResName, AtomName_MR
        atom2_key_tuple = (right_psf[2], right_psf[3], right_psf[4])
        current_restraint_key = tuple(sorted((atom1_key_tuple, atom2_key_tuple)))
        
        if current_restraint_key in filtered_column_keys:
            filtered_idx = filtered_column_keys.index(current_restraint_key)
            column_mapping_original_to_filtered[original_idx] = filtered_idx


    noe_data_results = {
        'nof_values': n_frames,
        'pairs_data': [], # Each entry corresponds to a column in filtered_column_keys
        'filtered_keys_definition': filtered_column_keys, # The definition of each column
        'column_mapping_original_to_filtered': column_mapping_original_to_filtered
    }

    for i in range(num_expected_cols):
        pair_summary = {
            'filtered_col_index': i,
            'restraint_definition_key': filtered_column_keys[i], # The tuple defining this restraint group
            'sum_r3': sum_r3[i],
            'sum_r6': sum_r6[i]
        }
        noe_data_results['pairs_data'].append(pair_summary)

    return noe_data_results


def run_noe_averaging(noe_data, formatted_atom_pairs):
    """Generate NOE analysis report.
    noe_data: Contains sums (from process_noe_data) and mapping.
    formatted_atom_pairs: The original list from find_atom_ids_for_pairs, used to define unique display lines.
    """
    if not noe_data or not noe_data.get('pairs_data'):
        print(f"\n{Colors.FAIL}ERROR: No NOE data available to analyze.{Colors.ENDC}")
        return

    print(f"\n{Colors.HEADER}NOE Analysis Results:{Colors.ENDC}")
    print(f"{Colors.OKBLUE}{'=' * 100}{Colors.ENDC}") 
    print(f"{Colors.BOLD}  EXP     <r^-3>    <r^-6>   Deviation    Violation          Atoms Involved{Colors.ENDC}") 
    print(f"{Colors.OKBLUE}{'=' * 100}{Colors.ENDC}") 

    unique_display_entries = {} 

    if not formatted_atom_pairs:
        print(f"{Colors.WARNING}formatted_atom_pairs is empty. Cannot generate detailed NOE analysis table.{Colors.ENDC}")
        return

    for original_idx, pair_desc_full_line in enumerate(formatted_atom_pairs):
        if "ERROR PARSING" in pair_desc_full_line:
            continue

        try:
            initial_parts = pair_desc_full_line.split("***")
            if len(initial_parts) < 2: 
                print(f"{Colors.WARNING}Malformed formatted_atom_pair line (missing '***'): '{pair_desc_full_line}'. Skipping.{Colors.ENDC}")
                continue

            dist_part_split = initial_parts[0].strip().split()
            if len(dist_part_split) < 3: 
                print(f"{Colors.WARNING}Malformed distance part in formatted_atom_pair line: '{initial_parts[0].strip()}'. Skipping.{Colors.ENDC}")
                continue
            distance_str = dist_part_split[2] 
            
            interaction_part_full = initial_parts[1].strip() 
            is_line_marked_degenerate = "(Degenerate)" in interaction_part_full # Still track this for internal logic if needed later
            
            interaction_halves = interaction_part_full.split("<==>")
            if len(interaction_halves) < 2:
                print(f"{Colors.WARNING}Malformed interaction part (missing '<==>'): '{interaction_part_full}'. Skipping.{Colors.ENDC}")
                continue

            left_half_full_desc = interaction_halves[0].strip() 
            right_half_full_desc_maybe_degen = interaction_halves[1].strip()
            right_half_full_desc = right_half_full_desc_maybe_degen.split(" (Degenerate)")[0].strip()

            l_parts = left_half_full_desc.split() # [ID, CHAIN, RESSEQ, RESNAME, ATOMNAME_MR]
            r_parts = right_half_full_desc.split()

            if len(l_parts) < 5 or len(r_parts) < 5:
                print(f"{Colors.WARNING}Could not parse atom details from: '{interaction_part_full}'. Expected 5 parts per side. Got L:{len(l_parts)}, R:{len(r_parts)}. Skipping.{Colors.ENDC}")
                continue

            # Key components for display grouping: (ResSeq, ResName, AtomName_MR)
            atom1_key_details = (l_parts[2], l_parts[3], l_parts[4]) 
            atom2_key_details = (r_parts[2], r_parts[3], r_parts[4])
            display_key = (distance_str, tuple(sorted((atom1_key_details, atom2_key_details))))

            # Construct the shortened atom description for display
            # Format: RESSEQ RESNAME ATOMNAME_MR
            atom1_desc_short = f"{l_parts[2]:>4} {l_parts[3]:>3} {l_parts[4]:>4}" 
            atom2_desc_short = f"{r_parts[2]:>4} {r_parts[3]:>3} {r_parts[4]:>4}"
            final_atom_details_for_print = f"{atom1_desc_short} <==> {atom2_desc_short}"


            if display_key not in unique_display_entries:
                unique_display_entries[display_key] = {
                    'atom_details_for_print': final_atom_details_for_print, # Store the shortened version
                    'is_group_degenerate': is_line_marked_degenerate, # Store if *any* original line for this group was degenerate
                    'first_original_formatted_pair_idx': original_idx 
                }
            else:
                if is_line_marked_degenerate: # If any line contributing to this unique display is degenerate, mark the group
                    unique_display_entries[display_key]['is_group_degenerate'] = True
        
        except Exception as e:
            print(f"{Colors.FAIL}Error parsing formatted_atom_pair line for display: '{pair_desc_full_line}'. Error: {e}{Colors.ENDC}")
            continue

    stats = {
        'total_displayed_noes': 0,
        'violation_count': 0,
        'total_violation_sum': 0.0,
        'max_violation': 0.0,
        'by_severity': defaultdict(int)
    }

    for display_key, data in unique_display_entries.items():
        exp_distance = float(display_key[0]) 
        line_to_print_atom_part = data['atom_details_for_print'] 
        # is_group_degenerate = data['is_group_degenerate'] # Available if needed, but not printed as "(Degenerate)"
        
        representative_original_idx = data['first_original_formatted_pair_idx'] 

        avg_r3_dist = 999.0 
        avg_r6_dist = 999.0 
        
        if noe_data.get('column_mapping_original_to_filtered') is None:
            print(f"{Colors.WARNING}column_mapping_original_to_filtered is missing in noe_data. Cannot map restraints for averaging.{Colors.ENDC}")
        elif representative_original_idx in noe_data['column_mapping_original_to_filtered']:
            idx_in_noe_pairs_data = noe_data['column_mapping_original_to_filtered'][representative_original_idx]
            
            if idx_in_noe_pairs_data < len(noe_data['pairs_data']):
                noe_sums = noe_data['pairs_data'][idx_in_noe_pairs_data]
                current_sum_r3 = noe_sums['sum_r3']
                current_sum_r6 = noe_sums['sum_r6']

                if noe_data['nof_values'] > 0:
                    avg_r3_dist = (current_sum_r3 / noe_data['nof_values'])**(-1/3) if current_sum_r3 > 1e-9 else 999.0 
                    avg_r6_dist = (current_sum_r6 / noe_data['nof_values'])**(-1/6) if current_sum_r6 > 1e-9 else 999.0
            else:
                print(f"{Colors.WARNING}Mapped index {idx_in_noe_pairs_data} out of bounds for noe_data.pairs_data (len {len(noe_data['pairs_data'])}). Restraint key for display: {display_key}{Colors.ENDC}")
        else:
            print(f"{Colors.WARNING}Original index {representative_original_idx} (key: {display_key}) not found in noe_data's column_mapping. Cannot retrieve averaged distances.{Colors.ENDC}")

        deviation = avg_r6_dist - exp_distance if avg_r6_dist != 999.0 else 999.0
        violation = max(0, deviation) if deviation != 999.0 else 999.0
        
        stats['total_displayed_noes'] += 1
        if violation > 1e-3 and violation != 999.0: 
            stats['violation_count'] += 1
            stats['total_violation_sum'] += violation
            stats['max_violation'] = max(stats['max_violation'], violation)
            
            if violation > 3.0: stats['by_severity']['>3.0Å'] += 1
            elif violation > 2.0: stats['by_severity']['>2.0Å'] += 1
            elif violation > 1.0: stats['by_severity']['>1.0Å'] += 1
            elif violation > 0.5: stats['by_severity']['>0.5Å'] += 1 
            else: stats['by_severity']['0-0.5Å'] +=1

        v_color = Colors.ENDC
        if violation == 999.0: v_color = Colors.FAIL
        elif violation > 3.0: v_color = Colors.FAIL
        elif violation > 1.0: v_color = Colors.WARNING
        elif violation > 1e-3 : v_color = Colors.OKGREEN 
        
        # Print the line with the modified atom part and no (Degenerate) flag
        print(
            f" {exp_distance:5.2f}    {avg_r3_dist:6.3f}    {avg_r6_dist:6.3f}    "
            f"{deviation:7.3f}    {v_color}{violation:7.3f}{Colors.ENDC}    "
            f"{line_to_print_atom_part}" 
        )

    print(f"\n{Colors.OKBLUE}{'=' * 100}{Colors.ENDC}") # Adjusted width
    if stats['total_displayed_noes'] > 0:
        violation_pct = (stats['violation_count'] / stats['total_displayed_noes']) * 100 if stats['total_displayed_noes'] > 0 else 0
        avg_viol = stats['total_violation_sum'] / stats['violation_count'] if stats['violation_count'] > 0 else 0.0
        
        summary_color = Colors.OKGREEN
        if violation_pct > 30: summary_color = Colors.FAIL
        elif violation_pct > 15: summary_color = Colors.WARNING

        print(f"{Colors.BOLD}SUMMARY: {stats['total_displayed_noes']} unique NOE restraints displayed.{Colors.ENDC}")
        print(f"{summary_color}Violations (>0.001Å): {stats['violation_count']} ({violation_pct:.1f}%){Colors.ENDC}")
        print(f"Average violation (for violated restraints): {avg_viol:.3f} Å")
        print(f"Maximum violation: {stats['max_violation']:.3f} Å")
        print("Violation severity distribution:")
        severity_order = ['0-0.5Å', '>0.5Å', '>1.0Å', '>2.0Å', '>3.0Å']
        sorted_severity_keys = [key for key in severity_order if key in stats['by_severity']]
        for severity_range in sorted_severity_keys: 
            print(f"  {severity_range}: {stats['by_severity'][severity_range]} restraints")
        for sev_key, count in stats['by_severity'].items():
            if sev_key not in sorted_severity_keys:
                 print(f"  {sev_key}: {count} restraints")
    else:
        print(f"{Colors.WARNING}No NOEs were displayed in the analysis.{Colors.ENDC}")


def identify_files(args):
    """Identify input files by extension"""
    mr_file = psf_file = dcd_file = None
    for arg in args:
        if arg.lower().endswith('.mr'):
            mr_file = arg
        elif arg.lower().endswith('.psf'):
            psf_file = arg
        elif arg.lower().endswith('.dcd'):
            dcd_file = arg
    if not mr_file or not psf_file:
        print(f"{Colors.FAIL}Error: Both .mr and .psf files are required.{Colors.ENDC}")
        sys.exit(1)
    return mr_file, psf_file, dcd_file

def main():
    if len(sys.argv) < 3:
        print(f"{Colors.FAIL}Usage: python3 {sys.argv[0]} <mr_file> <psf_file> [<dcd_file>]{Colors.ENDC}")
        sys.exit(1)

    mr_file_path, psf_file_path, dcd_file_path = identify_files(sys.argv[1:])
    
    print(f"{Colors.OKCYAN}Processing MR file: {mr_file_path}, PSF file: {psf_file_path}{Colors.ENDC}")
    if dcd_file_path: print(f"{Colors.OKCYAN}DCD file provided: {dcd_file_path}{Colors.ENDC}")

    # First print the atom pair table and save to file
    # formatted_atom_pairs still contains the exp_distance for internal use by other functions.
    formatted_atom_pairs, has_unknown = find_atom_ids_for_pairs(mr_file_path, psf_file_path)
    
    print(f"\n{Colors.OKBLUE}Atom Pair Table (from MR and PSF):{'='*43}{Colors.ENDC}")
    if not formatted_atom_pairs:
        print(f"{Colors.WARNING}No atom pairs could be formatted. Check MR/PSF files and parsing logic.{Colors.ENDC}")
    
    # Print loop for the first table (console) - remove distance
    for line in formatted_atom_pairs:
        line_to_print_console = line
        if "***" in line:
            parts = line.split("***", 1)
            first_part_elements = parts[0].strip().split() # [id1, id2, distance_val_str]
            if len(first_part_elements) == 3: # If it has the distance
                 # Format: id1 id2 (8 spaces for distance placeholder) *** details
                line_to_print_console = f"{first_part_elements[0]:>4} {first_part_elements[1]:>4}{'':4s}***{parts[1]}"
            # If first_part_elements is not 3 (e.g. error line), print as is.
        
        if "????" in line or "ERROR" in line.upper(): 
            print(f"{Colors.FAIL}{line_to_print_console}{Colors.ENDC}")
        else:
            print(line_to_print_console)
    print(f"{Colors.OKBLUE}{'='*77}{Colors.ENDC}")
    
    table1_out_path = "table1.out"
    with open(table1_out_path, "w") as table_file:
        table_file.write("="*85 + "\n")
        # Write to file, visually removing distance
        for line in formatted_atom_pairs:
            line_to_write_file = line
            if "***" in line:
                parts = line.split("***", 1)
                first_part_elements = parts[0].strip().split()
                if len(first_part_elements) == 3:
                    line_to_write_file = f"{first_part_elements[0]:>4} {first_part_elements[1]:>4}{'':4s}***{parts[1]}"

            table_file.write(line_to_write_file + "\n")
        table_file.write("="*85 + "\n")
    print(f"\n{Colors.OKGREEN}Full atom pair table saved to {table1_out_path}{Colors.ENDC}")
    
    if has_unknown:
        print(f"\n{Colors.FAIL}Unknown atoms or parsing errors encountered during PSF lookup. Aborting before CARMA.{Colors.ENDC}")
        sys.exit(1)
    if not formatted_atom_pairs:
        print(f"\n{Colors.FAIL}No valid atom pairs generated. Aborting before CARMA.{Colors.ENDC}")
        sys.exit(1)


    # Write atom IDs for CARMA distance calculation (only the ID pairs)
    # CARMA input format for -dist: atomid1 atomid2
    # The formatted_atom_pairs lines start with these IDs.
    carma_input_dist_file = "carmainput_for_distcalc.txt"
    try:
        with open(carma_input_dist_file, "w") as f_carma_in:
            # Create a set to write only unique atom ID pairs for CARMA
            unique_id_pairs_for_carma = set()
            for line in formatted_atom_pairs:
                if "ERROR" in line.upper(): continue # Skip error lines
                parts = line.split()
                if len(parts) >= 2: # Check if there are at least two parts (ID1, ID2)
                    atom_id1, atom_id2 = parts[0], parts[1]
                    # Ensure they are numbers and form a sorted tuple for uniqueness
                    try:
                        # Ensure atom_id1 and atom_id2 are not '????'
                        if atom_id1 == "????" or atom_id2 == "????":
                            continue
                        id_pair = tuple(sorted((int(atom_id1), int(atom_id2))))
                        if id_pair not in unique_id_pairs_for_carma:
                            f_carma_in.write(f"{atom_id1} {atom_id2}\n")
                            unique_id_pairs_for_carma.add(id_pair)
                    except ValueError: # Handles cases where atom_id1 or atom_id2 are not integers
                        # This can happen if the line was an error message not starting with IDs
                        # print(f"{Colors.WARNING}Skipping non-integer atom IDs for CARMA input: {atom_id1}, {atom_id2} from line '{line}'{Colors.ENDC}")
                        pass # Silently skip if not valid IDs



        carma_distances_file = "carma.distances"
        # filtered_keys_for_noe will hold the definitions of columns in carma.distances after potential filtering
        # It's a list of keys, where each key defines a restraint group.
        filtered_keys_for_noe = None 

        if dcd_file_path:
            run_carma_calc = True # Assume we need to run CARMA if DCD is provided
            if os.path.exists(carma_distances_file):
                resp = input(f"\n{Colors.WARNING}File '{carma_distances_file}' already exists. Regenerate with CARMA? (y/n, default y): {Colors.ENDC}").lower().strip()
                if resp in ('n', 'no'):
                    run_carma_calc = False
                    print(f"{Colors.OKCYAN}Skipping CARMA calculation, using existing '{carma_distances_file}'.{Colors.ENDC}")
            
            if run_carma_calc:
                print(f"\n{Colors.OKBLUE}Starting CARMA calculation at {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}{Colors.ENDC}")
                # Using the simplified carma_input_dist_file
                carma_cmd = ["carma", "-v", "-atmid" ,"ALLID","-dist","table1.out", psf_file_path, dcd_file_path]
                print(f"Executing: {' '.join(carma_cmd)}")
                try:
                    subprocess.run(carma_cmd, check=True) # check=True will raise CalledProcessError if CARMA fails
                    print(f"{Colors.OKGREEN}CARMA calculation finished successfully at {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}{Colors.ENDC}")
                except subprocess.CalledProcessError as e:
                    print(f"{Colors.FAIL}CARMA execution failed with error: {e}{Colors.ENDC}")
                    print(f"{Colors.FAIL}Cannot proceed with NOE analysis without {carma_distances_file}.{Colors.ENDC}")
                    sys.exit(1) # Exit if CARMA fails and was supposed to run
                except FileNotFoundError:
                    print(f"{Colors.FAIL}CARMA command not found. Please ensure CARMA is installed and in your PATH.{Colors.ENDC}")
                    sys.exit(1)


        # Proceed with NOE analysis if carma.distances exists (either newly created or pre-existing)
        if os.path.exists(carma_distances_file):
            print(f"\n{Colors.OKBLUE}Starting NOE processing using '{carma_distances_file}' at {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}{Colors.ENDC}")
            
            perform_filtering = True # Default to True: filter the carma.distances file
            
            # Ask user about filtering carma.distances, even if it was just generated
            # This allows using a pre-filtered file or re-filtering one.
            filter_response = input(f"\n{Colors.WARNING}Do you want to (re-)filter '{carma_distances_file}' for NOE processing (recommended for degenerate restraints)? (y/n, default y): {Colors.ENDC}").lower().strip()
            
            if filter_response in ('n', 'no'):
                perform_filtering = False
                print(f"{Colors.OKCYAN}Using existing '{carma_distances_file}' as-is without re-filtering. Column structure will be deduced.{Colors.ENDC}")
                filtered_keys_for_noe = None # process_noe_data will deduce keys
            else: # Default is yes, or explicit yes
                print(f"{Colors.OKCYAN}Proceeding to filter '{carma_distances_file}'.{Colors.ENDC}")
                # Resolve degeneracies in carma.distances file and get definitions of the filtered columns
                # The `formatted_atom_pairs` is used here to understand which original restraints map to which columns
                # before filtering.
                filtered_keys_for_noe = resolve_and_filter_degeneracies(formatted_atom_pairs, carma_distances_file)
                if filtered_keys_for_noe is None: # Filtering was chosen but failed
                    print(f"{Colors.FAIL}Filtering of '{carma_distances_file}' failed. Cannot proceed with NOE analysis.{Colors.ENDC}")
                    sys.exit(1) # Exit if filtering fails
            
            # Now, process_noe_data will use filtered_keys_for_noe if filtering was done,
            # or it will be None if user skipped filtering, and process_noe_data will deduce.
            noe_data_struct = process_noe_data(formatted_atom_pairs, carma_distances_file, filtered_keys_for_noe)
            
            if noe_data_struct: # Check if process_noe_data was successful
                run_noe_averaging(noe_data_struct, formatted_atom_pairs)
            else:
                print(f"{Colors.FAIL}NOE data processing failed. Cannot run averaging.{Colors.ENDC}")
            
            print(f"\n{Colors.OKGREEN}NOE processing completed at {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}{Colors.ENDC}")

        elif dcd_file_path : # DCD was given but carma.distances still doesn't exist (e.g. CARMA failed silently or was skipped and file never existed)
             print(f"{Colors.FAIL}'{carma_distances_file}' is missing. CARMA might have failed or was skipped. Cannot proceed.{Colors.ENDC}")
        else: # No DCD, no carma.distances
             print(f"{Colors.WARNING}No DCD file provided and no existing '{carma_distances_file}'. Cannot perform NOE analysis.{Colors.ENDC}")

    finally:
        if os.path.exists(carma_input_dist_file):
            try:
                os.remove(carma_input_dist_file)
                # print(f"{Colors.OKCYAN}Cleaned up temporary file: {carma_input_dist_file}{Colors.ENDC}")
            except OSError as e:
                print(f"{Colors.WARNING}Could not remove temporary file {carma_input_dist_file}: {e}{Colors.ENDC}")


if __name__ == "__main__":
    main()
