import re
import sys
import os
import subprocess
import time

# Bidirectional mapping for Markley nomenclature (Markley et al., 1998)
markley_map = {
    # Alanine (ALA)
    "QB": ["HB1", "HB2", "HB3"],  

    # Valine (VAL)
    "QG1": ["HG11", "HG12", "HG13"],  
    "QG2": ["HG21", "HG22", "HG23"],  

    # Leucine (LEU)
    "QD1": ["HD11", "HD12", "HD13"],  
    "QD2": ["HD21", "HD22", "HD23"],  

    # Isoleucine (ILE)
    "QG1": ["HG12", "HG13"],  
    "QG2": ["HG21", "HG22", "HG23"],  
    "QD1": ["HD11", "HD12", "HD13"],  

    # Threonine (THR)
    "QG1": ["HG1"],  
    "HG2": "HG2",  
    "HG3": "HG3",  

    # Serine (SER)
    "QB": ["HB2", "HB3"],  
    "HG": "HG",  

    # Cysteine (CYS)
    "QB": ["HB2", "HB3"],  
    "HG": "HG",  

    # Methionine (MET)
    "QB": ["HB2", "HB3"],  
    "QG": ["HG2", "HG3"],  
    "QE": ["HE1", "HE2", "HE3"],  

    # Phenylalanine (PHE)
    "QB": ["HB2", "HB3"],  
    "QD": ["HD1", "HD2"],  
    "QE": ["HE1", "HE2"],  
    "HZ": "HZ",  

    # Tyrosine (TYR)
    "QB": ["HB2", "HB3"],  
    "QD": ["HD1", "HD2"],  
    "QE": ["HE1", "HE2"],  
    "HH": "HH",  

    # Tryptophan (TRP)
    "QB": ["HB2", "HB3"],  
    "QD1": "HD1",  
    "QE3": "HE3",  
    "QZ2": "HZ2",  
    "QZ3": "HZ3",  
    "QH2": "HH2",  

    # Glutamine (GLN)
    "QB": ["HB2", "HB3"],  
    "QG": ["HG2", "HG3"],  
    "QE": ["HE21", "HE22"],  

    # Glutamate (GLU)
    "QB": ["HB2", "HB3"],  
    "QG": ["HG2", "HG3"],  
    "QD": ["HD21", "HD22"],  

    # Asparagine (ASN)
    "QB": ["HB2", "HB3"],  
    "QD": ["HD21", "HD22"],  

    # Aspartate (ASP)
    "QB": ["HB2", "HB3"],  
    "QD": ["HD21", "HD22"],  

    # Lysine (LYS)
    "QB": ["HB2", "HB3"],  
    "QG": ["HG2", "HG3"],  
    "QD": ["HD2", "HD3"],  
    "QE": ["HE2", "HE3"],  
    "QZ": ["HZ1", "HZ2", "HZ3"],  

    # Arginine (ARG)
    "QB": ["HB2", "HB3"],  
    "QG": ["HG2", "HG3"],  
    "QD": ["HD2", "HD3"],  
    "QH": ["HH11", "HH12", "HH21", "HH22"],  

    # Histidine (HIS)
    "QB": ["HB2", "HB3"],  
    "QD": ["HD1"],  
    "QE": ["HE1"],  

    # Proline (PRO)
    "QB": ["HB2", "HB3"],  
    "QG": ["HG2", "HG3"],  
    "QD": ["HD2", "HD3"],  
}

# Generate reverse mappings (to allow bidirectional lookup)
reverse_markley_map = {}
for k, v in markley_map.items():
    if isinstance(v, list):
        for item in v:
            reverse_markley_map[item] = k
    else:
        reverse_markley_map[v] = k
markley_map.update(reverse_markley_map)


def parse_mr_file(mr_file_path):
    """Extract atom pairs in the form 'residue_number residue_name atom_type' from the .mr file."""
    atom_pairs = []
    with open(mr_file_path, 'r') as mr_file:
        for line in mr_file:
            if "Entry H atom name" in line:
                break
            match = re.findall(r"(\d+)\s+(\w+)\s+(\w+)", line)
            if match and len(match) >= 2:
                atom_pairs.append(((match[0][0], match[0][1], match[0][2]), (match[1][0], match[1][1], match[1][2])))
    return atom_pairs
    
def apply_markley_nomenclature(atom_type):
    """Convert atom names to Markley nomenclature if applicable."""
    return markley_map.get(atom_type, atom_type)
    
def find_atom_entry(psf_file_path, residue_number, residue_name, atom_type):

    """Find atom entry line for a given residue and atom type in the psf file.
    If not found with original name, try mapped names but return original format."""
    original_atom_type = atom_type
    atom_types_to_try = [atom_type]
    
    # Add mapped names to try if original isn't found
    if atom_type in markley_map:
        mapped = markley_map[atom_type]
        if isinstance(mapped, list):
            atom_types_to_try.extend(mapped)
        else:
            atom_types_to_try.append(mapped)
    
    with open(psf_file_path, 'r') as psf_file:
        psf_content = psf_file.readlines()
    
    for current_type in atom_types_to_try:
        for line in psf_content:
            match = re.match(r"^\s*(\d+)\s+\w\s+{}\s+{}\s+{}\s".format(residue_number, residue_name, current_type), line)
            if match:
                # Return the line with original atom type to maintain exact formatting
                parts = line.strip().split()
                parts[4] = original_atom_type
                return " ".join(parts)
    
    return f"{residue_number} {residue_name} {atom_type} Not found"
    
def find_atom_ids_for_pairs(mr_file_path, psf_file_path):
    """Find formatted atom IDs and entries for each atom pair listed in the .mr file within the .psf file."""
    atom_pairs = parse_mr_file(mr_file_path)
    formatted_output = []

    for atom1, atom2 in atom_pairs:
        atom1_mapped = apply_markley_nomenclature(atom1[2])
        atom2_mapped = apply_markley_nomenclature(atom2[2])

        atom1_entry = find_atom_entry(psf_file_path, atom1[0], atom1[1], atom1[2])
        atom2_entry = find_atom_entry(psf_file_path, atom2[0], atom2[1], atom2[2])
        
        # Extract atom IDs
        atom1_id = atom1_entry.split()[0] if "Not found" not in atom1_entry else "Not found"
        atom2_id = atom2_entry.split()[0] if "Not found" not in atom2_entry else "Not found"
        
        # Format atom IDs to be right-aligned in 4 spaces
        atom1_id = f"{atom1_id:>4}"
        atom2_id = f"{atom2_id:>4}"
        
        # Format the entries with fixed spacing
        atom1_entry_parts = atom1_entry.split()
        if len(atom1_entry_parts) >= 5:
            atom1_formatted = f"{atom1_entry_parts[0]:>4} {atom1_entry_parts[1]:<1} {atom1_entry_parts[2]:>4} {atom1_entry_parts[3]:>3} {atom1_entry_parts[4]:>4}"
        else:
            atom1_formatted = atom1_entry
            
        atom2_entry_parts = atom2_entry.split()
        if len(atom2_entry_parts) >= 5:
            atom2_formatted = f"{atom2_entry_parts[0]:>4} {atom2_entry_parts[1]:<1} {atom2_entry_parts[2]:>4} {atom2_entry_parts[3]:>3} {atom2_entry_parts[4]:>4}"
        else:
            atom2_formatted = atom2_entry
        
        # Combine everything with proper spacing
        formatted_line = f"{atom1_id} {atom2_id}        *** {atom1_formatted}  <==> {atom2_formatted} ***"
        formatted_output.append(formatted_line)

    return formatted_output

def identify_files(args):
    """Identify input files by their extensions regardless of order."""
    mr_file = None
    psf_file = None
    dcd_file = None
    
    for arg in args:
        if arg.lower().endswith('.mr'):
            mr_file = arg
        elif arg.lower().endswith('.psf'):
            psf_file = arg
        elif arg.lower().endswith('.dcd'):
            dcd_file = arg
    
    if not mr_file or not psf_file:
        print("Error: Both .mr and .psf files are required.")
        sys.exit(1)
    
    return mr_file, psf_file, dcd_file

def main():

    if len(sys.argv) < 3:
        print("Usage: python3 script.py <mr_file> <psf_file> [<dcd_file>]")
        print("Files can be provided in any order.")
        sys.exit(1)

    mr_file_path, psf_file_path, dcd_file_path = identify_files(sys.argv[1:])

    formatted_atom_pairs = find_atom_ids_for_pairs(mr_file_path, psf_file_path)

    output_filename = "carmaindist"
    with open(output_filename, "w") as output_file:
        for line in formatted_atom_pairs:
            output_file.write(line + "\n")
            print(line)

    if dcd_file_path:
        # Record and display CARMA start time
        start_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        print(f"\nCARMA is starting at {start_time}, please wait...")

        # Run CARMA command
        try:
            carma_command = ["carma", "-v", "-atmid", "ALLID", "-dist", output_filename, psf_file_path, dcd_file_path]
            result = subprocess.run(carma_command, capture_output=True, text=True)
            
            print("\nCARMA has finished execution.")
            print("CARMA Output:\n")
            print(result.stdout)
            print(result.stderr)
        except Exception as e:
            print(f"Error running CARMA: {e}")
    else:
        print("\nNo .dcd file provided - output file created but CARMA was not run.")

if __name__ == "__main__":
    main()
