import re
import sys
import subprocess
import time

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

def find_atom_entry(psf_file_path, residue_number, residue_name, atom_type):
    """Find atom entry line for a given residue and atom type in the psf file."""
    with open(psf_file_path, 'r') as psf_file:
        for line in psf_file:
            match = re.match(r"^\s*(\d+)\s+\w\s+{}\s+{}\s+{}\s".format(residue_number, residue_name, atom_type), line)
            if match:
                return line.strip()
    return f"{residue_number} {residue_name} {atom_type} Not found"

def find_atom_ids_for_pairs(mr_file_path, psf_file_path):
    """Find formatted atom IDs and entries for each atom pair listed in the .mr file within the .psf file."""
    atom_pairs = parse_mr_file(mr_file_path)
    formatted_output = []

    max_left_width = 0
    atom_entries = []

    for atom1, atom2 in atom_pairs:
        atom1_entry = find_atom_entry(psf_file_path, atom1[0], atom1[1], atom1[2])
        atom2_entry = find_atom_entry(psf_file_path, atom2[0], atom2[1], atom2[2])
        
        atom1_id = atom1_entry.split()[0] if "Not found" not in atom1_entry else "Not found"
        atom2_id = atom2_entry.split()[0] if "Not found" not in atom2_entry else "Not found"

        left_part = f"{atom1_id:<10}{atom2_id:<10} *** {atom1_entry}"
        max_left_width = max(max_left_width, len(left_part))
        atom_entries.append((left_part, atom2_entry))

    # Align "<==>" dynamically
    for left_part, atom2_entry in atom_entries:
        spaces_needed = max_left_width - len(left_part) + 5  # Add buffer space
        formatted_line = f"{left_part}{' ' * spaces_needed}<==>  {atom2_entry} ***"
        formatted_output.append(formatted_line)

    return formatted_output

if len(sys.argv) < 4:
    print("Usage: python3 script.py <mr_file> <psf_file> <dcd_file>")
    sys.exit(1)

mr_file_path = sys.argv[1]
psf_file_path = sys.argv[2]
dcd_file_path = sys.argv[3]

formatted_atom_pairs = find_atom_ids_for_pairs(mr_file_path, psf_file_path)

output_filename = "carmaindist"
with open(output_filename, "w") as output_file:
    for line in formatted_atom_pairs:
        output_file.write(line + "\n")
        print(line)

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
