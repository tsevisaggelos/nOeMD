import re
import sys

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
    """Find atom entry line for a given residue and atom type in the psf file with specified format."""
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

    for atom1, atom2 in atom_pairs:
        # Find the atom entries in the psf file for each atom in the pair
        atom1_entry = find_atom_entry(psf_file_path, atom1[0], atom1[1], atom1[2])
        atom2_entry = find_atom_entry(psf_file_path, atom2[0], atom2[1], atom2[2])
        
        # Extract the atom IDs from each entry if they exist, otherwise use "Not found" message
        atom1_id = atom1_entry.split()[0] if "Not found" not in atom1_entry else "Not found"
        atom2_id = atom2_entry.split()[0] if "Not found" not in atom2_entry else "Not found"

        # Format the output with fixed-width fields for consistent alignment
        formatted_line = (
            f"{atom1_id:<8} {atom2_id:<8} *** {atom1_entry:<60} <==> {atom2_entry:<60} ***"
        )
        formatted_output.append(formatted_line)

    return formatted_output


# Check for command-line arguments or prompt for input
if len(sys.argv) >= 3:
    file1_path = sys.argv[1]
    file2_path = sys.argv[2]
elif len(sys.argv) == 2:
    print("Please provide both .mr and .psf files.")
    sys.exit(1)
else:
    file1_path = input("Please enter the path to the .mr or .psf file: ")
    file2_path = input("Please enter the path to the other file (.mr or .psf): ")

# Determine file types based on extensions
if file1_path.endswith(".mr") and file2_path.endswith(".psf"):
    mr_file_path = file1_path
    psf_file_path = file2_path
elif file1_path.endswith(".psf") and file2_path.endswith(".mr"):
    mr_file_path = file2_path
    psf_file_path = file1_path
else:
    print("Error: Please provide one .mr file and one .psf file.")
    sys.exit(1)

# Run the function and display results
formatted_atom_pairs = find_atom_ids_for_pairs(mr_file_path, psf_file_path)
for line in formatted_atom_pairs:
    print(line)

