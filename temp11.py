import re


amino_acids = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL"
}


mr_file = input("Enter the name of the MR file: ")
psf_file = input("Enter the name of the PSF file: ")


ch1 = mr_file.split(".")
ch2 = psf_file.split(".")

if len(ch1) == 2 and ch1[1] == "mr" and len(ch2) == 2 and ch2[1] == "psf":
    with open(mr_file) as file_1, open(psf_file) as file_2:
        
        mr_lines = file_1.readlines()
        psf_lines = file_2.readlines()

    
    psf_dict = {}
    
    
    pattern = re.compile(r"(" + "|".join(amino_acids) + r")\s+([A-Za-z0-9]{1})([A-Za-z0-9]{1})([A-Za-z0-9]{1})")
    
    
    for line in psf_lines:
        match = pattern.search(line)
        if match:
            
            atom_key = f"{match.group(1)} {match.group(2)}{match.group(3)}{match.group(4)}"
            psf_dict[atom_key] = line.strip()

    
    print("\nMatching atoms with assigned numbers from MR file:")
    for line_num, line in enumerate(mr_lines, start=1):
        match = pattern.search(line)
        if match:
            
            mr_atom_key = f"{match.group(1)} {match.group(2)}{match.group(3)}{match.group(4)}"
            if mr_atom_key in psf_dict:
                
                mr_details = line.strip()
                psf_details = psf_dict[mr_atom_key]
                mr_assigned_num = line.split()[0]
                psf_assigned_num = psf_details.split()[0]
                
                print(f"{mr_assigned_num}   {psf_assigned_num}        ***     {mr_details}     <==>     {psf_details}  ***")
            else:
                
                print(f"Atom \"{mr_atom_key}\" not found in line {line_num}: {line.strip()}")

else:
    print("File format error!")
