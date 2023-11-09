import os
import shutil
import glob
import subprocess
from tqdm import tqdm
from rdkit import Chem

def write_gjf(mol, name='NoName', options='', nproc=1, vmem=1):
    """
    Prepares a Gaussian input file for each conformer of mol.
    
    Parameters:
        mol: rdkit.Chem.rdchem.Mol
            RDKit molecule object
        name: str
            Name of the Gaussian job
        options: str
            Additional inputs (e.g. GD3BJ dispersion, polar)
        nproc: int
            Number of CPU cores
        vmem: int
            Amount of RAM in GB
    """
    os.makedirs(name, exist_ok=True)

    files = glob.glob(f"{name}/*.gjf")
    for file in os.listdir(name):
        file_path = os.path.join(name, file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print(f"Error occurred while deleting file {file_path}: {e}")
    
    for q in range(mol.GetNumConformers()):
        file_name = f"{name}/{name}_{q}.gjf"
        with open(file_name, 'w') as f:
            f.write(f"%chk={name}_{q}.chk\n%nprocshared={nproc}\n%mem={vmem}GB\n#p {options}\n\n{name}\n\n0 1\n")
            for i, atom in enumerate(mol.GetAtoms()):
                positions = mol.GetConformer(q).GetAtomPosition(i)
                f.write(f"{atom.GetSymbol()}\t{positions.x}\t{positions.y}\t{positions.z}\n")
            f.write('\n\n')
    return file_name


def run_gjf(g_path, name='NoName'):
    """
    Executes Gaussian for calculations as written by write_gjf function.
    
    Parameters:
        g_path: str
            Path to the Gaussian executable
        name: str
            Name of the Gaussian job
    """
    g_folder = os.path.dirname(g_path)
    g_exec = os.path.basename(g_path)
    
    # Check if the Gaussian executable exists
    if not os.path.isfile(g_path):
        raise FileNotFoundError(f"The Gaussian executable {g_path} does not exist.")
    
    files = glob.glob(f"{name}/*.gjf")
    if not files:
        print(f"No Gaussian input files found in {name}/")
        return
    
    print('Executing Gaussian jobs:')
    for file in tqdm(files):
        temp_gjf_path = os.path.join(g_folder, 'temp.gjf')
        shutil.copy(file, temp_gjf_path)
        
        try:
            subprocess.run([g_path, 'temp.gjf'], cwd=g_folder, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error occurred while running Gaussian on {file}: {e}")
            continue  # Continue with the next file if an error occurs
        
        shutil.copy(os.path.join(g_folder, 'temp.out'), f"{file[:-4]}.out")

def read_energy(file):
    """
    Reads the single-point energy from a Gaussian output file.
    
    Parameters:
        file: str
            Path to the Gaussian output file
    Returns:
        List of energies in Kcal mol-1
    """
    with open(file) as f:
        return [float(line.split(' = ')[1].split(' A.U.')[0]) * 627.5 for line in f if 'SCF Done:' in line]


def read_g_log(name='NoName'):
    """
    Reads the single-point energy from all Gaussian output files in a directory.
    
    Parameters:
        name: str
            Name of the Gaussian job
    Returns:
        List of energies in Kcal mol-1
    """
    files = glob.glob(f"{name}/*.out") or glob.glob(f"{name}/*.log")
    if not files:
        print(f'No Gaussian output file(s) named "{name}" were found')
        return
    
    energies = [read_energy(file)[-1] for file in files]
    print(len(energies))
    return energies

def atomic_number_to_symbol(atomic_number):
    '''
    Simply returns the atomic symbol for a given element number (lazy!)
    '''
    element_symbols = [
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Ni", "Co", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
        "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
        "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
        "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
        "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
    ]
    
    if 1 <= atomic_number <= len(element_symbols):
        return element_symbols[atomic_number - 1]
    else:
        return "Invalid Atomic Number"
    
def get_geometries(file_path):
    '''
    basic function that extracts any and all geometries present in a Gaussian output file.
    
    These can be the output of an optimisation (e.g. all_geometries[-0], the last geometry)
    They might be all geometries in an optimisation of some sort
    
    Args:
    file_path - the file you want to read
    
    Returns:
    all_geometries - the extracted geometry data from the "standard orientation" tables.
    '''

    all_geometries = [] # List to store the extracted lines

    with open(file_path, 'r') as file:
        file_text = file.read()
        match_start = [m.start() for m in re.finditer('Standard orientation',file_text)]
        match_finish = [n.start() for n in re.finditer('Rotational constants',file_text)]

    for m in range(len(match_start)):
        all_geometries.append(file_text[match_start[m]:match_finish[m]])
        
    #print('I found ' + str(len(all_geometries)) + ' geometries in total in ' + file_path + '\n')

    return all_geometries

def generate_coords(extracted_geometry):
    '''
    takes extracted geometry from get_geometries and polishes it so we can reuse it in a new input .gjf file
    
    Args:
    extracted_geometry - our extracted geometry (just a single geometry!) from get_geometries.
    
    Returns
    coords - .gjf formatted coordinates.
    '''
    
    lines = extracted_geometry.strip().split('\n')
    data_lines = lines[5:-1]  # Skipping header and footer lines

    atomic_data = []
    coords = [] # empty array for coordinate data
    # Extract atomic number, element symbol, and coordinates
    for line in data_lines:
        columns = line.split()
        atomic_number = int(columns[1])
        element_symbol = atomic_number_to_symbol(atomic_number)
        x, y, z = float(columns[3]), float(columns[4]), float(columns[5])
        atomic_data.append((element_symbol, x, y, z))

    # Print the extracted data
    for symbol, x, y, z in atomic_data:
        coords.append(f"{symbol}         {x} {y} {z}")
        
    return(coords)

def update_mol_geometry(filename,mol,opt_conf_id = -1):
    '''
    Function that takes an rdkit mol object, generates conformers (if needed), loads geometry from a gaussian ".log"
    file and updates the mol object with the new geometry.
    
    Usage case - a mol object is created and some basic optimisation performed, before using Gaussian to perform a 
    high level DFT calculation. The output of the DFT is read back in to this code so we can perform some analysis.
    
    Args:
    filename - the path and filename of the .log file in question
    mol - the mol object, with our without 3D conformer information
    opt_conf_id - the geometry to return (typically the final one, -1, but could be intermediate geometries)
    
    Returns:
    mol - the updated mol object with the new geometry.
    '''
    
    opt_geom = generate_coords(get_geometries(filename)[opt_conf_id]) # read the file, extract the geometry we want
    numeric_data = []
    for line in opt_geom:
        line_data = line.split()[1:]  # Exclude the first element (atom symbol)
        numeric_values = [float(value) for value in line_data]
        numeric_data.append(numeric_values)

    new_xyz_data = np.array(numeric_data)
    new_mol = mol
    rdDepictor.Compute2DCoords(new_mol,canonOrient=True)
    conf = new_mol.GetConformer()
    for i in range(new_mol.GetNumAtoms()):
        x,y,z = new_xyz_data[i]
        conf.SetAtomPosition(i,Point3D(x,y,z))

    return new_mol
