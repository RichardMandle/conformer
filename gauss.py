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


