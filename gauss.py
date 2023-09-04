import os
import shutil
import glob
import subprocess

from tqdm import tqdm

from rdkit import Chem


def write_gjf(mol,name='NoName',options='',nproc=1,vmem=1):
    """
    This function prepares a Gaussian input file for each conformer of mol
    
    User can supply functional (e.g. AM1, HF, B3LYP, MP2) and basis set 
    (if needed, e.g. 6-31G, 6-311G(d,p) etc.). Options requests additional
    inputs (e.g. GD3BJ dispersion, polar). Method lets us do other calculation
    types (e.g. opt, freq)
    
    On the offchance the user only supplies a bunch of conformers we default to
    a single point calculation at the AM1 method. 
    """

    # delete existing directory of 'name' and contents, if it exists
    if os.path.exists(name):
        shutil.rmtree(name)
        print('Overwriting old ' + name + ' directory...')
        
    #make directory for .gjf files
    print('Making new ' + name + ' directory...')
    os.mkdir(name)
    
    # loop over all conformers
    for q in range(0,mol.GetNumConformers()):
        file_name = name + '_' + str(q) + '.gjf'         # Write out Gaussian input file
        with open(name + '/' + file_name, 'w') as f:
            f.write('%chk=' + name + '_' + str(q) + '.chk\n') # write the header
            f.write('%nprocshared=' + str(nproc) + '\n') # number of CPU cores
            f.write('%mem=' + str(vmem) + 'GB\n') # ammount of RAM
            f.write('#p ' + options + '\n\n') # route information
            f.write(name + '\n\n') # job name
            f.write('0 1\n')

            for i, atom in enumerate(mol.GetAtoms()): # add atomic coordinates
                positions = mol.GetConformer(q).GetAtomPosition(i)
                f.write(atom.GetSymbol() + '\t' + str(positions.x) + '\t' + str(positions.y) + '\t' + str(positions.z))
                f.write('\n')
            f.write('\n\n')
    #print('Wrote ' + str(q+1) + ' Gaussian input files for different conformations of ' + name) # pass output to notebook
    return file_name

def run_gjf(g_path,name='NoName'):
    """
    This function executes Gaussian G16 to perform calculations as written by write_gjf function.
    Copies the file in question to 'temp.gjf' in the Gaussian directory; executes, then copies back temp.out 
    to the correct place and names it accordingly.
    
    Loops over all *.gjf files.
    
    It would be useful to adjust this so it works for unix too:
    import os
    if os.name == 'nt':
        g_exec = 'g*.exe'
    
    if os.name == 'posix':
        g_exec = 'g*'???
        
    that last part is clearly flawed, but we can implement something.
    
    how to handle loading of gaussian module on HPC? do it in "if os.name == posix"?
    
    """
    g_folder = '/'.join(g_path.split('/')[0:-1])+'/'
    g_exec = g_path.split('/')[-1]
        
    path = r'*.gjf'  # *** PAY ATTENTION TO Gauss FORMAT! ***
    files = glob.glob(name+ '/' +path)
    
    ReturnDir = os.getcwd()
    print('Executing Gaussian jobs:')
    for i in tqdm(range(0, len(files))):
        shutil.copyfile(files[i], g_folder + 'temp.gjf')    # copy to gaussian folder
        os.chdir(g_folder)                                  # change to gaussian folder
        subprocess.run(g_exec + ' temp.gjf')              # launch gaussian (should be version agnostic)

        FilePath = os.path.join(ReturnDir, files[i][:-4] + '.out')
        shutil.copyfile('temp.out', FilePath)  # move those outputs
        os.chdir(ReturnDir) # go back to files
        
    #print('I ran ' + str(i + 1) + ' Gaussian jobs before returning to ' + ReturnDir)
    return()
    
def read_energy(GaussianOutputFile='NoName'):
    """
    This function reads Gaussian output and takes the single-point energy ONLY
    If you've done geometry optimisation this will be ignored!
    
    The code that looks for energy searches for a string ("SCF DONE: E(") - it is 
    not clear how robust this is, but it works for E(RAM1) calculations fine.
    
    Sometimes gaussian gives you .out files, sometimes .log. This is configured
    for .out files, because thats what my laptop is giving me today (17/1/23)
    """
    energy=[]
    g_output=open(GaussianOutputFile,'r')
    for line in g_output: # find the line that contains the SCF energy and store
        if 'SCF Done:' in line: # *** possible weakpoint in code ***
            Scraped = str([line.strip('\n')])
            Scraped = Scraped.split('A.U.',1)[0]
            Scraped = Scraped.split(' = ',1)[1]
            energy.append((float(Scraped) * 627.5)) # extract the energy, convert to Kcal mol-1 (1Ha = 627.5 Kcal mol-1)
    g_output.close()
    
    return(energy)    
    
def read_g_log(name='NoName'):
    """
    This function reads Gaussian output and takes the single-point energy ONLY
    If you've done geometry optimisation this will be ignored!
    
    The code that looks for energy searches for a string ("SCF DONE: E(") - it is 
    not clear how robust this is, but it works for E(RAM1) calculations fine.
    
    Sometimes gaussian gives you .out files, sometimes .log. This is configured
    for .out files, because thats what my laptop is giving me today (17/1/23)
    """

    import glob 
    
    # make a list of *.out files we want to read  and their path 
    path = r'*.out' # *** PAY ATTENTION TO Gauss FORMAT! ***
    files = glob.glob(name+'/'+path)
    energy = []
    
    for i in range(0,len(files)):
        GaussianOutputFile = files[i] # loop over all *.out files and read into memory
        energy.append(read_energy(GaussianOutputFile)[-1]) # return final energy, incase its a geometry opt job

    if len(energy) == 1:
        print('No Gaussian output file(s) named "' + name + '" were found')
        return # don't return anything if we don't have an input
        
    print(len(energy))    
    return(energy)
