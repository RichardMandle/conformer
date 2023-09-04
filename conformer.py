import sys
import os
import shutil
import glob
import time
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, PyMol

import numpy as np

import matplotlib.pyplot as plt

from rdkit.Chem.Lipinski import RotatableBondSmarts
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions 
from rdkit.Chem import Descriptors3D as d3d

# import our own files here
import maths
import gauss as gau

def find_dimer_angle_by_group(mol,Pattern, confId=0):
    """
    This function returns the angle between two vectors which are defined by the
    x,y,z coordinates of the atom IDs in the molecule 'mol' given in the variable 'Pattern'
    
    Pattern is itself determined by a substring match to the SMILES/SMARTS code of 'mol'
    """

    Vec1 =  mol.GetConformer(confId).GetAtomPosition((Pattern[0])[0]) - mol.GetConformer(confId).GetAtomPosition((Pattern[0])[1])
    Vec2 =  mol.GetConformer(confId).GetAtomPosition((Pattern[1])[0]) - mol.GetConformer(confId).GetAtomPosition((Pattern[1])[1])
    Dimer_Angle = np.degrees(maths.angle_between(Vec1, Vec2))
    
    return(Dimer_Angle)
    
def find_dimer_angle(mol,atom_idx_list, confId=0):
    """
    This function returns the angle between two vectors which are defined by the
    x,y,z coordinates of the 4x atom ids in atom_idx_list.
    
    """

    Vec1 =  mol.GetConformer(confId).GetAtomPosition(atom_idx_list[0]) - mol.GetConformer(confId).GetAtomPosition(atom_idx_list[1])
    Vec2 =  mol.GetConformer(confId).GetAtomPosition(atom_idx_list[2]) - mol.GetConformer(confId).GetAtomPosition(atom_idx_list[3])
    Dimer_Angle = np.degrees(maths.angle_between(Vec1, Vec2))
    
    return(Dimer_Angle)

def get_conformer_methods(Print=False):
    
    embeded_methods = ['ETKDGv2',
                       'ETDG',
                       'ETKDG',
                       'ETKDGv3',
                       'srETKDGv3',
                      ]
    
    if print == False:
        print('Available conformer search methods:')
        for method in embeded_methods:
            print(method)
            
    return(embeded_methods)

def conf_gen(mol,
             num_of_conformer=500,
             embeded_method='ETKDGv2', 
             name='NoName',
             max_iter=5000,
             min_energy_MMFF=500,
             opt = False):
    
    # a list of embeded methods; if the embeded_method isn't here, we'll
    # simply go for a standard one.

    embeded_methods = get_conformer_methods(Print=False)
    
    if embeded_method not in embeded_methods:
        print('I\'ve not heard of the ' +  embeded_method + ' method, so I\'ll use the ' + embeded_methods[0] + ' instead.')
        embeded_method = embeded_methods[0] # set to first one.
    
    print('\nusing the ' + embeded_method + ' method')
    
    ps = getattr(AllChem, embeded_method)()
    ps.pruneRmsThresh = 0.33

    ps.useExpTorsionAnglePrefs = True
    ps.useBasicKnowledge = True
    ps.useRandomCoords=False
    ps.randomSeed = 0xf00d
    
    AllChem.EmbedMultipleConfs(mol, num_of_conformer, ps)
    
    if opt == True:
        AllChem.MMFFOptimizeMoleculeConfs(mol,maxIters=max_iter) # we can, but don't need to, run MMFF optimizer
    
    if num_of_conformer != mol.GetNumConformers():
        print('After pruning with an RMSD threshold of ' + str(ps.pruneRmsThresh) 
              + '\n' + str(num_of_conformer) + ' conformers were pruned to ' 
              + str(mol.GetNumConformers()) + ' conformer' + ('s'*(mol.GetNumConformers()>1))) # this bit at the end is to get the plural right...
    
    # write conformations
    writer = Chem.SDWriter(name + '_conformers.sdf')
    for cid in range(mol.GetNumConformers()):
        writer.write(mol, confId=cid)

    energy = [] # Conformer Energy in kcal mol-1
    mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94')
    for i in range(0,mol.GetNumConformers()):
        ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=i)
        ff.Minimize()
        energy.append(ff.CalcEnergy())

    with open(name + '_conformer_energies_.txt', 'w') as file:
        file.write('\n'.join(str(energ) for energ in energy))
    
    return(mol,energy)
    
def rad_gyr_gen(mol, energy, angle, name):
    '''
    Generates a list of all the radii of gyration for the different 
    conformers in 'mol'
    '''

    rmslist = []
    AllChem.AlignMolConformers(mol, RMSlist=rmslist)
    
    radii = []
    for i in range(0,mol.GetNumConformers()):
        radii.append(d3d.RadiusOfGyration(mol, confId = i, useAtomicMasses = True))
    radii = np.array(radii, dtype = float)
    print(radii.mean(), radii.std(), len(radii))
    
    fig0, ax0 = plt.subplots()
    plt.scatter(angle[1:],radii)

    plt.xlabel('Bend Angle / Degrees')
    plt.ylabel('Radius of Gyration ') 
    
    plt.savefig(name + 'Gr_vs_angle.png', bbox_inches='tight')
    
    return radii
    
def find_bond_groups(smi):
    """
    Find groups of contiguous rotatable bonds and return them sorted by decreasing size
    Taken from RDKit cookbook (https://rdkit.org/docs/Cookbook.html) on 23/01/2023
    """
    mol = Chem.MolFromSmiles(smi)
    rot_atom_pairs = mol.GetSubstructMatches(RotatableBondSmarts)
    rot_bond_set = set([mol.GetBondBetweenAtoms(*ap).GetIdx() for ap in rot_atom_pairs])
    rot_bond_groups = []
    while (rot_bond_set):
        i = rot_bond_set.pop()
        connected_bond_set = set([i])
        stack = [i]
        while (stack):
            i = stack.pop()
            b = mol.GetBondWithIdx(i)
            bonds = []
            for a in (b.GetBeginAtom(), b.GetEndAtom()):
                bonds.extend([b.GetIdx() for b in a.GetBonds() if (
                    (b.GetIdx() in rot_bond_set) and (not (b.GetIdx() in connected_bond_set)))])
            connected_bond_set.update(bonds)
            stack.extend(bonds)
        rot_bond_set.difference_update(connected_bond_set)
        rot_bond_groups.append(tuple(connected_bond_set))
    return tuple(sorted(rot_bond_groups, reverse = True, key = lambda x: len(x)))


def count_rot_bonds(smi):
    """
    retrun the largest number of contiguous rotatable bonds
    
    """
    
    bond_groups = find_bond_groups(smi)
    largest_n_cont_rot_bonds = len(bond_groups[0]) if bond_groups else 0

    #print('Number of contiguous rotatable bonds is: ' + str(largest_n_cont_rot_bonds))
    #print('Suggested number of conformers to screen: ' + str(3**largest_n_cont_rot_bonds))
    return (largest_n_cont_rot_bonds)   

def est_num_confs(smi,library_size='standard'):
    """
    Simple function to calculate the number of conformers we might want to study
    A few are predefined, or the user can pass a custom 'multiplier'
    the number of conformers is simply:

        Number of contiguous rotatable bonds ** Multiplier
    
    Probably worth changing from the LARGEST number of contiguous rotatable bonds to 
    the product of the square of the number of contiguous rotatable bonds?
    """
    
    largest_n_cont_rot_bonds= count_rot_bonds(smi)
    
    lib_names = ['tiny','small','standard','large','huge']
    lib_sizes = [1,2,3,4,5]

    if type(library_size) == int:
        ConfMultiplier = int(library_size)
        print('User passed an integer; generating ' + str(ConfMultiplier) + ' conformers')
        n_confs = ConfMultiplier
        
    if type(library_size) == str:
        if library_size in lib_names:
            ConfMultiplier = (lib_sizes[lib_names.index(library_size)])
            print('User requested the ' + library_size + ' multiplier (' + str(ConfMultiplier) + str(')'))
        if library_size not in lib_names:
            library_size = 'standard'
            ConfMultiplier = (lib_sizes[lib_names.index('standard')])
            print('Unexpected input; defaulting to standard multiplier (' + str(ConfMultiplier) + str(')'))

        n_confs = ConfMultiplier**largest_n_cont_rot_bonds
  
    print('The number of conformers to screen is ' + str(n_confs))
    
    return(n_confs)
    
def conf_search(mol,
                vec_def_method='atoms',
                atom_idx_list=[],
                embeded_method = 'ETKDGv2',
                vector1 = 'N#C',
                vector2 = 'N#C',
                n_conf=50,
                name='NoName',
                g_path = None,
                Gauss=False,
                options='',
                cores=1,
                ram=1,
                write_only=False):
    """
    Generates n_conf number of 3D conformers of the smiles string 'smi' and calcualtes the angle between
    nitrile groups.
    
    Conformer cordinates are written to .sdf by conf_gen
    Energies are written to .txt by another function (check)
    Intervector angles are written to .txt by THIS FUNCTION
    
    mol = input molecule as rdkit mol object
    g_path = path to Gaussian executable
    n_conf = number of conformers to generate
    name   = the name of the system, used for saving
    Gauss  = if is True, send confs to Gaussian for energy/optimisation (default = False)?
    
    if Gauss = True
        functional = calculation type to use in Gaussian (e.g. AM1, PM3, HF, B3LYP, MP2 etc.)
        method     = Gaussian job type (single point energy is default, options like opt, freq, polar etc.)
        basis      = Basis set for Gaussian calculations (e.g. 6-31G etc.)
        options    = Any other Gaussian options (e.g. GD3BJ, integral=ultrafine etc.)
    
    write_only lets us leave out the actual calculation and just prep the gaussian input files.
    Code is inflexible for unsymmetrical systems, but would be an easy fix.
    """
    
    mol,energy = conf_gen(mol,n_conf,embeded_method,name)
    
    angle = [] # Conformer Angle in Degrees
    if vec_def_method == 'atoms': 
        if len(atom_idx_list) != 4:
            print('\nUh oh!\nYour list of atoms (atom_idx_list) is too' + 
            ' short'*(len(atom_idx_list)<4) + 
            ' long'*(len(atom_idx_list)>4) + 
            ' with a total of ' + str(len(atom_idx_list)) + ' entries!')
            return
            
        print('\nCalculating bend-angle based on a vector between atoms ' +
            str(atom_idx_list[0:2]) + ' and ' + str(atom_idx_list[2:4]) + '\n')
            
        for i in range(0,mol.GetNumConformers()):
            angle.append(find_dimer_angle(mol, atom_idx_list, confId=i))
            
    if vec_def_method == 'smi' or vec_def_method == 'sma':
        Pattern = []
        if (vector1) == (vector2): # if they are the same, skip.
            if vec_def_method == 'sma':
                Pattern = mol.GetSubstructMatches(Chem.MolFromSmarts(vector1))
            if vec_def_method == 'smi':
                Pattern = mol.GetSubstructMatches(Chem.MolFromSmiles(vector1))     
                
        if (vector1) != (vector2):
            for vector in [vector1,vector2]:
                if vec_def_method == 'sma':
                    Pattern += (mol.GetSubstructMatches(Chem.MolFromSmarts(vector)))
                if vec_def_method == 'smi':
                    Pattern += (mol.GetSubstructMatches(Chem.MolFromSmiles(vector)))  
        print('matchesd pattern to atoms: ' + str(Pattern))        
        print('\nCalculating bend-angle based on a vector between ' +
            ((vec_def_method=='smi') * 'smi') + ((vec_def_method=='sma') * 'sma') +
            ' matches for ' + vector1 + ' and ' + vector2)
            
        for i in range(0,mol.GetNumConformers()):
            angle.append(find_dimer_angle_by_group(mol, Pattern, confId=i))
    
    with open(str(name) + '_conformer_angles_.txt', 'w') as file:
        file.write('\n'.join(str(angl) for angl in angle))
    
    if Gauss: # here, if we call Gauss=True (or Gauss = anything) we send the 
        gau.write_gjf(mol,name,options,nproc=cores,vmem=ram) # write the Gaussian input files
        if write_only == False:
            if not g_path:
                print('\nUh oh!\nPath to Gaussian executable not set!')
            gau.run_gjf(g_path,name) # Execute Gaussian
        energy = gau.read_g_log(name) # Read in Gaussian energies
        
    return(mol,energy,angle,name)

def calc_bend_hist(angle,energy,fit_gaussian=True,show_plot=False,name='NoName',Temp=298,BinSteps=10):
    """
    This function calculates the Boltzmann probability of each conformer
    and creates a histogram of probability vs bend angle
    
    Some statistical data on the degree of bend is also returned
    """

    Kb = 0.008314 # Boltzmann's constant

    deltaE_KJ = (energy[1:] - np.min(energy[1:]))* 4.184

    probability = np.exp(-deltaE_KJ/(Kb*Temp)) / np.sum(np.exp(-deltaE_KJ/(Kb*Temp)))

    y, x = np.histogram(angle[1:], weights = probability,bins=np.arange(0,190,BinSteps))
    fig = plt.figure()
    plt.bar(x[:-1], y,width=BinSteps*0.9)

    plt.xlabel('Bend Angle / Degrees')
    plt.ylabel('Probability at ' + str(Temp) + ' K')
    
    plt.xlim([0,180]) # set the limits on x and y
    plt.ylim([0,1])
    
    plt.title(name)
    
    print('Mean bend-angle is: ' + str(np.round(np.mean(angle[1:]),1)) + ' degrees, with standard deviation of: ' + str(np.round(np.std(angle[1:]),1)) + ' degrees')
    print('Median bend-angle is: ' + str(np.round(np.median(angle[1:]),1)) + ' degrees')
    print('Major distribution centred at: ' + str(x[np.argmax(y)]))
    
    if fit_gaussian == True:
        gaussian_x,gaussian_y,x_max,fwhm = analyse_angles(y, x)
        plt.plot(gaussian_x,gaussian_y,color='r')

        # Plot the dashed line along the FWHM
        plt.axvline(x_max, 
                    linestyle='--', 
                    color='g', 
                    alpha=0.66, 
                    label='Gaussian Max.')

        # Calculate the half-max and 1.5x max y values
        half_max = maths.gaussian(x_max, np.max(y), x_max, np.sqrt(fwhm / (2 * np.sqrt(2 * np.log(2)))))/2

        # Plot the vertical line from half-max to 1.5x max y
        plt.plot([x_max-fwhm/2,x_max+fwhm/2],[half_max,half_max],  
                 linestyle='--', 
                 color='g', 
                 alpha=0.66, 
                 label='FWHM')
        

    plt.savefig(name + 'prob_vs_angle.png', bbox_inches='tight')
    
    with open(name + '_conformer_probability_.txt', 'w') as file:
        file.write('\n'.join(str(probabilit) for probabilit in probability))
    
    if show_plot == True:
        plt.show()
    else:
        return fig
        
    return(deltaE_KJ,probability,x,y)

def analyse_angles(hist,bin_edges):
    """
    Args:
    hist      - histogrammed angle data
    bin_edges - the bin edges for the histogram
    figure    - the figure handle to plot to
    """
    
    # Compute the histogram bin centres
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Compute the total counts and weighted mean
    total_counts = np.sum(hist)
    mean = np.sum(bin_centers * hist) / total_counts

    variance = np.sum(((bin_centers - mean) ** 2) * hist) / total_counts  # Compute the weighted variance

    a_fit, b_fit, c_fit = [np.max(hist), mean, np.sqrt(variance)] # Fit the Gaussian function using the estimated parameters
    
    x_max = np.round(b_fit,2) 
    fwhm = np.round((2 * np.sqrt(2 * np.log(2)) * c_fit),2)   # Calculate the x-maximum and FWHM of the Gaussian

    x = np.linspace(np.min(bin_edges), np.max(bin_edges), 1000)     # Generate x-values for the Gaussian curve

    y = maths.gaussian(x, a_fit, b_fit, c_fit)     # Compute the y-values for the Gaussian curve

    print('A Gaussian fit to this data has a maximum at: ' + str(x_max) + ' with a FWHM of ' + str(fwhm))
    
    # Return the x&y values of the fit along with x-maximum and FWHM
    return x,y, x_max, fwhm

 
def pymol_grid_img(mol,style='yes',ray='no',start=0,end=-1):

    """
    Function for calling pymol and making a nice image
    having style = 'yes' will give a standard image in the UoL-SMP style
    raytacing defaults to off (ray='no') but can be called ('yes') or run externally in PyMol
    
    Make sure you update the pymol path in the os.system call
    
    you can specify start = i and end = m, where i and m are integers, to
    output only conformers between the ith and mth to Pymol. Python
    indexes from zero (i.e. hte first conformer is no 0, the 2nd is 1 etc.)
    """
    
    if end == -1:
        end = mol.GetNumConformers()
        
    os.system('cmd /k "\"C:\Program Files\Pymol\PyMOLWin.exe\" -R"')
    v= PyMol.MolViewer()
    v.DeleteAll()
    
    for i in range(start,end):
        v.ShowMol(mol, confId=i,name=name+'_conf-%d'%i,showOnly=False)
    
    if style.lower()=='yes' or style.lower() == 'y':
        # call grid mode; coloring options
        # putting multiple v.server.do calls on one line speeds it up a lot.
        v.server.do('set grid_mode, on;'+
                    'set_color oxygen, [1.0,0.4,0.4];'+
                    'set_color nitrogen, [0.5,0.5,1.0];'+
                    'set_color carbon, [0.5,0.5,0.5];'+
                    'set_color hydrogen, [1.0,1.0,1.0];'+
                    'util.cbag;recolor;'+
                    'remove solvent;'+
                    'as spheres;'+
                    'bg white;'+
                    'set light_count,8;'+
                    'set spec_count,1;'+
                    'set shininess, 10;'+
                    'set specular, 0.25;'+
                    'set ambient,0;'+
                    'set direct,0;'+
                    'set reflect,1.5;'+
                    'set ray_shadow_decay_factor, 0.1;'+
                    'set ray_shadow_decay_range, 2;'+
                    'unset depth_cue')

        if (ray.lower() == 'y' or ray.lower() == 'yes'):
            print('I am going to raytrace')
            v.server.do('ray') # you can raytrace if you want

    return() 