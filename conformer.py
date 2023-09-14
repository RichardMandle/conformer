import sys
import os
import shutil
import glob
import time
import subprocess
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem, Draw

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
    
def get_default_search_settings():
    default_settings = {
        'rms_threshold': 0.33,
        'use_torsion_pref': True,
        'use_knowledge': True,
        'use_random_coords': False,
        'MMFF Variant': "MMFF94",
        'random_seed': 61453,
        'opt': False,
        'max_opt_iter': 5000,
        'min_energy_MMFF': 500,
        'gaussian_job': False,
        'gaussian_options': "",
        'gaussian_nproc': "",
        'gaussian_vmem': ""
    }
    return default_settings
def get_conformer_methods(Print=False):
    
    embeded_methods = ['ETKDGv2',
                       'ETDG',
                       'ETKDG',
                       'ETKDGv3',
                       'srETKDGv3',
                      ]
    
    if print == True:
        print('Available conformer search methods:')
        for method in embeded_methods:
            print(method)
            
    return(embeded_methods)
    
def get_3d_descriptors(mol,angle = None):
    '''
    use d3d to get 3d descriptors for mol
    
    args:
    mol - a mol object with 3d conformer information
    
    returns:
    lets see
    '''
    
    # psuedo code
    
    property_dict = {}
    
    for i in range(0,mol.GetNumConformers()):
        property_dict[i] = {
        'asphericity': d3d.Asphericity(mol,confId=i),
        'eccentricity': d3d.Eccentricity(mol,confId=i),
        'ISF': d3d.InertialShapeFactor(mol,confId=i),
        'NPR1': d3d.NPR1(mol,confId=i),
        'NPR2': d3d.NPR2(mol,confId=i),
        'PMI1': d3d.PMI1(mol,confId=i),
        'PMI2': d3d.PMI2(mol,confId=i),
        'PMI3': d3d.PMI3(mol,confId=i),
        'ROG': d3d.RadiusOfGyration(mol,confId=i),
        'spherocity': d3d.SpherocityIndex(mol,confId=i),
        'bend angle': angle[i] if angle != None else None
        }
        
    
    return property_dict

def conf_gen(mol,
             num_of_conformer=500,
             embeded_method='ETKDGv2', 
             name='NoName',
             rms_thresh = 0.33,
             use_torsion_pref = True,
             use_knowledge = True,
             use_random_coords = False,
             random_seed = 61453,
             opt = False,
             max_iter=5000,
             min_energy_MMFF=500,
             MMFF_variant="MMFF94"):

    '''
    This function makes the call to rdkit to do a conformer generation
    Args:
    mol - the mol object to generate conformers for
    num_of_conformers - the number of conformers we want to generate
    name - the name of the job (used only when then writing output to terminal)
    rms_thresh - the threshold of RMS used to prune identical conformers.
    use_torsion_pref - use knowledge of preferred torsions in conformer generate
    use_knowledge - use basic chemical knowledge (e.g. flat benzene)
    use_random_coords - generates conformers by randomising the coordinates?
    randokm_seed - ability to set a different random seed (give integer)
    opt - calls a geometry optimizer using MMFF in rdkit (default = False)
    max_iter - if opt=True, the maximum iterations we'll use in our optimisation 
    min_energy_MMFF - the minimum energy to accept for a conformer if opt = True (in kcal mol?)
    
    Returns:
    mol - the mol object with conformers embedded
    energy - the energy of conformers within the mol object (above)
    '''
    
    embeded_methods = get_conformer_methods(Print=False)
    if embeded_method not in embeded_methods:
        print('I\'ve not heard of the ' +  embeded_method + ' method, so I\'ll use the ' + embeded_methods[0] + ' instead.')
        embeded_method = embeded_methods[0] # set to first one.
    
    print('\nusing the ' + embeded_method + ' method')
    
    ps = getattr(AllChem, embeded_method)()
    
    ps.pruneRmsThresh = rms_thresh
    ps.useExpTorsionAnglePrefs = use_torsion_pref
    ps.useBasicKnowledge = use_knowledge
    ps.useRandomCoords=use_random_coords
    ps.randomSeed = random_seed

    AllChem.EmbedMultipleConfs(mol, num_of_conformer, ps)
    
    if opt == True:
        AllChem.MMFFOptimizeMoleculeConfs(mol,maxIters=max_iter) # we can, but don't need to, run MMFF optimizer
    
    if num_of_conformer != mol.GetNumConformers():
        print('After pruning with an RMSD threshold of ' + str(ps.pruneRmsThresh) 
              + '\n' + str(num_of_conformer) + ' conformers were pruned to ' 
              + str(mol.GetNumConformers()) + ' conformer' + ('s'*(mol.GetNumConformers()>1))) # this bit at the end is to get the plural right...

    energy = [] # Conformer Energy in kcal mol-1
    mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=MMFF_variant)
    for i in range(0,mol.GetNumConformers()):
        ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=i)
        ff.Minimize()
        energy.append(ff.CalcEnergy())

    with open(name + '_conformer_energies_.txt', 'w') as file:
        file.write('\n'.join(str(energ) for energ in energy))
    
    return(mol,energy)
        
def conf_analysis(mol_conf, 
                energy,
                vec_def_method='atoms',
                atom_idx_list=[],
                vector1 = 'N#C',
                vector2 = 'N#C',
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
    
    mol_conf = input molecule conformers as rdkit mol object
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
    
    """

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
            
        for i in range(0,mol_conf.GetNumConformers()):
            angle.append(find_dimer_angle(mol_conf, atom_idx_list, confId=i))
            
    if vec_def_method == 'smi' or vec_def_method == 'sma':
        Pattern = []
        if (vector1) == (vector2): # if they are the same, skip.
            if vec_def_method == 'sma':
                Pattern = mol_conf.GetSubstructMatches(Chem.MolFromSmarts(vector1))
            if vec_def_method == 'smi':
                Pattern = mol_conf.GetSubstructMatches(Chem.MolFromSmiles(vector1))     
                
        if (vector1) != (vector2):
            for vector in [vector1,vector2]:
                if vec_def_method == 'sma':
                    Pattern += (mol_conf.GetSubstructMatches(Chem.MolFromSmarts(vector)))
                if vec_def_method == 'smi':
                    Pattern += (mol_conf.GetSubstructMatches(Chem.MolFromSmiles(vector)))  
        print('matched pattern to atoms: ' + str(Pattern))        
        print('\nCalculating bend-angle based on a vector between ' +
            ((vec_def_method=='smi') * 'smi') + ((vec_def_method=='sma') * 'sma') +
            ' matches for ' + vector1 + ' and ' + vector2)
            
        for i in range(0,mol_conf.GetNumConformers()):
            angle.append(find_dimer_angle_by_group(mol_conf, Pattern, confId=i))
    
    with open(str(name) + '_conformer_angles_.txt', 'w') as file:
        file.write('\n'.join(str(angl) for angl in angle))
    
    if Gauss: # here, if we call Gauss=True (or Gauss = anything) we send the 
        gau.write_gjf(mol_conf,name,options,nproc=cores,vmem=ram) # write the Gaussian input files
        if write_only == False:
            if not g_path:
                print('\nUh oh!\nPath to Gaussian executable not set!')
            gau.run_gjf(g_path,name) # Execute Gaussian
        energy = gau.read_g_log(name) # Read in Gaussian energies
        
    return(energy,angle)

def calc_bend_hist(angle,energy,fit_gaussian=True,show_plot=False,name='NoName',Temp=298,BinSteps=10):
    """
    This function calculates the Boltzmann probability of each conformer
    and creates a histogram of probability vs bend angle
    
    Some statistical data on the degree of bend is also returned
    
    args:
    angle - array of angles for each conformer
    energy - array of energies of each conformer (in kcal / mol-1)
    fit_gaussian (bool) - fit a gaussian to the histogram, or not?
    show_plot - do we want to show the plot (true), or have it returned as a figure object (false)
    name - the name of the material, used in saving
    Temp - the temperature at which to calculate the probability
    BinSteps - the size of the bins used in the histogram
    
    returns:
    """

    Kb = 0.008314 # Boltzmann's constant
    deltaE_KJ = (energy - np.min(energy))* 4.184
    probability = np.exp(-deltaE_KJ/(Kb*Temp)) / np.sum(np.exp(-deltaE_KJ/(Kb*Temp)))

    y, x = np.histogram(angle, weights = probability,bins=np.arange(0,190,BinSteps))
    
    px = 1/plt.rcParams['figure.dpi']
    fig = plt.figure(figsize=(600*px, 400*px))
    plt.bar(x[:-1], y,width=BinSteps*0.9)
    plt.xlabel('Bend Angle / Degrees')
    plt.ylabel('Probability at ' + str(Temp) + ' K')
    plt.xlim([0,180])
    plt.ylim([0,1])
    plt.title(name)
    
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
        
    if fit_gaussian == False:
        x_max = np.nan
        fwhm = np.nan
        
    plt.savefig(name + 'prob_vs_angle.png', bbox_inches='tight')
    
    with open(name + '_conformer_probability_.txt', 'w') as file:
        file.write('\n'.join(str(probabilit) for probabilit in probability))
    
    if show_plot == True:
        plt.show()
    else:
        return (fig,deltaE_KJ,probability,x,y,x_max,fwhm)
        
    return(deltaE_KJ,probability,x,y)

def analyse_angles(hist,bin_edges):
    """
    This function just fits a gaussian on top of a histogram and returns the x/y values for plotting
    and also the X-at-Y-max (x_max) and FWHM.
    
    Args:
    hist      - histogrammed angle data
    bin_edges - the bin edges for the histogram
    figure    - the figure handle to plot to
    
    Returns:
    x           - linspace data used for plotting the Gaussian fit
    y           - The Y- of the gaussian function w.r.t. x&y
    x_max       - the X value at maximum Y of the Gaussian function
    fwhm        - the FWHM of the gaussian function.
    """
    
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    total_counts = np.sum(hist)
    mean = np.sum(bin_centers * hist) / total_counts

    variance = np.sum(((bin_centers - mean) ** 2) * hist) / total_counts  # Compute the weighted variance

    a_fit, b_fit, c_fit = [np.max(hist), mean, np.sqrt(variance)] # Fit the Gaussian function using the estimated parameters
    
    x_max = np.round(b_fit,2) 
    fwhm = np.round((2 * np.sqrt(2 * np.log(2)) * c_fit),2)   # Calculate the x-maximum and FWHM of the Gaussian

    x = np.linspace(np.min(bin_centers), np.max(bin_centers), 1000)     # Generate x-values for the Gaussian curve

    y = maths.gaussian(x, a_fit, b_fit, c_fit)     # Compute the y-values for the Gaussian curve

    # Return the x&y values of the fit along with x-maximum and FWHM
    return x,y, x_max, fwhm