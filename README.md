Repo for LC conformational analysis tool.
<br><br>
Current version - V 1.0:<br> 
Future Plans:
<br> Constrints, ability to fix certain torsions
<br> Ability to use other computational engines (MOPAC, ORCA, TINKER, OpenMM)
<br> 
<br><br>

## Background
Many liquid crystalline materials can adopt multiple conformations, however for various reasons this is often neglected in computational treatment which typically focuses on a single conformation. This software tool is designed as a simple means to rapidly generate multiple conformers of semi-flexible molecules (such as liquid crystals) using rules-based methods embedded in rdkit. While primarily intended as a grapical means to explore the overall shape of liquid crystalline dimers and oligomers, the conformer tool is a easy way to access the powerful Rdkit conformer generation engine and can be applied to any system of interest.

## Installation 
* This tool is written in Python, so if you are on Windows you'll need to have a working Python install. You could use Anaconda for simplicity.
* With python instaled, navigate to the conformer directory. 
* Install dependent packages with conda (rdkit, numpy, matplotlib, tkinter, tqdm; you can do this with <br><br>```pip install -r requirements.txt```<br><br>

## Usage
* Launch the gui:<br><br>``` python gui.py ```<br><br>
  ![image](https://github.com/RichardMandle/conformer/assets/101199234/3dc3c54c-4ebf-4749-941c-f4d8f712f35a)
  <br>
* Load your molecule as a .mol, .mol2 or .sdf file (can be saved from ChemDraw) or even type in the smiles/smarts string 
* You can save your session here (including generated conformers) or reload a previous session 
* If you want to calculate a bend angle distribution you need to select two pairs of atoms that define the bend angle. You can do this by clicking each atom in the pair:
  ![image](https://github.com/RichardMandle/conformer/assets/101199234/0db54cfb-d9e2-4695-90eb-617f79ef9c5e)
  <br>
* the orientation of the vectors WILL affect the angle, currently. You want to select four atoms, these being two pairs of atoms that define a vector, making sure you click the "furthest" atom first for each pair (e.g. 'N' of CN). These get highlighted by different colours. Clicked the wrong one? If you click a 5th atom or the 'clear' button and it'll reset:
  ![image](https://github.com/RichardMandle/conformer/assets/101199234/a2e8454d-e471-4f14-9a9c-d385aa14f983)
  <br>
* In the "conformer search" tab set your desired options. Number of conformers to screen, the desired method used to generate conformers, the title of your job, and the evaluation method. You can specify the groups that define the vectors as smiles/smarts strings here if you prefer. We evaluate the energy of each conformer; you can use the inbuilt "MMFF" method which is fast. The advanced settings button lets you customise the conformer search options even more:
  ![image](https://github.com/RichardMandle/conformer/assets/101199234/19f3f7e6-7384-4cc1-b5dc-729b5ac35fec)
  <br>
*  Or you can run an external calculation through Gaussian (you need to set the path to the Gaussian executable in 'Settings'). This has a few extra options, method, number of CPU cores, ammount of RAM. If you select an "opt" job then the conformers will be optimised with your chosen method (after the RMS pruning etc.)
* Once the "Start Conformer Search" button is pressed you the conformer search begins. The window text will change from "Conformer generation in progress" to "Conformer generation complete" when finished; any errors will appear in the terminal.
  <br>
* Once complete we can navigate to the "Analysis #1" tab to visualise the first of our results as a histogram of probability of a given bend angle. A slider controls the bin width of the histogram, and you can specify the temperature (in K) used in the probability calculation. You can fit a Gaussian to the data, save the data (or histogram data) as .csv, or even just straight to .png. Some basic stats are also printed (mean angle, median angle etc.):
  ![image](https://github.com/RichardMandle/conformer/assets/101199234/32e13df2-90cd-438f-867a-1dc1d35b7447)
  <br>
* You can change some options on the plot and save the data (or histogram data) to .csv for proper plotting using other tools. You can save the plot shown in the Conformer gui tool if you like with "Save Plot as PNG".
* What at it you clicked the wrong atoms when defining the bend angle? No problem! Just navigate back to the first tab (Load/Save Molecule), reselect your atoms. Go to tab 2, click "recalculate bend angle". Done, and quickly too.
  <br>
* On the "Analysis #2" tab we can plot various 3d descriptors from rdkit against bend angle (useful) or against each other (less useful I think). Basic options control ploting, saving to .csv (I'll add .png later), and optional fitting with a 1,2,3,4 order polynomial:
  ![image](https://github.com/RichardMandle/conformer/assets/101199234/d5cd44c0-3020-4778-857b-3e7275a25ad3)

* Again, you can export data as .csv or save the plot as .png.
  <br>
* In the "View Conformers" tab we can do a quick inspection of the 3D geometry of our various conformers. A slider lets us select different conformers, and we can save the information to .sdf file for later retrieval.
* There is an interface with PyMol (which is very slow). It is _much faster_ to save the conformers as a .sdf file and reload into PyMol that way. You can use any visualisation software of course, but PyMol, VMD and QuteMol make especially nice images.
* If you export directly to PyMol, which is not really recommended (slow), then the check boxes enable raytracing, grid_image and default image style of the University of Leeds Soft Matter Physics group (UoL-SMP) change the appearance of the resulting image.
  ![image](https://github.com/RichardMandle/conformer/assets/101199234/836b0426-48cc-4995-833b-6bb46cff8a71)

* The PyMol interface can make nice images; of course, this is an interactive process, but to give an idea of the sort of output the 'default' options (raytrace=True,grid_image=True,UoL-SMP=True) give see below for the LC trimer CBO6B6OCB:
![image](https://github.com/RichardMandle/conformer/assets/101199234/1dd3e07a-7881-4e12-87d9-d9a6744baded)


## Citation
This work is described in Mandle _et al_ "Rapid Conformational Analysis of Semi-Flexible Liquid Crystals", however the following papers should also be cited if you make use of this softare tool:

* Riniker, S.; Landrum, G.A.J.J.o.c.i.; modeling. Better informed distance geometry: Using what we know to improve conformation generation. 2015, 55, 2562-2574.
* Vainio, M.J.; Johnson, M.S.J.J.o.c.i.; modeling. Generating conformer ensembles using a multiobjective genetic algorithm. 2007, 47, 2462-2474.
* Wang, S.; Witek, J.; Landrum, G.A.; Riniker, S. Improving conformer generation for small rings and macrocycles based on distance geometry and experimental torsional-angle preferences. J Chem Inf Model 2020, 60, 2044-2058.
* Landrum, G. Rdkit: Open-source cheminformatics. 2006.


Any comments/questions/requests email r<dot>mandle<at>leeds<dot>ac<dot>uk
