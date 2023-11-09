![image](https://github.com/RichardMandle/conformer/assets/101199234/f3f93b59-ba1c-43b6-9bca-8e42b25973c3)# conformer
Repo for LC conformational analysis tool.
<br><br>
Current version - V 0.2:<br> backend code working <br> usable GUI <br>
Future Plans:
<br> implement reading external .log files
<br> progress bar might be useful for large/slow jobs
<br><br>
## Usage
* Install dependent packages with conda (rdkit, numpy, matplotlib, tkinter, tqdm; you can do this with <br><br>```pip install -r requirements.txt```<br><br>
* Launch the gui:<br><br>``` python gui.py ```<br><br>
![image](https://github.com/RichardMandle/conformer/assets/101199234/60019c39-1fda-4d30-9aca-585cdaa89a9e)
* Load your molecule as a .mol file (can be saved from ChemDraw) or type in the smiles/smarts string <br>
* You can save your session here (including generated conformers) or reload a previous session <br>
* Time to select the atoms that define the two vectors we'll use to calcualte the angle. You can enter the SMILES or SMARTS strings of representative groups:<br>
* ![image](https://github.com/RichardMandle/conformer/assets/101199234/78074ac2-32fc-4a0d-97dc-b3b650078dc8)

* ... Or you can click on atoms on the first tab to define the vectors (the orientation of the vectors WILL affect the angle, currently. You want to select four atoms, these being two pairs of atoms that define a vector, making sure you click the "furthest" atom first for each pair (e.g. 'N' of CN). These get highlighted by different colours. Clicked the wrong one? If you click a 5th atom it'll just reset:
* ![image](https://github.com/RichardMandle/conformer/assets/101199234/ec2a0521-ae0f-42fa-812c-7c4d49b20e80)
* In the "conformer search" tab set your desired options. Number of conformers to screen, the desired method used to generate conformers, the title of your job, and the evaluation method. We evaluate the energy of each conformer; you can use the inbuilt "MMFF" method which is very fast. The advanced settings button lets you customise the conformer search options even more:
* ![image](https://github.com/RichardMandle/conformer/assets/101199234/269a83f8-c337-489d-832b-901f3a162854)
* ... Or you can run an external calculation through Gaussian (you need to set the path to the Gaussian executable in 'Settings'). This has a few extra options, method, number of CPU cores, ammount of RAM. If you select an "opt" job then the conformers will be optimised with your chosen method (after the RMS pruning etc.)
* ![image](https://github.com/RichardMandle/conformer/assets/101199234/7874985f-c73e-4313-816b-629c05f9cfa7)
* On the "Analysis #1" tab we can visualise the first of our results as a histogram of probability of a given bend angle. A slider controls the bin width of the histogram, and you can specify the temperature (in K) used in the probability calculation. You can fit a Gaussian to the data, save the data (or histogram data) as .csv, or even just straight to .png. Some basic stats are also printed (mean angle, median angle etc.):
![image](https://github.com/RichardMandle/conformer/assets/101199234/6039b39b-a79a-47e3-87b5-0ed551c3ea6c)
* On the "Analysis #2" tab we can plot various 3d descriptors from rdkit against bend angle (useful) or against each other (less useful I think). Basic options control ploting, saving to .csv (I'll add .png later), and optional fitting with a 1,2,3,4 order polynomial:
* ![image](https://github.com/RichardMandle/conformer/assets/101199234/9df0db9a-c73e-4673-82bf-d6575c55846c)
* In the "View Conformers" tab we can do a quick inspection of the 3D geometry of our various conformers. A slider lets us select different conformers, and we can save the information to .sdf file for later retrieval. There is an interface with PyMol (which is fairly slow, slower than just loading the .sdf file into PyMol) which makes nice images. Check boxes for raytracing, grid_image and default image style (UoL-SMP) change the appearance of the resulting image. a box has some information on the selected conformer
*![image](https://github.com/RichardMandle/conformer/assets/101199234/a49a0959-666c-4f8f-a86b-ac6a8e53238a)
* The PyMol interface can make nice images; of course, this is an interactive process, but to give an idea of the sort of output the 'default' options (raytrace=True,grid_image=True,UoL-SMP=True) give see below for the LC trimer CBO6B6OCB:
* ![image](https://github.com/RichardMandle/conformer/assets/101199234/e28d11ce-0466-46ed-bbd9-a7823529e810)

Voila!

Any comments/questions/requests email r<dot>mandle<at>leeds<dot>ac<dot>uk



