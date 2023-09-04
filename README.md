![image](https://github.com/RichardMandle/conformer/assets/101199234/f3f93b59-ba1c-43b6-9bca-8e42b25973c3)# conformer
Repo for LC conformational analysis tool.
<br><br>
Current version - V 0.1:<br> backend code working <br> initial GUI working<br>
Immediate plans:<br>improve plotting tab; add options for fits, temperature control etc.<br> add option for systematic rotor search (needs new selection tool)
<br> improve handling of Gaussian files, read output geometries incase the user calls "OPT"
<br> implement HPC usage case; write .gjf files, build SGE script (.sh), execute through SGE etc.
<br> implement reading external .log files
<br><br>
## Usage
* Install dependent packages with conda (rdkit, numpy, matplotlib, tkinter, tqdm, others????) <br>
* Launch the gui with $python gui.py <br>
![image](https://github.com/RichardMandle/conformer/assets/101199234/60019c39-1fda-4d30-9aca-585cdaa89a9e)
* Load your molecule as a .mol file (can be saved from ChemDraw) <br>
* Time to select the atoms that define the two vectors we'll use to calcualte the angle. You can enter the SMILES or SMARTS strings of representative groups:<br>
* ![image](https://github.com/RichardMandle/conformer/assets/101199234/b8270dd9-8c1c-4df5-943b-aea6201d35f8)
* ... Or you can click on atoms on the first tab to define the vectors (the orientation of the vectors WILL affect the angle, currently. You want to select four atoms, these being two pairs of atoms that define a vector, making sure you click the "furthest" atom first for each pair (e.g. 'N' of CN). These get highlighted by different colours. Clicked the wrong one? If you click a 5th atom it'll just reset:
* ![image](https://github.com/RichardMandle/conformer/assets/101199234/ec2a0521-ae0f-42fa-812c-7c4d49b20e80)
* In the "conformer search" tab set your desired options. Number of conformers to screen, the desired method used to generate conformers, the title of your job, and the evaluation method. We evaluate the energy of each conformer; you can use the inbuilt "MMFF" method which is very fast:
* ![image](https://github.com/RichardMandle/conformer/assets/101199234/bbb90fb4-42cf-4f20-8e65-89a81e0331ad)
* ... Or you can run an external calculation through Gaussian (you need to set the path to the Gaussian executable in 'Settings'). This has a few extra options, method, number of CPU cores, ammount of RAM. Don't run an OPT job as the code doesn't (currently) check for changes in geometry:
* ![image](https://github.com/RichardMandle/conformer/assets/101199234/7874985f-c73e-4313-816b-629c05f9cfa7)
* On the analysis tab we can visualise the results as a histogram. The slider controls the number of bins (or is it bin width?). I'll add more function here soon.
![Uploading image.png…]()
<br><br>
Voila!




