import tkinter as tk
from tkinter import ttk, filedialog

from tqdm import tqdm
import numpy as np

import csv
import json
import io
from collections import defaultdict

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

from PIL import ImageTk, Image

from rdkit import Chem
from rdkit.Chem import Draw, AllChem, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

#import our own modules
import draw
import conformer

class ConformerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Conformer GUI")
        self.notebook = ttk.Notebook(self.root)

        self.tab1 = ttk.Frame(self.notebook)
        self.tab2 = ttk.Frame(self.notebook)
        self.tab3 = ttk.Frame(self.notebook)
        self.tab4 = ttk.Frame(self.notebook)

        self.notebook.add(self.tab1, text="Load Molecule")
        self.notebook.add(self.tab2, text="Conformer Search")
        self.notebook.add(self.tab3, text="Analysis")
        self.notebook.add(self.tab4, text="Settings")

        self.notebook.pack(expand=True, fill="both")
        
        #
        self.canvas = tk.Canvas(self.tab1, width=800, height=300)
        self.canvas.pack()

        self.canvas.bind("<Button-1>", self.get_nearest_atom)
        
        # load the last saved path settings
        entry_gaussian = tk.Entry(self.tab4, width=40)
        entry_pymol = tk.Entry(self.tab4, width=40)
                
        #mol/image options
        self.loaded_mol = None
        self.drawer = None
        self.highlights = []
        self.settings = {}
        
        
        # set label options
        self.selected_label = tk.Label(root, text="Selected Atoms:")
        self.selected_label.pack()
        self.atom_label = tk.Label(root, text="")
        self.atom_label.pack()
        self.selected_atoms_label = tk.Label(root, text="")
        self.selected_atoms_label.pack()
        
        # Load and update settings from JSON
        self.gaussian_path, self.pymol_path = self.load_settings()  # Pass the entries to the method
        if entry_gaussian:
            print('Gaussian Path is: ' + self.gaussian_path)
        if entry_pymol:
            print('PyMol Path is: ' + self.pymol_path)
            
        #create tab contents last
        self.create_tab1_content()
        self.create_tab2_content()
        self.create_tab3_content()
        self.create_tab4_content(entry_gaussian, entry_pymol)  # Pass the entries to the method

    def create_tab1_content(self):
        self.mol_image = tk.Label(self.tab1, text=" \n ")
        self.mol_image.pack()

        open_button = tk.Button(self.tab1, text="Open File Dialog", command=self.open_file_dialog)
        open_button.pack()

        # Additional label to display various information
        self.smiles_label = tk.Label(self.tab1, text=" ")
        self.smiles_label.pack()
        
        self.rotatable_bonds_label = tk.Label(self.tab1, text=" ")
        self.rotatable_bonds_label.pack()
        
    def draw_molecule(self):
        if self.loaded_mol is not None:
            # for debug
            #for i, atom in enumerate(self.loaded_mol.GetAtoms()):
            #    print(str(atom.GetSymbol()) +' ' + str(atom.GetIdx()))
                
            AllChem.Compute2DCoords(self.loaded_mol)
            self.drawer = rdMolDraw2D.MolDraw2DCairo(800, 300)
            
            if len(self.highlights) == 0:
                self.drawer.DrawMolecule(self.loaded_mol)
             
            if len(self.highlights) > 0:
                colors = [(0.0039, 0.8, 0.7176, 0.3),  # define some colours to use
                          (0.9843, 0.0078, 0.4980, 0.3), 
                          (1.0, 0.9843, 0.0, 0.3), 
                          (0.2235, 0.2353, 1.0, 0.3)]
                          
                athighlights = defaultdict(list)
                arads = {}
                
                for a in self.loaded_mol.GetAtoms():
                    if a.GetIdx() in [x -1 for x in self.highlights]:
                        aid = a.GetIdx()
                        athighlights[aid].append(colors[[x -1 for x in self.highlights].index(a.GetIdx())])
                        arads[aid] = 0.3
                        
                self.drawer.DrawMoleculeWithHighlights(self.loaded_mol,"",dict(athighlights),{},arads,{})
                
            self.drawer.FinishDrawing()
            
            img = self.drawer.GetDrawingText()
            img = Image.open(io.BytesIO(img))
            self.photo = ImageTk.PhotoImage(img)

            self.canvas.create_image(0, 0, image=self.photo, anchor=tk.NW) # Display the above image on the canvas
            
            # Update selected atoms label
            selected_atoms_text = ""
            for i, atom_idx in enumerate(self.highlights):
                selected_atoms_text += f"Atom{i+1}: index {atom_idx-1} / type {self.loaded_mol.GetAtomWithIdx(atom_idx-1).GetSymbol()} \n" 
                
            self.selected_atoms_label.config(text=selected_atoms_text)     

    def get_nearest_atom(self, event):
        if self.loaded_mol is None:
            return

        x, y = event.x, event.y
        min_distance = float(10) # set a cutoff for how far away we can be (in pixels??)
        nearest_atom = None

        for i in range(len(self.loaded_mol.GetAtoms())):
            drawn_atom_coords = self.drawer.GetDrawCoords(i)
            atom_x, atom_y = drawn_atom_coords.x, drawn_atom_coords.y
            distance = np.sqrt((atom_x - x)**2 + (atom_y - y)**2)
            if distance < min_distance:
                min_distance = distance
                nearest_atom = i

        if nearest_atom is not None:
            self.atom_label.config(text=f"Last Selected Atom: {nearest_atom + 1}")
            
            self.highlights.append(nearest_atom+1)
            if len(self.highlights) == 5: # reset list if it has 5 entries
                self.highlights = [nearest_atom+1]
                
            self.draw_molecule()
            
    def open_file_dialog(self):
        filetypes = (
            ("MOL ", "*.mol"),
            ("MOL2", "*.mol2"),
            ("PDB", "*.pdb"),
            ("smiles", "*.smi"),
            ("smarts", "*.sma"),
            ("Text files", "*.txt"),
            ("All files", "*.*")
        )
        file_path = filedialog.askopenfilename(filetypes=filetypes)

        if file_path:
            try:
                mol = None
                file_extension = file_path.split(".")[-1].lower()

                if file_extension == "mol":
                    mol = Chem.MolFromMolFile(file_path)
                elif file_extension == "mol2":
                    mol = Chem.MolFromMol2File(file_path)
                elif file_extension == "pdb":
                    mol = Chem.MolFromPDBFile(file_path)
                elif file_extension == "smi":
                    with open(file_path, "r") as f:
                        smi = f.readline().strip()
                        mol = Chem.MolFromSmiles(smi)
                elif file_extension == "sma":
                    with open(file_path, "r") as f:
                        smi = f.readline().strip()
                        mol = Chem.MolFromSmarts(smi)
                elif file_extension == "txt":
                    with open(file_path, "r") as f:
                        mol = Chem.MolFromMolBlock(f.read())

                if mol is not None:
                    self.loaded_mol = Chem.AddHs(mol)  # Store the loaded mol object and add hydrogens
                    print("Molecule successfully loaded!\n")

                    # Update the tab label with some information
                    smiles = Chem.MolToSmiles(mol)
                    self.smiles_label.config(text="SMILES: " + smiles)
                    
                    rot_bonds = conformer.count_rot_bonds(smiles)
                    self.rotatable_bonds_label.config(text= str(rot_bonds) + ' Rotatble Bonds')
                    
                    est_num_confs = conformer.count_rot_bonds(smiles)
                    self.rotatable_bonds_label.config(text= str(3**est_num_confs) + ' Estimated conformers to screen')
                    
                    self.draw_molecule()
                
            except Exception as e:
                print("Error:", str(e))
                
    ## TAB 2 ##
    def create_tab2_content(self):
        tab2_label = tk.Label(self.tab2, text="Conformer Search Options")
        tab2_label.pack(pady=20)

        # Create a grid for the settings
        vector_settings_grid = ttk.LabelFrame(self.tab2, text="Angle Definition Settings")
        vector_settings_grid.pack(padx=10, pady=10, fill="both", expand=True)
        
        # Create a dropdown menu for angle definition method
        angle_method_label = tk.Label(vector_settings_grid, text="Angle Definition Method:")
        angle_method_label.grid(row=0, column=0, padx=10, pady=5, sticky="e")

        self.angle_methods = ["by atom index", "by SMILES match", "by SMARTS match"]
        self.selected_angle_method = tk.StringVar()
        
        angle_method_dropdown = ttk.Combobox(vector_settings_grid, textvariable=self.selected_angle_method, values=self.angle_methods)
        angle_method_dropdown.grid(row=0, column=1, padx=10, pady=5, sticky="w")
        angle_method_dropdown.set(self.angle_methods[0])
        
        # Initialize SMILES and SMARTS entry widgets
        smiles_group1_label = tk.Label(vector_settings_grid, text="SMILES Group 1:")
        self.smiles_group1_entry = tk.Entry(vector_settings_grid, width=30)
        smiles_group2_label = tk.Label(vector_settings_grid, text="SMILES Group 2:")
        self.smiles_group2_entry = tk.Entry(vector_settings_grid, width=30)

        smarts_group1_label = tk.Label(vector_settings_grid, text="SMARTS Group 1:")
        self.smarts_group1_entry = tk.Entry(vector_settings_grid, width=30)
        smarts_group2_label = tk.Label(vector_settings_grid, text="SMARTS Group 2:")
        self.smarts_group2_entry = tk.Entry(vector_settings_grid, width=30)
        
        def update_settings():
            angle_method = self.selected_angle_method.get()
            print(angle_method)
            # Hide all widgets
            for widget in vector_settings_grid.winfo_children():
                widget.grid_remove()

            # Show angle method dropdown
            angle_method_label.grid(row=0, column=0, padx=10, pady=5, sticky="e")
            angle_method_dropdown.grid(row=0, column=1, padx=10, pady=5, sticky="w")

            # Display appropriate labels and entries based on angle method
            if angle_method == "by atom index":
                if len(self.highlights) == 4:
                    # Show highlights information
                    highlights_label = tk.Label(vector_settings_grid, text=f'Atom group #1: {self.highlights[0]}, {self.highlights[1]}')
                    highlights_label.grid(row=1, column=0, columnspan=2, padx=10, pady=5)
    
                    highlights_label2 = tk.Label(vector_settings_grid, text=f'Atom group #2: {self.highlights[2]}, {self.highlights[3]}')
                    highlights_label2.grid(row=2, column=0, columnspan=2, padx=10, pady=5)
                
                if len(self.highlights) != 4:
                    highlights_label = tk.Label(vector_settings_grid, text='Please select atoms using the tool on tab #1')
                    highlights_label.grid(row=1, column=0, columnspan=2, padx=10, pady=5)
            
            elif angle_method == "by SMILES match":
                smiles_group1_label.grid(row=1, column=0, padx=10, pady=5, sticky="e")
                self.smiles_group1_entry.grid(row=1, column=1, padx=10, pady=5, sticky="w")

                smiles_group2_label.grid(row=2, column=0, padx=10, pady=5, sticky="e")
                self.smiles_group2_entry.grid(row=2, column=1, padx=10, pady=5, sticky="w")

            elif angle_method == "by SMARTS match":
                smarts_group1_label.grid(row=1, column=0, padx=10, pady=5, sticky="e")
                self.smarts_group1_entry.grid(row=1, column=1, padx=10, pady=5, sticky="w")

                smarts_group2_label.grid(row=2, column=0, padx=10, pady=5, sticky="e")
                self.smarts_group2_entry.grid(row=2, column=1, padx=10, pady=5, sticky="w")
            
        self.selected_angle_method.trace("w", lambda *args: update_settings())  # Set up the trace to call update_settings when the angle definition method changes

        update_settings()  # Call update_settings function initially
        
        settings_grid = ttk.LabelFrame(self.tab2, text="Conformer Search Settings")
        settings_grid.pack(padx=10, pady=10, fill="both", expand=True)
        # Create labels and entry boxes for different settings
        settings_labels = [ "Number of Conformers:", "Job Name:"]
        self.settings_entries = []
        for i, label_text in enumerate(settings_labels):
            label = tk.Label(settings_grid, text=label_text)
            label.grid(row=i , column=0, padx=10, pady=5, sticky="e")

            entry = tk.Entry(settings_grid, width=30)
            entry.grid(row=i , column=1, padx=10, pady=5, sticky="w")
            self.settings_entries.append(entry)

        # Create a dropdown menu for the conformer and evaluation methods selection
        method_label = tk.Label(settings_grid, text="Conformer Method:")
        method_label.grid(row=len(settings_labels) , column=0, padx=10, pady=5, sticky="e")

        self.methods = conformer.get_conformer_methods()
        self.selected_method = tk.StringVar()

        method_dropdown = ttk.Combobox(settings_grid, textvariable=self.selected_method, values=self.methods)
        method_dropdown.grid(row=len(settings_labels) , column=1, padx=10, pady=5, sticky="w")
        method_dropdown.set(self.methods[0])

        eval_method_label = tk.Label(settings_grid, text="Evaluation Method:")
        eval_method_label.grid(row=len(settings_labels) + 1 , column=0, padx=10, pady=5, sticky="e")

        self.eval_methods = ["MMFF (internal)", "Gaussian (external)"]
        self.selected_eval_method = tk.StringVar()

        eval_method_dropdown = ttk.Combobox(settings_grid, textvariable=self.selected_eval_method, values=self.eval_methods)
        eval_method_dropdown.grid(row=len(settings_labels) + 1 , column=1, padx=10, pady=5, sticky="w")
        eval_method_dropdown.set(self.eval_methods[0])

        # Create additional settings widgets for Gaussian options
        self.additional_settings_widgets = []
        options_label = tk.Label(settings_grid, text="Gaussian Job Options:")
        options_entry = tk.Entry(settings_grid, width=30)
        cores_label = tk.Label(settings_grid, text="Number of CPU cores:")
        cores_entry = tk.Entry(settings_grid, width=30)
        ram_label = tk.Label(settings_grid, text="RAM / GB:")
        ram_entry = tk.Entry(settings_grid, width=30)

        def update_additional_settings(*args):
            # Clear the previous additional settings widgets
            for widget in self.additional_settings_widgets:
                widget.grid_forget()

            # Check if the selected method is "Gaussian (external)"
            if self.selected_eval_method.get() == "Gaussian (external)":
                options_label.grid(row=len(settings_labels) + 2, column=0, padx=10, pady=5, sticky="e")
                options_entry.grid(row=len(settings_labels) + 2, column=1, padx=10, pady=5, sticky="w")
                self.additional_settings_widgets.extend([options_label, options_entry])

                cores_label.grid(row=len(settings_labels) + 3, column=0, padx=10, pady=5, sticky="e")
                cores_entry.grid(row=len(settings_labels) + 3, column=1, padx=10, pady=5, sticky="w")
                self.additional_settings_widgets.extend([cores_label, cores_entry])

                ram_label.grid(row=len(settings_labels) + 4, column=0, padx=10, pady=5, sticky="e")
                ram_entry.grid(row=len(settings_labels) + 4, column=1, padx=10, pady=5, sticky="w")
                self.additional_settings_widgets.extend([ram_label, ram_entry])

                # Add buttons for writing Gaussian input and reading output
                write_gaussian_input_button = tk.Button(settings_grid, text="Write Gaussian Input Only", command=print('to do #1'))
                write_gaussian_input_button.grid(row=len(settings_labels) + 5, column=0, padx=10, pady=5, sticky="e")

                read_gaussian_output_button = tk.Button(settings_grid, text="Read Gaussian Output Files", command=print('to do #2'))
                read_gaussian_output_button.grid(row=len(settings_labels) + 5, column=1, padx=10, pady=5, sticky="w")

                self.additional_settings_widgets.extend([write_gaussian_input_button, read_gaussian_output_button])

        # Set up the trace to call update_additional_settings when the evaluation method changes
        self.selected_eval_method.trace("w", update_additional_settings)
        
        # Create a button for advanced settings
        advanced_settings_button = tk.Button(settings_grid, text="Advanced Settings", command=self.show_advanced_settings)
        advanced_settings_button.grid(row=len(settings_labels) + 7, columnspan=2, pady=10)

        # Create a button to trigger the conformer search
        search_button = tk.Button(settings_grid, text="Start Conformer Search", command=self.start_conformer_search)
        search_button.grid(row=len(settings_labels) + 6, columnspan=2, pady=20)
        
    def show_advanced_settings(self):
        # Create a popup window for advanced settings
        advanced_settings_window = tk.Toplevel(self.root)
        advanced_settings_window.title("Advanced Settings")

        # Create and layout the advanced settings widgets
        options_frame = ttk.LabelFrame(advanced_settings_window, text="Advanced Options")
        options_frame.grid(row=0, column=0, padx=10, pady=10)

        # Create and add the widgets for the advanced options
        advanced_options = [
            ("RMS Threshold for Pruning", "0.33", "0.33"),
            ("Random Seed", "61453", "61453"),
            ("Use Torsion Preferences", ["True", "False"], "True"),
            ("Use Chemical Knowledge", ["True", "False"], "True"),
            ("Use Random Coords", ["True", "False"], "False"),
            ("Do MMFF Optimization", ["True", "False"], "False"),
        ]

        mmff_options_frame = None  # Initialize as None
        option_vars = {}  # Create a dictionary to store option variables
        mmff_option_vars = {}  # Create a dictionary to store MMFF option variables

        for row, (label_text, values, default_value) in enumerate(advanced_options):
            label = tk.Label(options_frame, text=(label_text + ':'))
            label.grid(row=row, column=0, padx=10, pady=5, sticky="e")

            option_var = None  # Initialize option_var as None

            if isinstance(values, list):
                option_var = tk.StringVar(value=default_value)
                option_menu = ttk.Combobox(options_frame, textvariable=option_var, values=values)
                option_menu.grid(row=row, column=1, padx=10, pady=5, sticky="w")
            else:
                option_var = tk.StringVar(value=default_value)
                entry = tk.Entry(options_frame, textvariable=option_var, width=30)
                entry.grid(row=row, column=1, padx=10, pady=5, sticky="w")

            # Store the option variable in the dictionary
            option_vars[label_text] = option_var

            if label_text == "Do MMFF Optimization:":
                def show_mmff_options():
                    nonlocal mmff_options_frame
                    if mmff_options_frame is not None:
                        mmff_options_frame.destroy()
                    if option_var.get() == "True":
                        mmff_options_frame = ttk.LabelFrame(advanced_settings_window, text="MMFF Options")
                        mmff_options_frame.grid(row=row + 1, column=0, padx=10, pady=10)

                        mmff_advanced_options = [
                            ("Energy Threshold (kcal/mol)", "100"),
                            ("Max Iterations", "500"),
                        ]

                        for mmff_row, (mmff_label_text, mmff_default_value) in enumerate(mmff_advanced_options):
                            mmff_label = tk.Label(mmff_options_frame, text=mmff_label_text)
                            mmff_label.grid(row=mmff_row, column=0, padx=10, pady=5, sticky="e")

                            mmff_entry_var = tk.StringVar(value=mmff_default_value)
                            mmff_entry = tk.Entry(mmff_options_frame, textvariable=mmff_entry_var, width=30)
                            mmff_entry.grid(row=mmff_row, column=1, padx=10, pady=5, sticky="w")

                            mmff_option_vars[mmff_label_text] = mmff_entry_var

                option_var.trace_add("write", lambda *args: show_mmff_options())
                
                    
        def save_advanced_settings():
            # Retrieve and save the advanced settings here
            # Access the values of the widgets using their corresponding variables
            for label_text, option_var in option_vars.items():
                self.settings[label_text] = option_var.get()

            if mmff_options_frame:
                mmff_settings = {}
                for mmff_label_text, mmff_entry_var in mmff_option_vars.items():
                    mmff_settings[mmff_label_text] = mmff_entry_var.get()
                self.settings["MMFF Options"] = mmff_settings

            print("Advanced Settings:", self.settings)

        save_button = tk.Button(advanced_settings_window, text="Save Settings", command=save_advanced_settings)
        save_button.grid(row=len(advanced_options) + 1, column=0, columnspan=2, pady=10)
        
        
    def start_conformer_search(self):
        # Gather inputs from the settings entries
        num_conformers = int(self.settings_entries[0].get())
        job_name = self.settings_entries[1].get()
        selected_method = self.selected_method.get()
        selected_eval_method = self.selected_eval_method.get()     
        
        gaussian_job = False # if its a Gaussian job then we'll set this flag to true
        gaussian_options = ""
        gaussian_nproc = ""
        gaussian_vmem = ""
        
        if selected_eval_method == self.eval_methods[-1]: # if its the gaussian method then do this
            gaussian_job = True 
            gaussian_options = self.additional_settings_widgets[1].get()
            gaussian_nproc = self.additional_settings_widgets[3].get()
            gaussian_vmem = self.additional_settings_widgets[5].get()
        
        rms_threshold = float(0.33) # TO DO - this is turning into a tuple, why?? TO DO
        use_torsion_pref = True
        use_knowledge = True
        use_random_coords = False
        random_seed = int(61453)
        opt = False
        max_opt_iter=5000
        min_energy_MMFF=500
            
        if self.settings:
            print('Using Advanced Settings')
            rms_threshold = float(self.settings['RMS Threshold for Pruning'])
            use_torsion_pref = self.settings['Use Torsion Preferences'] == True
            use_knowledge = self.settings['Use Chemical Knowledge'] == True
            use_random_coords = self.settings['Use Random Coords'] == True
            random_seed = int(self.settings['Random Seed'])
            opt = self.settings['Do MMFF Optimization'] == True
            if opt == True:
                max_opt_iter=float(self.settings['Max Iterations'])
                min_energy_MMFF= float(self.settings['Energy Threshold (kcal/mol)'])
        
        atoms=[]
        vec1 = ''
        vec2 = ''
        print(self.selected_angle_method.get())
        if self.selected_angle_method.get() == 'by atom index':
            vector_definition_method = 'atoms'
            atoms = [x-1 for x in self.highlights] #remember, rdkit indexes from zero...

        if self.selected_angle_method.get() == 'by SMILES match':
            vector_definition_method = 'smi'
            vec1 = self.smiles_group1_entry.get()
            vec2 = self.smiles_group2_entry.get()
            
        if self.selected_angle_method.get() == 'by SMARTS match':
            vector_definition_method = 'sma'
            vec1 = self.smarts_group1_entry.get()
            vec2 = self.smarts_group2_entry.get()

        # Call conf.conf_search with the gathered inputs
        # need to update with additional parameters from advanced settings
        self.mol_conf,self.energy,self.angle,name = conformer.conf_search(
            mol=self.loaded_mol,  
            vec_def_method=vector_definition_method,
            atom_idx_list=atoms,
            g_path = self.gaussian_path,
            embeded_method=selected_method,
            vector1=vec1,
            vector2=vec2,
            n_conf=num_conformers,
            name=job_name,
            rms_thresh = rms_threshold,
            use_torsion_pref = use_torsion_pref,
            use_knowledge = use_knowledge,
            use_random_coords = use_random_coords,
            random_seed = random_seed,
            opt = opt,
            max_iter= max_opt_iter,
            min_energy_MMFF= min_energy_MMFF,
            Gauss=gaussian_job,
            options=(gaussian_options),
            cores=(gaussian_nproc),
            ram=(gaussian_vmem),
            write_only = False # TO DO - this is hard coded, remember to expose it!!
        )
        if self.angle:
            print('Successful completion of conformer search')
            print('Number of conformers screened: ' + str(len(self.angle)))
        
    ## TAB 3 ##
    def create_tab3_content(self):
        tab3_label = tk.Label(self.tab3, text="Histograms and Plots")
        tab3_label.pack(pady=20)

        # Create canvas widgets for histograms and plots
        self.canvas_histogram = tk.Canvas(self.tab3, width=250, height=200)
        self.canvas_histogram.pack(side=tk.LEFT, padx=10, pady=10)

        self.canvas_radii = tk.Canvas(self.tab3, width=250, height=200)
        self.canvas_radii.pack(side=tk.LEFT, padx=10, pady=10)

        bin_label = tk.Label(self.tab3, text="Width of Histogram Bins:")
        bin_label.pack()

        self.bin_slider = tk.Scale(self.tab3, from_=1, to=100, orient="horizontal")
        self.bin_slider.set(10)  # Set an initial value
        self.bin_slider.pack()

        # Entry widget for temperature input
        temp_label = tk.Label(self.tab3, text="Temperature / K:")
        temp_label.pack()
        self.temp_entry = tk.Entry(self.tab3)
        self.temp_entry.insert(0, "298")  # default to 298 K, why not!
        self.temp_entry.pack()

        # Dropdown menu for fitting Gaussian
        fit_gaussian_label = tk.Label(self.tab3, text="Fit Gaussian?")
        fit_gaussian_label.pack()
        self.fit_gaussian_var = tk.StringVar(self.tab3)
        self.fit_gaussian_var.set("Yes")  # Set an initial value
        fit_gaussian_option = tk.OptionMenu(self.tab3, self.fit_gaussian_var, "Yes", "No")
        fit_gaussian_option.pack()

        # Button to generate histograms and plots
        generate_button = tk.Button(self.tab3, text="Generate Histograms and Plots", command=self.generate_histograms_and_plots)
        generate_button.pack()

        # Button to save data as CSV
        save_data_button = tk.Button(self.tab3, text="Save Data as CSV", command=self.save_data_as_csv)
        save_data_button.pack()

        # Button to save the plot as PNG
        save_plot_button = tk.Button(self.tab3, text="Save Plot as PNG", command=self.save_plot_as_png)
        save_plot_button.pack()

    def generate_histograms_and_plots(self):
        num_bins = self.bin_slider.get()  # Get the selected number of bins from the slider
        temperature = float(self.temp_entry.get())  # Get user-input temperature
        fit_gaussian = (self.fit_gaussian_var.get() == "Yes")  # Check fit Gaussian option

        #deltaE_KJ, probability, x, y 
        self.returned_fig = conformer.calc_bend_hist(self.angle, 
                                                self.energy, 
                                                fit_gaussian=fit_gaussian, 
                                                name=self.settings_entries[1].get(), 
                                                Temp=temperature, 
                                                BinSteps=num_bins)

        self.returned_radii_fig = conformer.rad_gyr_gen(self.mol_conf,
                                                        self.energy, 
                                                        self.angle, 
                                                        name=self.settings_entries[1].get(),
                                                        show_plot=False)


        # Destroy the previous canvas widgets to remove the previous plots
        for widget in self.canvas_histogram.winfo_children():
            widget.destroy()

        for widget in self.canvas_radii.winfo_children():
            widget.destroy()

        canvas_img = FigureCanvasTkAgg(self.returned_fig, master=self.canvas_histogram)
        canvas_img.get_tk_widget().pack()
        
        canvas_img = FigureCanvasTkAgg(self.returned_radii_fig, master=self.canvas_radii)
        canvas_img.get_tk_widget().pack()

        # Plot radius of gyration vs angle on the canvas_plot
        # to do

    def save_data_as_csv(self):
        data = [['Angle / Deg', [angle for angle in self.angle]], ['dE / Kj Mol', [energy for energy in self.energy]]]

        # Create a file dialog for saving the CSV file
        file_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV Files", "*.csv")])

        if file_path:
            with open(file_path, 'w', newline='') as csvfile:
                csvwriter = csv.writer(csvfile)
                csvwriter.writerows(data)
            
    def save_plot_as_png(self):
        file_path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG Image", "*.png")])
        if file_path:
            self.returned_fig.savefig(file_path)


    def create_tab4_content(self, entry_gaussian, entry_pymol):

        label_gaussian = tk.Label(self.tab4, text="Path to Gaussian:")
        entry_gaussian = tk.Entry(self.tab4, width=40)
        button_browse_gaussian = tk.Button(self.tab4, text="Browse",
                                            command=lambda: self.browse_gaussian(entry_gaussian))
        
        label_pymol = tk.Label(self.tab4, text="Path to PyMOL:")
        entry_pymol = tk.Entry(self.tab4, width=40)
        button_browse_pymol = tk.Button(self.tab4, text="Browse",
                                        command=lambda: self.browse_pymol(entry_pymol))
        
        label_gaussian.pack(pady=10)
        entry_gaussian.pack(pady=5)
        button_browse_gaussian.pack()
        
        label_pymol.pack(pady=10)
        entry_pymol.pack(pady=5)
        button_browse_pymol.pack()
        
        button_save = tk.Button(self.tab4, text="Save Settings",
                                command=lambda: self.save_settings(entry_gaussian.get(), entry_pymol.get()))
        button_save.pack(pady=10)
        
        button_load = tk.Button(self.tab4, text="Load Settings",
                                command=lambda: self.load_settings_and_update_ui(entry_gaussian, entry_pymol))
        button_load.pack(pady=5)

    def browse_gaussian(self, entry):
        file_path = filedialog.askopenfilename(filetypes=[("Executable files", "*.exe")])
        if file_path:
            entry.delete(0, tk.END)
            entry.insert(0, file_path)

    def browse_pymol(self, entry):
        file_path = filedialog.askopenfilename(filetypes=[("Executable files", "*.exe")])
        if file_path:
            entry.delete(0, tk.END)
            entry.insert(0, file_path)  

    def save_settings(self, gaussian_path, pymol_path):
        settings = {
            "gaussian": gaussian_path,
            "pymol": pymol_path
        }
        with open("settings.json", "w") as f:
            json.dump(settings, f)

    def load_settings(self):
        try:
            with open("settings.json", "r") as f:
                settings = json.load(f)
                return settings.get("gaussian", ""), settings.get("pymol", "")
        except FileNotFoundError:
            return "", ""

    def load_settings_and_update_ui(self, entry_gaussian, entry_pymol):
        gaussian_path, pymol_path = self.load_settings()
        entry_gaussian.delete(0, tk.END)
        entry_gaussian.insert(0, gaussian_path)
        entry_pymol.delete(0, tk.END)
        entry_pymol.insert(0, pymol_path)
        
        
if __name__ == "__main__":
    root = tk.Tk()
    app = ConformerGUI(root)
    root.mainloop()
