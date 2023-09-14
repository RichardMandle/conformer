import tkinter as tk
from tkinter import ttk, filedialog, OptionMenu, BooleanVar, simpledialog  

from tqdm import tqdm
import numpy as np

import csv
import json
import io
import os
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
        self.tab5 = ttk.Frame(self.notebook)
        self.settings_tab = ttk.Frame(self.notebook)

        self.notebook.add(self.tab1, text="Load Molecule")
        self.notebook.add(self.tab2, text="Conformer Search")
        self.notebook.add(self.tab3, text="Analysis #1")
        self.notebook.add(self.tab4, text="Analysis #2")
        self.notebook.add(self.tab5, text="View Conformers")
        self.notebook.add(self.settings_tab, text="Settings")
        
        self.notebook.bind("<<NotebookTabChanged>>", self.on_tab_change)
        
        self.notebook.pack(expand=True, fill="both")
        
        self.canvas = tk.Canvas(self.tab1, width=500, height=300)
        self.canvas.pack()
        self.canvas.bind("<Button-1>", self.get_nearest_atom)
        
        # load the last saved path settings
        entry_gaussian = tk.Entry(self.settings_tab, width=40)
        entry_pymol = tk.Entry(self.settings_tab, width=40)
                
        #mol/image options
        self.loaded_mol = None
        self.property_dict = None
        self.drawer = None
        self.highlights = []
        self.settings = {}
        
        self.probability = []
        self.angle = []
        self.energy = []
        
        #Define BooleanVar variables for pymol
        self.pymol_ray_trace_var = tk.BooleanVar()
        self.pymol_ray_trace_var.set(False)
        self.pymol_draw_grid_var = tk.BooleanVar()
        self.pymol_draw_grid_var.set(True)
        self.pymol_use_uol_smp_var = tk.BooleanVar()
        self.pymol_use_uol_smp_var.set(True)
        
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
        self.create_tab4_content()
        self.create_tab5_content()
        self.generate_and_display_first_conformer()
        self.create_settings_tab_content(entry_gaussian, entry_pymol)  # Pass the entries to the method
    
    def on_tab_change(self, event):
        # Check if the current tab is the one you want to trigger the canvas update
        if self.notebook.index(self.notebook.select()) == 2:  # Tab 3
            self.generate_histograms_and_plots()  # Call the function to update canvas
        
        if self.notebook.index(self.notebook.select()) == 3:  # Tab 4
            # Recreate tab4 content if self.loaded_mol is available
            if hasattr(self, 'loaded_mol') and self.loaded_mol is not None:
                # Clear the tab4 frame to recreate the content
                for widget in self.tab4.winfo_children():
                    widget.destroy()
                self.create_tab4_content()
            else:
                # Display a message if self.loaded_mol is not available
                self.clear_tab4_content()
                print("No mol object found")       
                
        if self.notebook.index(self.notebook.select()) == 4:  # Tab 4
            self.generate_and_display_first_conformer()  # Call the function to update the conformer canvas
    
    
    def create_tab1_content(self):
        self.mol_image = tk.Label(self.tab1, text=" \n ")
        self.mol_image.pack()

        open_button = tk.Button(self.tab1, text="Open File Dialog", command=self.open_file_dialog)
        open_button.pack()

    def draw_molecule(self):
        if self.loaded_mol is not None:
                
            AllChem.Compute2DCoords(self.loaded_mol)
            self.drawer = rdMolDraw2D.MolDraw2DCairo(500, 300)
            
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
            #self.atom_label.config(text=f"Last Selected Atom: {nearest_atom + 1}")
            
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
                    self.draw_molecule()
                
            except Exception as e:
                print("Error:", str(e))
                
    ## TAB 2 ##
    def create_tab2_content(self):
        tab2_label = tk.Label(self.tab2, text="Conformer Search Options")
        tab2_label.pack(pady=5)

        vector_settings_grid = ttk.LabelFrame(self.tab2, text="Angle Definition Settings")
        vector_settings_grid.pack(padx=10, pady=5, fill="both", expand=True)
        
        angle_method_label = tk.Label(vector_settings_grid, text="Angle Definition Method:")
        angle_method_label.grid(row=0, column=0, padx=10, pady=5, sticky="e")

        self.angle_methods = ["by atom index", "by SMILES match", "by SMARTS match","None"]
        self.selected_angle_method = tk.StringVar()
        
        angle_method_dropdown = ttk.Combobox(vector_settings_grid, textvariable=self.selected_angle_method, values=self.angle_methods)
        angle_method_dropdown.grid(row=0, column=1, padx=10, pady=5, sticky="w")
        angle_method_dropdown.set(self.angle_methods[0])
        
        smiles_group1_label = tk.Label(vector_settings_grid, text="SMILES Group 1:")
        self.smiles_group1_entry = tk.Entry(vector_settings_grid, width=30)
        smiles_group2_label = tk.Label(vector_settings_grid, text="SMILES Group 2:")
        self.smiles_group2_entry = tk.Entry(vector_settings_grid, width=30)

        smarts_group1_label = tk.Label(vector_settings_grid, text="SMARTS Group 1:")
        self.smarts_group1_entry = tk.Entry(vector_settings_grid, width=30)
        smarts_group2_label = tk.Label(vector_settings_grid, text="SMARTS Group 2:")
        self.smarts_group2_entry = tk.Entry(vector_settings_grid, width=30)
        
        def update_settings():
            # function to update conformer search settings.
            angle_method = self.selected_angle_method.get()
            for widget in vector_settings_grid.winfo_children():
                widget.grid_remove() # Hide all widgets on update

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
            
            elif angle_method == "None":
                print('No angle method selected; will compute energy only')
                
        self.selected_angle_method.trace("w", lambda *args: update_settings())  # Set up the trace to call update_settings when the angle definition method changes

        update_settings()  # Call update_settings function initially
        
        settings_grid = ttk.LabelFrame(self.tab2, text="Conformer Search Settings")
        settings_grid.pack(padx=10, pady=5, fill="both", expand=True)
        settings_labels = [ "Number of Conformers:", "Job Name:"]
        self.settings_entries = []
        for i, label_text in enumerate(settings_labels):
            label = tk.Label(settings_grid, text=label_text)
            label.grid(row=i , column=0, padx=10, pady=5, sticky="e")
            entry = tk.Entry(settings_grid, width=30)
            entry.grid(row=i , column=1, padx=10, pady=5, sticky="w")
            self.settings_entries.append(entry)

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

        # additional settings widgets for Gaussian options
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
        
        advanced_settings_button = tk.Button(settings_grid, text="Advanced Settings", command=self.show_advanced_settings)
        advanced_settings_button.grid(row=len(settings_labels) + 7, columnspan=2, pady=5)

        search_button = tk.Button(settings_grid, text="Start Conformer Search", command=self.start_conformer_search)
        search_button.grid(row=len(settings_labels) + 6, columnspan=2, pady=20)
        
    def show_advanced_settings(self):
        # this function makes a popup window for some more advanced settings
        advanced_settings_window = tk.Toplevel(self.root)
        advanced_settings_window.title("Advanced Settings")

        options_frame = ttk.LabelFrame(advanced_settings_window, text="Advanced Options")
        options_frame.grid(row=0, column=0, padx=10, pady=5)

        advanced_options = [
            ("RMS Threshold for Pruning", "0.33", "0.33"),
            ("Random Seed", "61453", "61453"),
            ("Use Torsion Preferences", ["True", "False"], "True"),
            ("Use Chemical Knowledge", ["True", "False"], "True"),
            ("Use Random Coords", ["True", "False"], "False"),
            ("MMFF Variant", ["MMFF94", "MMFF94s"], "MMFF94"),
            ("Do MMFF Optimization", ["True", "False"], "False"),
        ]

        mmff_options_frame = None   # Initialize as None
        option_vars = {}            # Create a dictionary to store option variables
        mmff_option_vars = {}       # Create a dictionary to store MMFF option variables

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

            if label_text == "Do MMFF Optimization":
                def show_mmff_options():
                    nonlocal mmff_options_frame
                    if mmff_options_frame is not None:
                        mmff_options_frame.destroy()
                    if option_var.get() == "True":
                        mmff_options_frame = ttk.LabelFrame(advanced_settings_window, text="MMFF Options")
                        mmff_options_frame.grid(row=row + 1, column=0, padx=10, pady=5)

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
            for label_text, option_var in option_vars.items():
                self.settings[label_text] = option_var.get()

            if mmff_options_frame:
                mmff_settings = {}
                for mmff_label_text, mmff_entry_var in mmff_option_vars.items():
                    mmff_settings[mmff_label_text] = mmff_entry_var.get()
                self.settings["MMFF Options"] = mmff_settings

        save_button = tk.Button(advanced_settings_window, text="Save Settings", command=save_advanced_settings)
        save_button.grid(row=len(advanced_options) + 1, column=0, columnspan=2, pady=5)
        
    def start_conformer_search(self):
        # Gather inputs from the settings entries
        num_conformers = int(self.settings_entries[0].get())
        job_name = self.settings_entries[1].get()
        selected_method = self.selected_method.get()
        selected_eval_method = self.selected_eval_method.get()     
        
        settings_dict = conformer.get_default_search_settings() # retrieve the search defaults from conformer.pymol
        
        if selected_eval_method == self.eval_methods[-1]: # if its the gaussian method then do this
            settings_dict['gaussian_job'] = True 
            settings_dict['gaussian_options'] = self.additional_settings_widgets[1].get()
            settings_dict['gaussian_nproc'] = self.additional_settings_widgets[3].get()
            settings_dict['gaussian_vmem'] = self.additional_settings_widgets[5].get()
      
        if self.settings:
            print('Using Advanced Settings')
            settings_dict['rms_threshold'] = float(self.settings['RMS Threshold for Pruning'])
            settings_dict['use_torsion_pref'] = self.settings['Use Torsion Preferences'] == 'True'
            settings_dict['use_knowledge'] = self.settings['Use Chemical Knowledge'] == 'True'
            settings_dict['use_random_coords'] = self.settings['Use Random Coords'] == 'True'
            settings_dict['MMFF Variant'] = self.settings["MMFF Variant"]
            settings_dict['random_seed'] = int(self.settings['Random Seed'])
            settings_dict['opt'] = ['Do MMFF Optimization'] == 'True'
            if settings_dict['opt'] == True:
                settings_dict['max_opt_iter']=float(self.settings['Max Iterations'])
                settings_dict['min_energy_MMFF']= float(self.settings['Energy Threshold (kcal/mol)'])
        
        # do the conformer search, store the conformers in self.mol_conf and energies in self.energy:
        self.mol_conf, self.energy = conformer.conf_gen(mol=self.loaded_mol,
                                                        embeded_method=selected_method,
                                                        num_of_conformer=num_conformers,
                                                        name=job_name,
                                                        rms_thresh = settings_dict['rms_threshold'],
                                                        use_torsion_pref = settings_dict['use_torsion_pref'],
                                                        use_knowledge = settings_dict['use_knowledge'],
                                                        use_random_coords = settings_dict['use_random_coords'],
                                                        random_seed = settings_dict['random_seed'],
                                                        MMFF_variant = settings_dict['MMFF Variant'],
                                                        opt = settings_dict['opt'],
                                                        max_iter= settings_dict['max_opt_iter'],
                                                        min_energy_MMFF= settings_dict['min_energy_MMFF'])
                
        atoms=[]
        vec1 = ''
        vec2 = ''

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
        
        if self.selected_angle_method.get() == 'None':
            # in the case that its 'None', just use some dummy atom coordinates.
            # The reason to do this is maybe a Gaussian calculation is requested, and currently
            # this is only allowed when an angle calculation is performed. TO DO - fix this by
            # refactoring the code to permit Gaussian jobs without angle information.
            atoms = [0,1,2,3]

        # Call conf.conf_search with the gathered inputs
        # need to update with additional parameters from advanced settings
        self.energy,self.angle = conformer.conf_analysis(
            mol_conf=self.mol_conf,
            energy = self.energy,
            vec_def_method=vector_definition_method,
            atom_idx_list=atoms,
            g_path = self.gaussian_path,
            vector1=vec1,
            vector2=vec2,
            name=job_name,
            Gauss=settings_dict['gaussian_job'],
            options=settings_dict['gaussian_options'],
            cores=settings_dict['gaussian_nproc'],
            ram=settings_dict['gaussian_vmem'],
            write_only = False # TO DO - this is hard coded, remember to expose it!!
        )
        if self.angle:
            print('\nSuccessful completion of conformer search')
            print('Number of conformers screened: ' + str(len(self.energy)))
            self.conformer_slider.config(to=(len(self.angle) - 1)) # After updating the data, update the conformer slider's range.
        
        if self.selected_angle_method.get() == 'None':
            self.angle = np.zeros_like(self.energy) # set the angles to zero, they weren't wanted!
            
    ## TAB 3 ##
    def create_tab3_content(self):
       
        tab3_label = tk.Label(self.tab3, text="Shape Analysis and Plotting")
        tab3_label.pack()

        self.canvas_histogram = tk.Canvas(self.tab3, width=200, height=200)
        self.canvas_histogram.pack(side=tk.LEFT, padx=5, pady=5)
        
        settings_frame = tk.LabelFrame(self.tab3, text="Plot Settings and Options")
        settings_frame.pack(side=tk.RIGHT, padx=5, pady=5)
        
        bin_label = tk.Label(settings_frame, text="Width of Histogram Bins:")
        bin_label.pack(pady=5)

        self.bin_slider = tk.Scale(settings_frame, from_=1, to=60, orient="horizontal")
        self.bin_slider.set(10)  # Set an initial value
        self.bin_slider.pack()

        temp_label = tk.Label(settings_frame, text="Temperature / K:")
        temp_label.pack(pady=5)
        self.temp_entry = tk.Entry(settings_frame)
        self.temp_entry.insert(0, "298")  # default to 298 K, why not!
        self.temp_entry.pack()

        fit_gaussian_label = tk.Label(settings_frame, text="Fit Gaussian?") # Gauss fit menu
        fit_gaussian_label.pack(pady=5)
        self.fit_gaussian_var = tk.StringVar(settings_frame)
        self.fit_gaussian_var.set("Yes")  # Set an initial value
        fit_gaussian_option = tk.OptionMenu(settings_frame, self.fit_gaussian_var, "Yes", "No")
        fit_gaussian_option.pack()

        save_data_button = tk.Button(settings_frame, text="Save Data as CSV", command=self.save_data_as_csv)
        save_data_button.pack()
        
        save_data_button = tk.Button(settings_frame, text="Save Hist. Data as CSV", command=self.save_hist_data_as_csv)
        save_data_button.pack()

        save_plot_button = tk.Button(settings_frame, text="Save Plot as PNG", command=self.save_plot_as_png)
        save_plot_button.pack()
        
        generate_button = tk.Button(settings_frame, text="Update Plot", command=self.generate_histograms_and_plots)
        generate_button.pack(pady=10)
        
        self.conformer_angle_stats = tk.Label(settings_frame, text="")
        self.conformer_angle_stats.pack(pady=10)


    def generate_histograms_and_plots(self):
        num_bins = self.bin_slider.get()                        # Get the selected number of bins from the slider
        temperature = float(self.temp_entry.get())              # Get user-input temperature
        fit_gaussian = (self.fit_gaussian_var.get() == "Yes")   # Check fit Gaussian option
        
        if self.angle == []:
            print('No conformer data present')
            return
            
        if self.angle != []:
            self.returned_fig,self.delta_energy,self.probability,self.hist_x,self.hist_y,xmax,fwhm = conformer.calc_bend_hist(self.angle, 
                                                    self.energy, 
                                                    fit_gaussian=fit_gaussian, 
                                                    name=self.settings_entries[1].get(), 
                                                    Temp=temperature, 
                                                    BinSteps=num_bins)

            for widget in self.canvas_histogram.winfo_children(): 
                widget.destroy() # Destroy the previous canvas widgets to remove the previous plots

            canvas_img = FigureCanvasTkAgg(self.returned_fig, master=self.canvas_histogram)
            canvas_img.get_tk_widget().pack()
            
            # Update the label with conformer information
            info_text = f"Mean Angle: {np.round(np.mean(self.angle),2)}°\nMedian Angle: {np.round(np.median(self.angle),2)}°\nFit X@Y(max): {xmax}°\nFit FWHM: {fwhm}"
            self.conformer_angle_stats.config(text=info_text)

    def save_data_as_csv(self):
        '''
        Just saves data as a CSV file for processing elsewhere.
        '''
        file_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV Files", "*.csv")])
        if file_path:
            with open(file_path, 'w', newline='') as csvfile:
                csvwriter = csv.writer(csvfile)
                csvwriter.writerow(['Angle / Deg', ' Energy / kj mol', 'dE / kj mol-1', 'Probability'])
                for angle, energy, delta_energy, probability in zip(self.angle, self.energy, self.delta_energy, self.probability):
                    csvwriter.writerow([angle, energy, delta_energy, probability])
        print('\nSaved data to ' + str(file_path))
        
    def save_hist_data_as_csv(self):
        '''
        Just saves histogram data as a CSV file for processing elsewhere.
        '''
        file_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV Files", "*.csv")])
        if file_path:
            with open(file_path, 'w', newline='') as csvfile:
                csvwriter = csv.writer(csvfile)
                csvwriter.writerow(['Angle / Deg', 'Probability'])
                for angle, energy in zip(self.hist_x, self.hist_y):
                    csvwriter.writerow([angle, energy])
        print('\nSaved data to ' + str(file_path))
                
    def save_plot_as_png(self):
        '''
        saves the current plot as .png
        '''
        file_path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG Image", "*.png")])
        if file_path:
            self.returned_fig.savefig(file_path)
            print('\nSaved figure to ' + str(file_path))
            
    ### TAB 4 ###
    def create_tab4_content(self):
        tab4_frame = ttk.Frame(self.tab4)
        tab4_frame.grid(row=0, column=0, sticky="nsew")

        fig, ax = plt.subplots(figsize=(600/plt.rcParams['figure.dpi'], 400/plt.rcParams['figure.dpi']))
        canvas = FigureCanvasTkAgg(fig, master=tab4_frame)
        canvas.get_tk_widget().grid(row=0, column=0, rowspan=2, sticky="nsew")
        
        grid1 = ttk.LabelFrame(self.tab4, text="Plots")
        grid1.grid(row=0, column=1, padx=10, pady=10, sticky="n")

        grid1_0 = ttk.LabelFrame(grid1, text="Plot Settings and Options")
        grid1_0.grid(row=0, column=0, padx=10, pady=10, sticky="w")

        grid1_1 = ttk.LabelFrame(grid1, text="Fit Settings and Options")
        grid1_1.grid(row=1, column=0, padx=10, pady=10, sticky="w")

        x_var = tk.StringVar(tab4_frame)
        y_var = tk.StringVar(tab4_frame)

        # Check if the loaded_mol object exists
        if not hasattr(self, 'loaded_mol') or self.loaded_mol is None:
            print("No mol object found")
            no_mol_label = tk.Label(tab4_frame, text="No mol object found. Please load a molecule in a previous tab.")
            no_mol_label.grid(row=0, column=0, padx=10, pady=10, sticky="n")
            return
        
        else:
            self.property_dict = conformer.get_3d_descriptors(self.loaded_mol,self.angle) # Get available properties from the property_dict
            available_properties = list(self.property_dict[0].keys())

            x_label = tk.Label(grid1_0, text="Select X Property:")
            x_label.grid(row = 0, column = 0)

            x_dropdown = tk.OptionMenu(grid1_0, x_var, *available_properties)
            x_dropdown.grid(row = 0, column = 1)

            y_label = tk.Label(grid1_0, text="Select Y Property:")
            y_label.grid(row = 1, column = 0)

            y_dropdown = tk.OptionMenu(grid1_0, y_var, *available_properties)
            y_dropdown.grid(row = 1, column = 1)
            
            x_var.set(available_properties[-1])
            y_var.set(available_properties[0])

            def plot_data():
                x_property = x_var.get()
                y_property = y_var.get()

                x_data = [self.property_dict[i][x_property] for i in range(len(self.property_dict))]
                y_data = [self.property_dict[i][y_property] for i in range(len(self.property_dict))]

                ax.clear()
                ax.scatter(x_data, y_data)
                ax.set_xlabel(x_property)
                ax.set_ylabel(y_property)
                ax.set_title(f'Scatter Plot: {x_property} vs {y_property}')
                plt.tight_layout()
                update_axis_limits()
                canvas.draw()

            plot_button = tk.Button(grid1_0, text="Plot Data", command=plot_data)
            plot_button.grid(row = 4, column = 0)

            x_axis_limit_label = tk.Label(grid1_0, text="X lim (min, max) / auto")
            x_axis_limit_label.grid(row = 2, column = 0)
            
            x_axis_limit_entry = tk.Entry(grid1_0, width=10)
            x_axis_limit_entry.grid(row = 2, column = 1)
            x_axis_limit_entry.insert(0, "auto")  # Initial X axis limits
            
            x_axis_limit_label = tk.Label(grid1_0, text="Y lim (min, max) / auto")
            x_axis_limit_label.grid(row = 3, column = 0)   
            
            y_axis_limit_entry = tk.Entry(grid1_0, width=10)
            y_axis_limit_entry.grid(row = 3, column = 1)
            y_axis_limit_entry.insert(0, "auto")  # Initial Y axis limits

            def update_axis_limits():
                x_limits = x_axis_limit_entry.get().split(',')
                y_limits = y_axis_limit_entry.get().split(',')
                
                if 'auto' not in x_limits: # let the user pass 'auto' to default to automatic axis limits
                    ax.set_xlim(float(x_limits[0]), float(x_limits[1]))
                    
                if 'auto' not in y_limits:
                    ax.set_ylim(float(y_limits[0]), float(y_limits[1]))
                
                plt.tight_layout()
                plt.draw()

            def save_csv():
                file_path = filedialog.asksaveasfilename(defaultextension='.csv', filetypes=[('CSV Files', '*.csv')])
                if file_path:
                    x_property = x_var.get()
                    y_property = y_var.get()
                    x_data = [self.property_dict[i][x_property] for i in range(len(self.property_dict))]
                    y_data = [self.property_dict[i][y_property] for i in range(len(self.property_dict))]

                    data = np.column_stack((x_data, y_data))
                    np.savetxt(file_path, data, delimiter=',', header=f'{x_property},{y_property}', comments='')

            save_csv_button = tk.Button(grid1_0, text="Save as CSV", command=save_csv)
            save_csv_button.grid(row = 4, column = 1)   

            # Curve Fitting Attempt...
            fitting_options = ["Linear", "Polynomial (2)"]
            fitting_var = tk.StringVar(tab4_frame)
            fitting_var.set(fitting_options[0])  # Set the default option

            fitting_label = tk.Label(grid1_1, text="Select Fitting Type:")
            fitting_label.grid(row = 0, column = 0)  

            fitting_dropdown = tk.OptionMenu(grid1_1, fitting_var, *fitting_options)
            fitting_dropdown.grid(row = 1, column = 0)  

            # Function to put our curvefit onto the plot
            def fit_and_plot_data():
                x_property = x_var.get()
                y_property = y_var.get()
                x_data = np.array([self.property_dict[i][x_property] for i in range(len(self.property_dict))])
                y_data = np.array([self.property_dict[i][y_property] for i in range(len(self.property_dict))])

                ax.clear()
                ax.scatter(x_data, y_data, label="Data")
                
                # Perform curve fitting based on the selected fitting type
                if fitting_var.get() == "Linear":
                    fit_params = np.polyfit(x_data, y_data, 1)
                    new_x_data = np.linspace(np.min(x_data) * 0.9, np.max(x_data) * 1.1, 1000)
                    fitted_curve = fit_params[0] * new_x_data + fit_params[1]
                    r_squared = 1 - (sum((y_data - (fit_params[0] * x_data + fit_params[1])) ** 2) / ((len(y_data) - 1) * np.var(y_data, ddof=1)))
                    ax.plot(new_x_data, fitted_curve, label=f"Linear Fit (R^2 = {r_squared:.4f})", color="red")
                    
                elif fitting_var.get() == "Polynomial (2)":
                    fit_params = np.polyfit(x_data, y_data, 2)
                    new_x_data = np.linspace(np.min(x_data) * 0.9, np.max(x_data) * 1.1, 1000)
                    fitted_curve = np.polyval(fit_params, new_x_data)
                    r_squared = 1 - (sum((y_data - np.polyval(fit_params, x_data)) ** 2) / ((len(y_data) - 1) * np.var(y_data, ddof=1)))
                    ax.plot(new_x_data, fitted_curve, label=f"2nd order Polynomial Fit (R^2 = {r_squared:.4f})", color="green")
                
                ax.set_xlabel(x_property)
                ax.set_ylabel(y_property)
                ax.set_title(f'Scatter Plot: {x_property} vs {y_property}')
                ax.legend()
                plt.tight_layout()
                update_axis_limits
                canvas.draw()

            fit_button = tk.Button(grid1_1, text="Fit Data", command=fit_and_plot_data)
            fit_button.grid(row = 2, column = 0) 

            #fit_and_plot_data()
                    
    ## TAB 5 ##
    def create_tab5_content(self):
        tab5_outer_frame = tk.Frame(self.tab5)
        tab5_outer_frame.grid(row=0, column=0, sticky="nsew")

        self.canvas_conformer = tk.Canvas(tab5_outer_frame, width=500, height=300)
        self.canvas_conformer.grid(row=0, column=0, padx=5, pady=5)
        
        tab5_inner_frame = tk.Frame(tab5_outer_frame)
        tab5_inner_frame.grid(row=0, column=1, sticky="n")

        settings_frame = tk.LabelFrame(tab5_inner_frame, text="Conformer Selection and Options:")
        settings_frame.grid(row=0, column=0, padx=5, pady=5)

        conformer_label = tk.Label(settings_frame, text="Select Conformer:")
        conformer_label.grid(row=0, column=0, padx=5, pady=5)

        self.conformer_slider = tk.Scale(settings_frame, from_=0, to=len(self.angle), orient="horizontal", command=self.display_selected_conformer)
        self.conformer_slider.grid(row=0, column=1, padx=5, pady=5)
        self.conformer_slider.set(0)  # Set an initial conformer
        
        save_conf_buton = tk.Button(settings_frame, text="Save conformers as .sdf", command=self.write_molecule)
        save_conf_buton.grid(row=1, column=0, columnspan=2, padx=5, pady=5)
        
        draw_pymol_button = tk.Button(settings_frame, text="Draw in PyMol (slow!)", command=self.draw_in_pymol)
        draw_pymol_button.grid(row=2, column=0, columnspan=2, padx=5, pady=5)

        self.pymol_ray_trace_var = BooleanVar()
        self.pymol_ray_trace_var.set(False)
        self.pymol_draw_grid_var = BooleanVar()
        self.pymol_draw_grid_var.set(True)

        pymol_options_label = tk.Label(settings_frame, text="PyMol Options:")
        pymol_options_label.grid(row=3, column=0, columnspan=2, padx=5, pady=5)

        pymol_ray_trace_check = tk.Checkbutton(settings_frame, text="PyMol RayTrace?", variable=self.pymol_ray_trace_var)
        pymol_ray_trace_check.grid(row=4, column=0, columnspan=2, padx=5, pady=5)

        pymol_draw_grid_check = tk.Checkbutton(settings_frame, text="Draw Grid Image?", variable=self.pymol_draw_grid_var)
        pymol_draw_grid_check.grid(row=5, column=0, columnspan=2, padx=5, pady=5)

        pymol_use_uol_smp_check = tk.Checkbutton(settings_frame, text="Use UoL-SMP style?", variable=self.pymol_use_uol_smp_var)
        pymol_use_uol_smp_check.grid(row=6, column=0, columnspan=2, padx=5, pady=5)

        stats_frame = tk.LabelFrame(tab5_inner_frame, text="Conformer Statistics and Properties:")
        stats_frame.grid(row=1, column=0, padx=5, pady=5)

        self.conformer_info_label = tk.Label(stats_frame, text="")
        self.conformer_info_label.grid(row=0, column=0, padx=5, pady=5)
    
    def write_molecule(self):
        '''
        Function for saving conformers as SDF
        '''
        root = tk.Tk()
        root.withdraw()  # Hide the root window

        filetypes = (
            ("SDF", "*.sdf"),
        )

        file_path = filedialog.asksaveasfilename(filetypes=filetypes)
        mol = self.mol_conf # 
            
        # write conformations
        writer = Chem.SDWriter(file_path + '.sdf')
        for cid in range(mol.GetNumConformers()):
            writer.write(mol, confId=cid)

    def draw_in_pymol(self):
        draw.pymol_draw(self.loaded_mol,path=self.pymol_path, ray=self.pymol_ray_trace_var.get(), grid=self.pymol_draw_grid_var.get(), style=self.pymol_use_uol_smp_var.get())
        
    def generate_and_display_first_conformer(self):
        if self.angle != [] and self.loaded_mol != None:
            self.display_selected_conformer(0)
        if self.angle == [] and self.loaded_mol != None:
            print('No conformer/angle data')
            
    def display_selected_conformer(self, _):
        conformer_index = self.conformer_slider.get()
        if 0 <= conformer_index < len(self.angle):
            conformer_angle = f"{self.angle[conformer_index]:.2f}°"
            conformer_energy = f"{self.delta_energy[conformer_index]:.2f} kJ/mol"
            conformer_probability = f"{self.probability[conformer_index]:.2f}"  # Probability is unitless

            # Update the label with conformer information
            info_text = f"Conformer Angle: {conformer_angle}\nConformer ΔE: {conformer_energy}\nConformer Probability: {conformer_probability}"
            self.conformer_info_label.config(text=info_text)

            # Draw and display the selected conformer image
            self.draw_and_display_conformer_image(Chem.RemoveHs(self.loaded_mol), conformer_index)
        else:
            self.conformer_info_label.config(text="Invalid conformer index")

    def draw_and_display_conformer_image(self, mol, conformer_index):
    
        self.canvas_conformer.delete("all") # Clear previous conformer image
        img = self.draw_molecule_image(mol, conformer_index)  # Draw the selected conformer image using RDKit drawing methods
        self.display_image_on_canvas(img) # Display the conformer image on the canvas

    def draw_molecule_image(self, mol, conformer_index):
        if mol is not None:
            # Create an image of the molecule with the selected conformer
            mol_with_conformer = Chem.Mol(mol)
            mol_with_conformer.RemoveAllConformers()
            conformer = mol.GetConformer(conformer_index)
            mol_with_conformer.AddConformer(conformer, assignId=True)

            img = Draw.MolToImage(mol_with_conformer, size=(400, 400))
            img = ImageTk.PhotoImage(Image.fromarray(np.array(img)))

            return img

    def display_image_on_canvas(self, img):
        self.canvas_conformer.create_image(0, 0, anchor=tk.NW, image=img)
        self.canvas_conformer.image = img

    ### settings tab ###
    def create_settings_tab_content(self, entry_gaussian, entry_pymol):

        label_gaussian = tk.Label(self.settings_tab, text="Path to Gaussian:")
        entry_gaussian = tk.Entry(self.settings_tab, width=40)
        button_browse_gaussian = tk.Button(self.settings_tab, text="Browse",
                                            command=lambda: self.browse_gaussian(entry_gaussian))
        
        label_pymol = tk.Label(self.settings_tab, text="Path to PyMOL:")
        entry_pymol = tk.Entry(self.settings_tab, width=40)
        button_browse_pymol = tk.Button(self.settings_tab, text="Browse",
                                        command=lambda: self.browse_pymol(entry_pymol))
        
        label_gaussian.pack(pady=5)
        entry_gaussian.pack(pady=5)
        button_browse_gaussian.pack()
        
        label_pymol.pack(pady=5)
        entry_pymol.pack(pady=5)
        button_browse_pymol.pack()
        
        button_save = tk.Button(self.settings_tab, text="Save Settings",
                                command=lambda: self.save_settings(entry_gaussian.get(), entry_pymol.get()))
        button_save.pack(pady=5)
        
        button_load = tk.Button(self.settings_tab, text="Load Settings",
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
