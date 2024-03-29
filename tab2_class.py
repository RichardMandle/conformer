import tkinter as tk
from tkinter import ttk 
import conformer
import os
import threading
from rdkit import Chem

class TabTwo(ttk.Frame):
    def __init__(self, parent, shared_data,  callback_handler, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.shared_data = shared_data 
        self.callback_handler=callback_handler
        ttk.Frame.__init__(self,parent, *args, **kwargs)
        self.create_tab2_content()

    def create_tab2_content(self):
        tab2_label = ttk.Label(self, text="Conformer Search Options")
        tab2_label.pack(pady=5)

        vector_settings_grid = ttk.LabelFrame(self, text="Angle Definition Settings")
        vector_settings_grid.pack(padx=10, pady=5, fill="both", expand=True)
        
        angle_method_label = ttk.Label(vector_settings_grid, text="Angle Definition Method:")
        angle_method_label.grid(row=0, column=0, padx=10, pady=5, sticky="e")

        self.shared_data.angle_methods = ["None","by atom index", "by SMILES match", "by SMARTS match"]
        self.selected_angle_method = tk.StringVar()
        
        angle_method_dropdown = ttk.Combobox(vector_settings_grid, textvariable=self.selected_angle_method, values=self.shared_data.angle_methods)
        angle_method_dropdown.grid(row=0, column=1, padx=10, pady=5, sticky="w")
        angle_method_dropdown.set(self.shared_data.angle_methods[0])
        
        smiles_group1_label = ttk.Label(vector_settings_grid, text="SMILES Group 1:")
        self.smiles_group1_entry = ttk.Entry(vector_settings_grid, width=30)
        smiles_group2_label = ttk.Label(vector_settings_grid, text="SMILES Group 2:")
        self.smiles_group2_entry = ttk.Entry(vector_settings_grid, width=30)

        smarts_group1_label = ttk.Label(vector_settings_grid, text="SMARTS Group 1:")
        self.smarts_group1_entry = ttk.Entry(vector_settings_grid, width=30)
        smarts_group2_label = ttk.Label(vector_settings_grid, text="SMARTS Group 2:")
        self.smarts_group2_entry = ttk.Entry(vector_settings_grid, width=30)
        
        def update_settings():
            # function to update conformer search settings.
            angle_method = self.selected_angle_method.get()
            for widget in vector_settings_grid.winfo_children():
                widget.grid_remove() # Hide all widgets on update

            angle_method_label.grid(row=0, column=0, padx=10, pady=5, sticky="e")
            angle_method_dropdown.grid(row=0, column=1, padx=10, pady=5, sticky="w")

            # Display appropriate labels and entries based on angle method
            if angle_method == "by atom index":
                if len(self.shared_data.highlights) == 4:
                    # Show highlights information
                    highlights_label = ttk.Label(vector_settings_grid, text=f'Atom group #1: {self.shared_data.highlights[0]}, {self.shared_data.highlights[1]}')
                    highlights_label.grid(row=1, column=0, columnspan=2, padx=10, pady=5)
    
                    highlights_label2 = ttk.Label(vector_settings_grid, text=f'Atom group #2: {self.shared_data.highlights[2]}, {self.shared_data.highlights[3]}')
                    highlights_label2.grid(row=2, column=0, columnspan=2, padx=10, pady=5)
                
                if len(self.shared_data.highlights) != 4:
                    highlights_label = ttk.Label(vector_settings_grid, text='Please select atoms using the tool on tab #1')
                    highlights_label.grid(row=1, column=0, columnspan=2, padx=10, pady=5)
            
            elif angle_method == "by SMILES match" or angle_method == "by SMARTS match":
                group1_label.grid(row=1, column=0, padx=10, pady=5, sticky="e")
                self.group1_entry.grid(row=1, column=1, padx=10, pady=5, sticky="w")

                group2_label.grid(row=2, column=0, padx=10, pady=5, sticky="e")
                self.group2_entry.grid(row=2, column=1, padx=10, pady=5, sticky="w")

            elif angle_method == "None":
                print('No angle method selected; will compute energy only')
                
        self.selected_angle_method.trace("w", lambda *args: update_settings())  # Set up the trace to call update_settings when the angle definition method changes

        update_settings()  # Call update_settings function initially
        
        settings_grid = ttk.LabelFrame(self, text="Conformer Search Settings")
        settings_grid.pack(padx=10, pady=5, fill="both", expand=True)
        settings_labels = [ "Number of Conformers:", "Job Name:"]
        self.shared_data.settings_entries = []
        for i, label_text in enumerate(settings_labels):
            label = ttk.Label(settings_grid, text=label_text)
            label.grid(row=i , column=0, padx=10, pady=5, sticky="e")
            entry = ttk.Entry(settings_grid, width=30)
            entry.grid(row=i , column=1, padx=10, pady=5, sticky="w")
            self.shared_data.settings_entries.append(entry)

        method_label = ttk.Label(settings_grid, text="Conformer Method:")
        method_label.grid(row=len(settings_labels) , column=0, padx=10, pady=5, sticky="e")

        self.methods = conformer.get_conformer_methods()
        self.selected_method = tk.StringVar()

        method_dropdown = ttk.Combobox(settings_grid, textvariable=self.selected_method, values=self.methods)
        method_dropdown.grid(row=len(settings_labels) , column=1, padx=10, pady=5, sticky="w")
        method_dropdown.set(self.methods[0])

        eval_method_label = ttk.Label(settings_grid, text="Evaluation Method:")
        eval_method_label.grid(row=len(settings_labels) + 1 , column=0, padx=10, pady=5, sticky="e")

        self.eval_methods = ["MMFF (internal)", "Gaussian (external)"]
        self.selected_eval_method = tk.StringVar()

        eval_method_dropdown = ttk.Combobox(settings_grid, textvariable=self.selected_eval_method, values=self.eval_methods)
        eval_method_dropdown.grid(row=len(settings_labels) + 1 , column=1, padx=10, pady=5, sticky="w")
        eval_method_dropdown.set(self.eval_methods[0])

        # additional settings widgets for Gaussian options
        self.additional_settings_widgets = []
        options_label = ttk.Label(settings_grid, text="Gaussian Job Options:")
        options_entry = ttk.Entry(settings_grid, width=30)
        cores_label = ttk.Label(settings_grid, text="Number of CPU cores:")
        cores_entry = ttk.Entry(settings_grid, width=30)
        ram_label = ttk.Label(settings_grid, text="RAM / GB:")
        ram_entry = ttk.Entry(settings_grid, width=30)

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

        # Set up the trace to call update_additional_settings when the evaluation method changes
        self.selected_eval_method.trace("w", update_additional_settings)

        search_button = ttk.Button(settings_grid, text="Start Conformer Search", command=self.start_conformer_search)
        search_button.grid(row=len(settings_labels) + 6, columnspan=2, pady=5)

        search_button = ttk.Button(settings_grid, text="Recalculate Bend-Angles", command=self.conformer_analysis_logic)
        search_button.grid(row=len(settings_labels) + 7, columnspan=2, pady=5)
                
        advanced_settings_button = ttk.Button(settings_grid, text="Advanced Settings", command=self.show_advanced_settings)
        advanced_settings_button.grid(row=len(settings_labels) + 8, columnspan=2, pady=5)

        self.status_label = ttk.Label(settings_grid, text="")
        self.status_label.grid(row=len(settings_labels) + 9, columnspan=2, pady=5)
        
        
    def show_advanced_settings(self):
        # this function makes a popup window for some more advanced settings
        advanced_settings_window = tk.Toplevel(self)
        advanced_settings_window.title("Advanced Settings")

        options_frame = ttk.LabelFrame(advanced_settings_window, text="Advanced Options")
        options_frame.grid(row=0, column=0, padx=10, pady=5)

        advanced_options = [
            ("RMS Threshold for Pruning", "0.675", "0.675"),
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
            label = ttk.Label(options_frame, text=(label_text + ':'))
            label.grid(row=row, column=0, padx=10, pady=5, sticky="e")

            option_var = None  # Initialize option_var as None

            if isinstance(values, list):
                option_var = tk.StringVar(value=default_value)
                option_menu = ttk.Combobox(options_frame, textvariable=option_var, values=values)
                option_menu.grid(row=row, column=1, padx=10, pady=5, sticky="w")
            else:
                option_var = tk.StringVar(value=default_value)
                entry = ttk.Entry(options_frame, textvariable=option_var, width=30)
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
                            mmff_label = ttk.Label(mmff_options_frame, text=mmff_label_text)
                            mmff_label.grid(row=mmff_row, column=0, padx=10, pady=5, sticky="e")

                            mmff_entry_var = tk.StringVar(value=mmff_default_value)
                            mmff_entry = ttk.Entry(mmff_options_frame, textvariable=mmff_entry_var, width=30)
                            mmff_entry.grid(row=mmff_row, column=1, padx=10, pady=5, sticky="w")

                            mmff_option_vars[mmff_label_text] = mmff_entry_var

                option_var.trace_add("write", lambda *args: show_mmff_options())

                    
        def save_advanced_settings():
            for label_text, option_var in option_vars.items():
                self.shared_data.settings[label_text] = option_var.get()

            if mmff_options_frame:
                mmff_settings = {}
                for mmff_label_text, mmff_entry_var in mmff_option_vars.items():
                    mmff_settings[mmff_label_text] = mmff_entry_var.get()
                self.shared_data.settings["MMFF Options"] = mmff_settings

        save_button = ttk.Button(advanced_settings_window, text="Save Settings", command=save_advanced_settings)
        save_button.grid(row=len(advanced_options) + 1, column=0, columnspan=2, pady=5)
    
    def update_status_label(self, message):
        # updates the status label
        self.status_label.config(text=message)    
        
    def start_conformer_search(self):
        # start a seperate thread for conformer generation so we don't freeze the GUI
        conformer_thread = threading.Thread(target=self.conformer_generation_logic)
        conformer_thread.start()
        self.update_status_label("Conformer Generation In Progress...\nPlease be patient")
        
    def conformer_generation_logic(self):
        local_mol = Chem.Mol(self.shared_data.loaded_mol) # COPY the mol object so we don't break it with reselecting things later.
        self.settings_dict = conformer.get_default_search_settings() # retrieve the search defaults from conformer.pymol

        if self.selected_eval_method.get()  == 'Gaussian (external)': # if its the gaussian method then do this
            self.settings_dict['gaussian_job'] = True 
            self.settings_dict['gaussian_options'] = self.additional_settings_widgets[1].get()
            self.settings_dict['gaussian_nproc'] = self.additional_settings_widgets[3].get()
            self.settings_dict['gaussian_vmem'] = self.additional_settings_widgets[5].get()
      
        if self.shared_data.settings:
            print('Using Advanced Settings')
            self.settings_dict['rms_threshold'] = float(self.shared_data.settings['RMS Threshold for Pruning'])
            self.settings_dict['use_torsion_pref'] = self.shared_data.settings['Use Torsion Preferences'] == 'True'
            self.settings_dict['use_knowledge'] = self.shared_data.settings['Use Chemical Knowledge'] == 'True'
            self.settings_dict['use_random_coords'] = self.shared_data.settings['Use Random Coords'] == 'True'
            self.settings_dict['MMFF Variant'] = self.shared_data.settings["MMFF Variant"]
            self.settings_dict['random_seed'] = int(self.shared_data.settings['Random Seed'])
            self.settings_dict['opt'] = self.shared_data.settings['Do MMFF Optimization'] == 'True'

            if self.settings_dict['opt'] == True:
                self.settings_dict['max_opt_iter']=int(self.shared_data.settings["MMFF Options"]['Max Iterations'])
                self.settings_dict['min_energy_MMFF']= float(self.shared_data.settings["MMFF Options"]['Energy Threshold (kcal/mol)'])
        
        # do the conformer search, store the conformers in self.shared_data.mol_conf and energies in self.shared_data.energy:
        self.shared_data.mol_conf, self.shared_data.energy = conformer.conf_gen(mol=local_mol,
                                                        embeded_method=self.selected_method.get(),
                                                        num_of_conformer=int(self.shared_data.settings_entries[0].get()),
                                                        name=self.shared_data.settings_entries[1].get(),
                                                        rms_thresh = self.settings_dict['rms_threshold'],
                                                        use_torsion_pref = self.settings_dict['use_torsion_pref'],
                                                        use_knowledge = self.settings_dict['use_knowledge'],
                                                        use_random_coords = self.settings_dict['use_random_coords'],
                                                        random_seed = self.settings_dict['random_seed'],
                                                        MMFF_variant = self.settings_dict['MMFF Variant'],
                                                        opt = self.settings_dict['opt'],
                                                        max_iter= self.settings_dict['max_opt_iter'],
                                                        min_energy_MMFF= self.settings_dict['min_energy_MMFF'])
                
        self.conformer_analysis_logic() # do the analysis
        
    def conformer_analysis_logic(self):
        # this function splits the analysis logic from the conformer generation logic so we can recompute as needed
        if not hasattr(self, 'settings_dict'):  # Check if self has settings_dict
            print('using default settings')
            self.settings_dict = conformer.get_default_search_settings()
        
        if self.shared_data.angle != []:
            self.shared_data.angle = [] # clear existing angle data.
            self.shared_data.probability = [] # clear existing probability data.
        
        atoms=[]
        vec1 = ''
        vec2 = ''
        if self.settings_dict['opt'] == True:
            print('Doing MMFF optimisation')
                            
        if self.selected_angle_method.get() == 'by atom index':
            vector_definition_method = 'atoms'
            atoms = [x-1 for x in self.shared_data.highlights] #remember, rdkit indexes from zero...

        if self.selected_angle_method.get() == 'by SMILES match' or self.selected_angle_method.get() == 'by SMARTS match':
            if self.selected_angle_method.get() == 'by SMARTS match':
                vector_definition_method = 'sma'
            if self.selected_angle_method.get() == 'by SMARTS match':
                vector_definition_method = 'smi'
            vec1 = self.group1_entry.get()
            vec2 = self.group2_entry.get()

        
        if self.selected_angle_method.get() == 'None':
            # in the case that its 'None', just use some dummy atom coordinates.
            # The reason to do this is maybe a Gaussian calculation is requested, and currently
            # this is only allowed when an angle calculation is performed. TO DO - fix this by
            # refactoring the code to permit external Gaussian jobs without angle information. 
            vector_definition_method = 'atoms'
            atoms = [0,1,2,3]

        if self.settings_dict['gaussian_job'] == True:
            self.update_status_label("Running external Gaussian jobs...\nSee terminal for remaining time")
        ###print(self.settings_dict)    
        self.update_status_label("Analysing Conformers...\nPlease be patient")
        self.shared_data.energy,self.shared_data.angle = conformer.conf_analysis(
            mol_conf=self.shared_data.mol_conf,
            energy = self.shared_data.energy,
            vec_def_method=vector_definition_method,
            atom_idx_list=atoms,
            g_path = self.shared_data.gaussian_path,
            vector1 = vec1,
            vector2 = vec2,
            name = self.shared_data.settings_entries[1].get(),
            Gauss = self.settings_dict['gaussian_job'],
            options = self.settings_dict['gaussian_options'],
            cores = self.settings_dict['gaussian_nproc'],
            ram = self.settings_dict['gaussian_vmem'],
            write_only = False # TO DO - this is hard coded if you ever want to expose it
        )
        # to do - when reloading we aren't getting angle data for all conformers...
        if self.shared_data.angle:
            print('Successfully screened ' + str(len(self.shared_data.energy)) + ' conformers!')
            if self.shared_data.angle != []:
                self.shared_data.delta_e , self.shared_data.probability = conformer.calc_bend_prob(self.shared_data.energy) # get initial value at arbitrary temperature
            if len(self.shared_data.angle) != len(self.shared_data.probability):
                print(f"Error: Mismatch in array lengths. Angles: {len(self.shared_data.angle)}, Probabilities: {len(self.shared_data.probability)}")
                #print([a for a in self.shared_data.angle])
                print('*****')
                #print([p for p in self.shared_data.probability])
                return
                
            self.shared_data.property_dict = conformer.get_3d_descriptors(self.shared_data.mol_conf,self.shared_data.angle,self.shared_data.probability) # Get available properties from the property_dicts
            self.callback_handler.call_callbacks("data_updated") # notify tabs to update
            self.update_status_label(f"Conformer Generation Complete!\nGenerated {str(len(self.shared_data.energy))} conformers!")