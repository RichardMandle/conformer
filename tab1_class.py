import tkinter as tk
from tkinter import filedialog,  ttk  
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageTk
import numpy as np
from collections import defaultdict
import io
import pickle # for saving/loading state

class TabOne(ttk.Frame):

    def __init__(self, parent, shared_data, callback_handler, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.shared_data = shared_data 
        self.callback_handler=callback_handler
        ttk.Frame.__init__(self,parent, *args, **kwargs)
        self.create_tab1_content()
        
    def create_ui_element(self, element_type, row, column, text='', command=None):
        element = element_type(self, text=text, command=command if command else None)
        element.grid(row=row, column=column, padx=5, pady=5, sticky="w")
        return element
         
    def create_tab1_content(self):

        labels = ["Load from File:", "Generate from SMILES/SMARTS"]
        for i, label_text in enumerate(labels):
            ttk.Label(self, text=label_text).grid(row=0, column=i+1, padx=0, pady=10, sticky="w")

        self.create_ui_element(ttk.Button, 1, 1, "Open Molecule from File", self.open_file_dialog)
        self.molecule_entry = ttk.Entry(self)
        self.molecule_entry.grid(row=1, column=2, padx=5, pady=5, sticky="w")

        self.create_ui_element(ttk.Button, 2, 1, "Save Conformer Session", self.save_session)
        self.create_ui_element(ttk.Button, 3, 1, "Load Conformer Session", self.load_session)
        self.create_ui_element(ttk.Button, 2, 2, "Generate Molecule", self.get_mol_from_text)

        self.canvas = tk.Canvas(self, width=775, height=400)
        self.canvas.grid(row=4, column=1, columnspan=4, padx=10, pady=10, sticky='ew')
        self.canvas.bind("<Button-1>", self.get_nearest_atom)
    
    def save_session(self):
        file_path = filedialog.asksaveasfilename(
                    defaultextension='.pkl',
                    filetypes=[('Pickle files', '*.pkl'), ('All files', '*.*')],
                    title="Save File")
                    
        data_to_save = {'loaded_mol' : self.shared_data.loaded_mol,
                        'property_dict' : self.shared_data.property_dict,
                        'highlights' : self.shared_data.highlights,
                        'settings' : self.shared_data.settings,
                        'probability' : self.shared_data.probability,
                        'angle' : self.shared_data.angle,
                        'energy' : self.shared_data.energy,
                        'delta_energy' : self.shared_data.delta_energy,
                        'gaussian_path' : self.shared_data.gaussian_path,
                        'pymol_path' : self.shared_data.pymol_path,
                        'mol_conf' : self.shared_data.mol_conf}
        try:
            with open(file_path, 'wb') as file:
                pickle.dump(data_to_save, file)

            print(f"Saved data has {self.shared_data.mol_conf.GetNumConformers()} conformers")
            print("Session saved successfully!")
            
        except Exception as e:
            print(f"An error occurred while saving: {e}")
            
    def load_session(self):
        file_path = filedialog.askopenfilename(
            defaultextension='.pkl',
            filetypes=[('Pickle files', '*.pkl'), ('All files', '*.*')],
            title="Open File")

        if file_path:
            try:
                with open(file_path, 'rb') as file:
                    data_loaded = pickle.load(file)
                print("Session loaded successfully!")
                              
                # Update the attributes with the loaded data
                self.shared_data.loaded_mol = data_loaded.get('loaded_mol', None)
                self.shared_data.property_dict = data_loaded.get('property_dict', {})
                self.shared_data.highlights = data_loaded.get('highlights', None)
                self.shared_data.settings = data_loaded.get('settings', {})
                self.shared_data.probability = data_loaded.get('probability', None)
                self.shared_data.angle = data_loaded.get('angle', None)
                self.shared_data.energy = data_loaded.get('energy', None)
                self.shared_data.delta_energy = data_loaded.get('delta_energy', None)
                self.shared_data.gaussian_path = data_loaded.get('gaussian_path', None)
                self.shared_data.pymol_path = data_loaded.get('pymol_path', None)
                self.shared_data.mol_conf = data_loaded.get('mol_conf', None)

                print(f"Loaded data has {data_loaded.get('mol_conf', None).GetNumConformers()} conformers")

            except Exception as e:
                print(f"An error occurred while loading data: {e}")    
                
        self.draw_molecule(clear_confs=False) # redraw the molecule from the conformer given...
        self.callback_handler.call_callbacks("data_updated") # notify tabs to update via callback
        
    def get_mol_from_text(self):
        text = self.molecule_entry.get()
        if text:
            mol = Chem.MolFromSmiles(text) # try SMILES first
            if mol:
                print('Loaded molecule from SMILES input')
            else:
                mol = Chem.MolFromSmarts(text) # try SMILES next
                if mol:
                    print('Loaded molecule from SMARTS input')
                    
        if mol is not None:
            self.shared_data.loaded_mol = Chem.AddHs(mol)  # Store the loaded mol object and add hydrogens
            self.draw_molecule(clear_confs=True)  
            
    def draw_molecule(self, clear_confs=True):
        if self.shared_data.loaded_mol is not None:
            # set clear_confs to delete conformer data; useful when loading a new molecule, but not when loading a prev. session
            AllChem.Compute2DCoords(self.shared_data.loaded_mol,clearConfs=clear_confs) # clearConfs = False to avoid overwriting conformer data we already have for this mol object.
            
            self.shared_data.drawer = rdMolDraw2D.MolDraw2DCairo(500, 300)
            
            if len(self.shared_data.highlights) == 0:
                self.shared_data.drawer.DrawMolecule(self.shared_data.loaded_mol)
            
            if len(self.shared_data.highlights) > 0:
                colors = [(0.0039, 0.8, 0.7176, 0.3),  # define some colours to use
                        (0.9843, 0.0078, 0.4980, 0.3), 
                        (1.0, 0.9843, 0.0, 0.3), 
                        (0.2235, 0.2353, 1.0, 0.3)]
                        
                athighlights = defaultdict(list)
                arads = {}
                
                for a in self.shared_data.loaded_mol.GetAtoms():
                    if a.GetIdx() in [x -1 for x in self.shared_data.highlights]:
                        aid = a.GetIdx()
                        athighlights[aid].append(colors[[x -1 for x in self.shared_data.highlights].index(a.GetIdx())])
                        arads[aid] = 0.3
                        
                self.shared_data.drawer.DrawMoleculeWithHighlights(self.shared_data.loaded_mol,"",dict(athighlights),{},arads,{})
                
            self.shared_data.drawer.FinishDrawing()
            
            img = self.shared_data.drawer.GetDrawingText()
            img = Image.open(io.BytesIO(img))
            self.photo = ImageTk.PhotoImage(img)
    
            self.canvas.create_image(0, 0, image=self.photo, anchor=tk.NW) # Display the above image on the canvas

    def get_nearest_atom(self, event):
        if self.shared_data.loaded_mol is None:
            return
    
        x, y = event.x, event.y
        min_distance = float(10) # set a cutoff for how far away we can be (in pixels??)
        nearest_atom = None
    
        for i in range(len(self.shared_data.loaded_mol.GetAtoms())):
            drawn_atom_coords = self.shared_data.drawer.GetDrawCoords(i)
            atom_x, atom_y = drawn_atom_coords.x, drawn_atom_coords.y
            distance = np.sqrt((atom_x - x)**2 + (atom_y - y)**2)
            if distance < min_distance:
                min_distance = distance
                nearest_atom = i
    
        if nearest_atom is not None:
            self.shared_data.highlights.append(nearest_atom+1)
            if len(self.shared_data.highlights) == 5: # reset list if it has 5 entries
                self.shared_data.highlights = [nearest_atom+1]
                
            self.draw_molecule(clear_confs=False)
            
    def open_file_dialog(self):
        file_path = filedialog.askopenfilename(filetypes=[("MOL ", "*.mol"), ("MOL2", "*.mol2"), ("SDF", "*.sdf")])
        if file_path:
            try:
                mol = None
                file_extension = file_path.split(".")[-1].lower()
                if file_extension == "mol":
                    mol = Chem.MolFromMolFile(file_path)
                elif file_extension == "mol2":
                    mol = Chem.MolFromMol2File(file_path)
                elif file_extension == "sdf":
                    suppl = Chem.SDMolSupplier(file_path)
                    for m in suppl:
                        mol = m # assume there is just 1 in the sdf... very weak code TO DO - fix
                
            except Exception as e:
                print("Error:", str(e))
                
            if mol is not None:
                self.shared_data.loaded_mol = Chem.AddHs(mol)
                print("Molecule successfully loaded!\n")
                self.draw_molecule(clear_confs=True)   