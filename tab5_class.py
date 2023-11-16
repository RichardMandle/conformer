import tkinter as tk
from tkinter import ttk, filedialog, BooleanVar
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageTk
import numpy as np
import draw

class TabFive(ttk.Frame):
    def __init__(self, parent, shared_data,  callback_handler, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.shared_data = shared_data 
        self.callback_handler = callback_handler
        ttk.Frame.__init__(self,parent, *args, **kwargs)
        self.create_tab5_content()

    def create_tab5_content(self):
        tab5_outer_frame = ttk.Frame(self)
        tab5_outer_frame.grid(row=0, column=0, sticky="nsew")
        
        self.canvas_conformer = tk.Canvas(tab5_outer_frame, width=550, height=360)
        self.canvas_conformer.grid(row=0, column=0, padx=5, pady=5)
        
        tab5_inner_frame = ttk.Frame(tab5_outer_frame)
        tab5_inner_frame.grid(row=0, column=1, sticky="n")

        settings_frame = ttk.LabelFrame(tab5_inner_frame, text="Conformer Selection and Options:")
        settings_frame.grid(row=0, column=0, padx=5, pady=5)

        conformer_label = ttk.Label(settings_frame, text="Select Conformer:")
        conformer_label.grid(row=0, column=0, padx=5, pady=5)

        self.conformer_slider_value_label = ttk.Label(settings_frame, text="Conformer: 0")
        self.conformer_slider_value_label.grid(row=1, column=0, padx=5, pady=5)

        self.conformer_slider = ttk.Scale(settings_frame, from_=0, to=int(len(self.shared_data.angle)-1), orient="horizontal", command=self.update_conformer_display)
        self.conformer_slider.grid(row=1, column=1, padx=5, pady=5)
        self.conformer_slider.set(0)
    
        save_conf_buton = ttk.Button(settings_frame, text="Save conformers as .sdf", command=self.write_molecule)
        save_conf_buton.grid(row=2, column=0, columnspan=2, padx=5, pady=5)
        
        draw_pymol_button = ttk.Button(settings_frame, text="Draw in PyMol (slow!)", command=self.draw_in_pymol)
        draw_pymol_button.grid(row=3, column=0, columnspan=2, padx=5, pady=5)

        self.pymol_ray_trace_var = BooleanVar()
        self.pymol_ray_trace_var.set(False)
        self.pymol_draw_grid_var = BooleanVar()
        self.pymol_draw_grid_var.set(True)

        pymol_options_label = ttk.Label(settings_frame, text="PyMol Options:")
        pymol_options_label.grid(row=4, column=0, columnspan=2, padx=5, pady=5)

        pymol_ray_trace_check = ttk.Checkbutton(settings_frame, text="PyMol RayTrace?", variable=self.shared_data.pymol_ray_trace_var)
        pymol_ray_trace_check.grid(row=5, column=0, columnspan=2, padx=5, pady=5)

        pymol_draw_grid_check = ttk.Checkbutton(settings_frame, text="Draw Grid Image?", variable=self.shared_data.pymol_draw_grid_var)
        pymol_draw_grid_check.grid(row=6, column=0, columnspan=2, padx=5, pady=5)

        pymol_use_uol_smp_check = ttk.Checkbutton(settings_frame, text="Use UoL-SMP style?", variable=self.shared_data.pymol_use_uol_smp_var)
        pymol_use_uol_smp_check.grid(row=7, column=0, columnspan=2, padx=5, pady=5)

        stats_frame = ttk.LabelFrame(tab5_inner_frame, text="Conformer Statistics and Properties:")
        stats_frame.grid(row=1, column=0, padx=5, pady=5)

        self.conformer_info_label = ttk.Label(stats_frame, text="")
        self.conformer_info_label.grid(row=0, column=0, padx=5, pady=5)
        
        self.callback_handler.register_callback("data_updated", self.update_conformer_slider_range)
    
    def update_conformer_slider_range(self):
        if self.shared_data.angle:
            max_index = len(self.shared_data.angle) - 1
        else:
            max_index = 0
        self.conformer_slider.config(to=max_index)
        self.conformer_slider.set(0)
        self.update_conformer_display(0)
        
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
        mol = self.shared_data.mol_conf 

        writer = Chem.SDWriter(file_path + '.sdf')
        for cid in range(mol.GetNumConformers()):
            writer.write(mol, confId=cid)

    def draw_in_pymol(self):
        draw.pymol_draw(self.shared_data.mol_conf,path=self.pymol_path, ray=self.pymol_ray_trace_var.get(), grid=self.pymol_draw_grid_var.get(), style=self.pymol_use_uol_smp_var.get())
        
    def generate_and_display_first_conformer(self):
        if self.shared_data.angle != [] and self.shared_data.mol_conf != None:
            self.display_selected_conformer(0)
        if self.shared_data.angle == [] and self.shared_data.mol_conf != None:
            print('No conformer/angle data')
            
    def update_conformer_display(self, val):
        if int(len(self.shared_data.angle)) >=1: #check its got something in it.
            self.conformer_slider_value_label.config(text=f"Conformer: {int(float(val))}")
            self.display_selected_conformer(val)
    
    def display_selected_conformer(self, _):
        if not self.shared_data.angle:  # Check if angle data is empty
            self.conformer_info_label.config(text="No conformers loaded")
            return
            
        conformer_index = int(self.conformer_slider.get())
        if 0 <= conformer_index < len(self.shared_data.angle):
            conformer_angle = f"{self.shared_data.angle[conformer_index]:.2f}°"
            conformer_energy = f"{self.shared_data.delta_energy[conformer_index]:.2f} kJ/mol"
            conformer_probability = f"{self.shared_data.probability[conformer_index]:.2f}"  # Probability is unitless

            info_text = f"Conformer Angle: {conformer_angle}\nConformer ΔE: {conformer_energy}\nConformer Probability: {conformer_probability}"
            self.conformer_info_label.config(text=info_text)

            self.draw_and_display_conformer_image(Chem.RemoveHs(self.shared_data.mol_conf), conformer_index)
        else:
            self.conformer_info_label.config(text="Invalid conformer index")

    def draw_and_display_conformer_image(self, mol, conformer_index):
    
        self.canvas_conformer.delete("all") # Clear previous conformer image
        img = self.draw_molecule_image(mol, conformer_index)  # Draw the selected conformer image using RDKit drawing methods
        self.display_image_on_canvas(img) # Display the conformer image on the canvas

    def draw_molecule_image(self, mol, conformer_index):
        if mol is not None:
            mol_with_conformer = Chem.Mol(mol)
            mol_with_conformer.RemoveAllConformers()
            conformer = mol.GetConformer(conformer_index)
            mol_with_conformer.AddConformer(conformer, assignId=True)

            img = Draw.MolToImage(mol_with_conformer, size=(550, 360))
            img = ImageTk.PhotoImage(Image.fromarray(np.array(img)))

            return img

    def display_image_on_canvas(self, img):
        self.canvas_conformer.create_image(0, 0, anchor=tk.NW, image=img)
        self.canvas_conformer.image = img