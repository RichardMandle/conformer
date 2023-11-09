import tkinter as tk
from tkinter import ttk, filedialog, OptionMenu, BooleanVar, simpledialog  

#import our own modules
import draw
import conformer

#import tabs
from tab1_class import TabOne
from tab2_class import TabTwo
from tab3_class import TabThree
from tab4_class import TabFour
from tab5_class import TabFive
from tab6_class import TabSix

class CallbackHandler:
    '''
    CallbackHandler class is used for updating data between tabs on change; fixes *annoying* behaviour between tabs 2,4,5. 
    There is probably an easier way to fix this but this seems to work OK
    '''
    def __init__(self):
        self.callbacks = {}
        
    def register_callback(self, event_name, callback):
        if event_name not in self.callbacks:
            self.callbacks[event_name] = []
        self.callbacks[event_name].append(callback)

    def call_callbacks(self, event_name, *args, **kwargs):
        for callback in self.callbacks.get(event_name, []):
            callback(*args, **kwargs)
            
class SharedData:
    '''
    This class holds all shared data we'll need in tabs, and also enables easy saving/reloading of conformer sessions
    '''
    def __init__(self):       
        self.loaded_mol = None 
        self.property_dict = None
        self.drawer = None
        self.highlights = []
        self.settings = {}
        self.probability = []
        self.angle = []
        self.energy = []
        self.delta_energy = []
        self.gaussian_path = []
        self.pymol_path = []
        self.mol_conf = None
                
        #Define BooleanVar variables for pymol
        self.pymol_ray_trace_var = tk.BooleanVar()
        self.pymol_ray_trace_var.set(False)
        self.pymol_draw_grid_var = tk.BooleanVar()
        self.pymol_draw_grid_var.set(True)
        self.pymol_use_uol_smp_var = tk.BooleanVar()
        self.pymol_use_uol_smp_var.set(True)
        

class ConformerGUI:
    '''
    Our main class for the GUI, basically just loads other tabs
    '''
    def __init__(self, root):
        self.shared_data = SharedData() # access the SharedData class
        self.callback_handler = CallbackHandler() 
        
        self.root = root
        self.root.title("Conformer GUI")
        self.notebook = ttk.Notebook(self.root)

        self.tab1 = TabOne(self.notebook, self.shared_data, self.callback_handler)
        self.tab2 = TabTwo(self.notebook, self.shared_data, self.callback_handler)
        self.tab3 = TabThree(self.notebook, self.shared_data, self.callback_handler)
        self.tab4 = TabFour(self.notebook, self.shared_data, self.callback_handler)
        self.tab5 = TabFive(self.notebook, self.shared_data, self.callback_handler)
        self.tab6 = TabSix(self.notebook, self.shared_data)

        self.notebook.add(self.tab1, text="Load Molecule")
        self.notebook.add(self.tab2, text="Conformer Search")
        self.notebook.add(self.tab3, text="Analysis #1")
        self.notebook.add(self.tab4, text="Analysis #2")
        self.notebook.add(self.tab5, text="View Conformers")
        self.notebook.add(self.tab6, text="Settings")
        
        self.notebook.bind("<<NotebookTabChanged>>", self.on_tab_change)
        
        self.notebook.pack(expand=True, fill="both")

        entry_gaussian = tk.Entry(self.tab6, width=40)
        entry_pymol = tk.Entry(self.tab6, width=40)

        tab_six_instance = self.tab6
        self.shared_data.gaussian_path, self.shared_data.pymol_path = gaussian_path, pymol_path = tab_six_instance.load_settings()  # Call load_settings on the TabSix instance  # Pass the entries to the method
        
        if entry_gaussian:
            print('Gaussian Path is: ' + self.shared_data.gaussian_path)
        if entry_pymol:
            print('PyMol Path is: ' + self.shared_data.pymol_path)

    
    def on_tab_change(self, event):
        if self.notebook.index(self.notebook.select()) == 2:  # Tab 3
            self.tab3.generate_histograms_and_plots()  # Call the function to update canvas
        
        if self.notebook.index(self.notebook.select()) == 3:  # Tab 4
            if self.shared_data.loaded_mol is not None:
                self.tab3.generate_histograms_and_plots() # run this so we definitely have some energy/angle data
                for widget in self.tab4.winfo_children():
                    widget.destroy()
                self.tab4.create_tab4_content()
            else:
                print("No mol object found (class warning)")       
                
        if self.notebook.index(self.notebook.select()) == 4:  # Tab 5
            self.tab3.generate_histograms_and_plots() # run this so we definitely have some energy/angle data
            self.tab5.generate_and_display_first_conformer()  # Call the function to update the conformer canvas
        
if __name__ == "__main__":
    root = tk.Tk()
    app = ConformerGUI(root)
    root.mainloop()
