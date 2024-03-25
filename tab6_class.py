import tkinter as tk
from tkinter import ttk, filedialog
import json

class TabSix(ttk.Frame):
    def __init__(self, parent, shared_data, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.shared_data = shared_data 
        
        ttk.Frame.__init__(self,parent, *args, **kwargs)
        self.create_tab6_content()

    def create_tab6_content(self):

        label_gaussian = ttk.Label(self, text="Path to Gaussian:")
        entry_gaussian = ttk.Entry(self, width=40)
        button_browse_gaussian = ttk.Button(self, text="Browse",
                                            command=lambda: self.browse_gaussian(entry_gaussian))
        
        label_pymol = ttk.Label(self, text="Path to PyMOL:")
        entry_pymol = ttk.Entry(self, width=40)
        button_browse_pymol = ttk.Button(self, text="Browse",
                                        command=lambda: self.browse_pymol(entry_pymol))
        
        label_gaussian.pack(pady=5)
        entry_gaussian.pack(pady=5)
        button_browse_gaussian.pack()
        
        label_pymol.pack(pady=5)
        entry_pymol.pack(pady=5)
        button_browse_pymol.pack()
        
        button_save = ttk.Button(self, text="Save Settings",
                                command=lambda: self.save_settings(entry_gaussian.get(), entry_pymol.get()))
        button_save.pack(pady=5)
        
        button_load = ttk.Button(self, text="Load Settings",
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
        entry_gaussian.delete(0, ttk.END)
        entry_gaussian.insert(0, gaussian_path)
        entry_pymol.delete(0, ttk.END)
        entry_pymol.insert(0, pymol_path)
        
        