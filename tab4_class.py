import tkinter as tk
from tkinter import ttk, filedialog
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import conformer

class TabFour(ttk.Frame):
    def __init__(self, parent, shared_data,  callback_handler, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.shared_data = shared_data 
        self.callback_handler = callback_handler
        self.current_figure = None # use to keep track of figures and close each time we replot to not waste memory
        self.create_tab4_content()

    def create_tab4_content(self):
        
        self.fig, self.ax = plt.subplots(figsize=(600/plt.rcParams['figure.dpi'], 400/plt.rcParams['figure.dpi']))    
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().grid(row=0, column=0, rowspan=2, sticky="nsew")
                
        toolbar_frame = tk.Frame(self)
        toolbar_frame.grid(row=3, column=0, sticky="ew")
        toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        toolbar.update()
        toolbar.pack(side=tk.TOP, fill=tk.X)
        
        grid1 = ttk.LabelFrame(self, text="Plots")
        grid1.grid(row=0, column=1, padx=10, pady=10, sticky="n")

        grid1_0 = ttk.LabelFrame(grid1, text="Plot Settings and Options")
        grid1_0.grid(row=0, column=0, padx=10, pady=10, sticky="w")

        grid1_1 = ttk.LabelFrame(grid1, text="Fit Settings and Options")
        grid1_1.grid(row=1, column=0, padx=10, pady=10, sticky="w")

        self.x_var = tk.StringVar(self)
        self.y_var = tk.StringVar(self)
        
        save_csv_button = tk.Button(grid1_0, text="Save as CSV", command=self.save_csv)
        save_csv_button.grid(row = 5, column = 0)

        save_png_button = tk.Button(grid1_0, text="Save as PNG", command=self.save_png)
        save_png_button.grid(row = 6, column = 0)               

        # Curve Fitting Attempt...
        fitting_options = ["Linear", "Polynomial (2)", "Polynomial (3)", "Polynomial (4)"]
        self.fitting_var = tk.StringVar(self)
        self.fitting_var.set(fitting_options[0])  # Set the default option

        fitting_label = tk.Label(grid1_1, text="Select Fitting Type:")
        fitting_label.grid(row = 0, column = 0)  

        fitting_dropdown = tk.OptionMenu(grid1_1, self.fitting_var, *fitting_options)
        fitting_dropdown.grid(row = 1, column = 0)  

        fit_button = tk.Button(grid1_1, text="Fit Data", command=self.fit_and_plot_data)
        fit_button.grid(row = 2, column = 0) 

        if self.shared_data.mol_conf is None:
            print("No conformer information found")
            no_mol_label = tk.Label(self, text="No mol-conformer object found. \nPlease load a molecule and generate conformers in a previous tab.")
            no_mol_label.grid(row=0, column=0, padx=10, pady=10, sticky="n")
            return
        
        else:
            available_properties = list(self.shared_data.property_dict[0].keys())

            x_label = tk.Label(grid1_0, text="Select X Property:")
            x_label.grid(row = 0, column = 0)

            x_dropdown = tk.OptionMenu(grid1_0, self.x_var, *available_properties)
            x_dropdown.grid(row = 1, column = 0)

            y_label = tk.Label(grid1_0, text="Select Y Property:")
            y_label.grid(row = 2, column = 0)

            y_dropdown = tk.OptionMenu(grid1_0, self.y_var, *available_properties)
            y_dropdown.grid(row = 3, column = 0)

            plot_button = tk.Button(grid1_0, text="Plot Data", command=lambda: self.plot_data())
            plot_button.grid(row = 4, column = 0)
        
            self.x_var.set(available_properties[-1])
            self.y_var.set(available_properties[0])
    
    def get_xy_properties(self):
        return  self.x_var.get(),self.y_var.get()
    
    def get_xy_data(self):
        x_property, y_property = self.get_xy_properties()
        x_data = [self.shared_data.property_dict[i][x_property] for i in range(len(self.shared_data.property_dict))]
        y_data = [self.shared_data.property_dict[i][y_property] for i in range(len(self.shared_data.property_dict))]
        return x_data, y_data

    def plot_data(self):
        self.ax.clear()
        x_property, y_property = self.get_xy_properties()
        x_data, y_data = self.get_xy_data()
        
        self.ax.scatter(x_data, y_data, label="Data")

        self.ax.set_xlabel(x_property)
        self.ax.set_ylabel(y_property)
        self.ax.set_title(f'Scatter Plot: {x_property} vs {y_property}')
        plt.tight_layout()
        plt.legend()
        self.canvas.draw()
        self.canvas.get_tk_widget().update() 
        
    def save_png(self):
        file_path = filedialog.asksaveasfilename(defaultextension='.png', filetypes=[('PNG Image', '*.png')])
        if file_path:
            self.fig.savefig(file_path)
        
    def save_csv(self):
        file_path = filedialog.asksaveasfilename(defaultextension='.csv', filetypes=[('CSV Files', '*.csv')])
        if file_path:
            x_property, y_property = self.get_xy_properties()
            x_data, y_data = self.get_xy_data()
            
            data = np.column_stack((x_data, y_data))
            np.savetxt(file_path, data, delimiter=',', header=f'{x_property},{y_property}', comments='')
            
    def fit_and_plot_data(self):
        if self.current_figure is not None:
            plt.close()
            
        self.plot_data() # call the plot the x/y data
        x_property, y_property = self.get_xy_properties() # get the data here so we can fit it
        x_data, y_data = self.get_xy_data()
        
        fitting_choice = self.fitting_var.get()
        degree_mapping = {"Linear": 1, "Polynomial (2)": 2, "Polynomial (3)": 3, "Polynomial (4)": 4}
        
        if fitting_choice in degree_mapping:
            self.perform_curve_fit(x_data, y_data, degree_mapping[fitting_choice])
        plt.legend()
        self.canvas.draw()

        
    def perform_curve_fit(self, x_data, y_data, degree):
        # this code computes the polynomial fit and plots it
        fit_params = np.polyfit(x_data, y_data, degree)
        new_x_data = np.linspace(np.min(x_data) * 0.9, np.max(x_data) * 1.1, 1000)
        fitted_curve = np.polyval(fit_params, new_x_data)
        r_squared = 1 - (sum((y_data - np.polyval(fit_params, x_data)) ** 2) / ((len(y_data) - 1) * np.var(y_data, ddof=1)))
        label = f"{degree}th order Polynomial Fit (R^2 = {r_squared:.4f})"
        if degree == 1:
            label = f"Linear Fit (R^2 = {r_squared:.4f})"
        self.ax.plot(new_x_data, fitted_curve, label=label, color="red" if degree == 1 else "green")
    