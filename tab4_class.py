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
        ttk.Frame.__init__(self,parent, *args, **kwargs)
        self.create_tab4_content()

    def create_tab4_content(self):

        fig, ax = plt.subplots(figsize=(600/plt.rcParams['figure.dpi'], 400/plt.rcParams['figure.dpi']))
        canvas = FigureCanvasTkAgg(fig, master=self)
        canvas.get_tk_widget().grid(row=0, column=0, rowspan=2, sticky="nsew")
                
        toolbar_frame = tk.Frame(self)
        toolbar_frame.grid(row=3, column=0, sticky="ew")
        toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
        toolbar.update()
        toolbar.pack(side=tk.TOP, fill=tk.X)
        
        grid1 = ttk.LabelFrame(self, text="Plots")
        grid1.grid(row=0, column=1, padx=10, pady=10, sticky="n")

        grid1_0 = ttk.LabelFrame(grid1, text="Plot Settings and Options")
        grid1_0.grid(row=0, column=0, padx=10, pady=10, sticky="w")

        grid1_1 = ttk.LabelFrame(grid1, text="Fit Settings and Options")
        grid1_1.grid(row=1, column=0, padx=10, pady=10, sticky="w")

        x_var = tk.StringVar(self)
        y_var = tk.StringVar(self)

        if self.shared_data.mol_conf is None:
            print("No conformer information found")
            no_mol_label = tk.Label(self, text="No mol-conformer object found. \nPlease load a molecule and generate conformers in a previous tab.")
            no_mol_label.grid(row=0, column=0, padx=10, pady=10, sticky="n")
            return
        
        else:
            self.shared_data.property_dict = conformer.get_3d_descriptors(self.shared_data.mol_conf,self.shared_data.angle,self.shared_data.probability) # Get available properties from the property_dict
            available_properties = list(self.shared_data.property_dict[0].keys())

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

                x_data = [self.shared_data.property_dict[i][x_property] for i in range(len(self.shared_data.property_dict))]
                y_data = [self.shared_data.property_dict[i][y_property] for i in range(len(self.shared_data.property_dict))]

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
                try:
                    x_limits = x_axis_limit_entry.get().split(',')
                    y_limits = y_axis_limit_entry.get().split(',')
                    
                    if len(x_limits) == 2 and 'auto' not in x_limits:
                        ax.set_xlim(float(x_limits[0]), float(x_limits[1]))
                    
                    if len(y_limits) == 2 and 'auto' not in y_limits:
                        ax.set_ylim(float(y_limits[0]), float(y_limits[1]))
                except ValueError:
                    print("Invalid number format. Please enter two numbers separated by a comma or 'auto'.")
                except Exception as e:
                    print(f"An error occurred: {e}")
                
            def save_png():
                file_path = filedialog.asksaveasfilename(defaultextension='.png', filetypes=[('PNG Image', '*.png')])
                if file_path:
                    fig.savefig(file_path)
                
            def save_csv():
                file_path = filedialog.asksaveasfilename(defaultextension='.csv', filetypes=[('CSV Files', '*.csv')])
                if file_path:
                    x_property = x_var.get()
                    y_property = y_var.get()
                    x_data = [self.shared_data.property_dict[i][x_property] for i in range(len(self.shared_data.property_dict))]
                    y_data = [self.shared_data.property_dict[i][y_property] for i in range(len(self.shared_data.property_dict))]

                    data = np.column_stack((x_data, y_data))
                    np.savetxt(file_path, data, delimiter=',', header=f'{x_property},{y_property}', comments='')

            save_csv_button = tk.Button(grid1_0, text="Save as CSV", command=save_csv)
            save_csv_button.grid(row = 4, column = 1)

            save_png_button = tk.Button(grid1_0, text="Save as PNG", command=save_png)
            save_png_button.grid(row = 4, column = 2)               

            # Curve Fitting Attempt...
            fitting_options = ["Linear", "Polynomial (2)", "Polynomial (3)", "Polynomial (4)"]
            fitting_var = tk.StringVar(self)
            fitting_var.set(fitting_options[0])  # Set the default option

            fitting_label = tk.Label(grid1_1, text="Select Fitting Type:")
            fitting_label.grid(row = 0, column = 0)  

            fitting_dropdown = tk.OptionMenu(grid1_1, fitting_var, *fitting_options)
            fitting_dropdown.grid(row = 1, column = 0)  

            # Function to put our curvefit onto the plot
            def fit_and_plot_data():
                x_property = x_var.get()
                y_property = y_var.get()
                x_data = np.array([self.shared_data.property_dict[i][x_property] for i in range(len(self.shared_data.property_dict))])
                y_data = np.array([self.shared_data.property_dict[i][y_property] for i in range(len(self.shared_data.property_dict))])

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

                elif fitting_var.get() == "Polynomial (3)":
                    fit_params = np.polyfit(x_data, y_data, 3)
                    new_x_data = np.linspace(np.min(x_data) * 0.9, np.max(x_data) * 1.1, 1000)
                    fitted_curve = np.polyval(fit_params, new_x_data)
                    r_squared = 1 - (sum((y_data - np.polyval(fit_params, x_data)) ** 2) / ((len(y_data) - 1) * np.var(y_data, ddof=1)))
                    ax.plot(new_x_data, fitted_curve, label=f"3rd order Polynomial Fit (R^2 = {r_squared:.4f})", color="green")


                elif fitting_var.get() == "Polynomial (4)":
                    fit_params = np.polyfit(x_data, y_data, 4)
                    new_x_data = np.linspace(np.min(x_data) * 0.9, np.max(x_data) * 1.1, 1000)
                    fitted_curve = np.polyval(fit_params, new_x_data)
                    r_squared = 1 - (sum((y_data - np.polyval(fit_params, x_data)) ** 2) / ((len(y_data) - 1) * np.var(y_data, ddof=1)))
                    ax.plot(new_x_data, fitted_curve, label=f"4th order Polynomial Fit (R^2 = {r_squared:.4f})", color="green")

                ax.set_xlabel(x_property)
                ax.set_ylabel(y_property)
                ax.set_title(f'Scatter Plot: {x_property} vs {y_property}')
                ax.legend()
                plt.tight_layout()
                update_axis_limits
                canvas.draw()

            fit_button = tk.Button(grid1_1, text="Fit Data", command=fit_and_plot_data)
            fit_button.grid(row = 2, column = 0) 