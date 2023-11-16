import tkinter as tk
from tkinter import ttk, filedialog
import conformer
import csv
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

class TabThree(ttk.Frame):
    def __init__(self, parent, shared_data,  callback_handler, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.shared_data = shared_data 
        self.callback_handler = callback_handler
        ttk.Frame.__init__(self,parent, *args, **kwargs)
        self.create_tab3_content()

    def create_tab3_content(self):
       
        tab3_label = ttk.Label(self, text="Shape Analysis and Plotting")
        tab3_label.pack()

        self.canvas_histogram = tk.Canvas(self, width=200, height=200)
        self.canvas_histogram.pack(side=tk.LEFT, padx=5, pady=5)
        
        settings_frame = ttk.LabelFrame(self, text="Plot Settings and Options")
        settings_frame.pack(side=tk.RIGHT, padx=5, pady=5)
        
        bin_label = ttk.Label(settings_frame, text="Width of Histogram Bins:")
        bin_label.pack(pady=5)

        self.bin_slider = ttk.Scale(settings_frame, from_=1, to=60, orient="horizontal", command=self.update_bin_label)
        self.bin_slider.set(10)  # Set an initial value
        self.bin_slider.pack()
        
        self.bin_value_label = ttk.Label(settings_frame, text="Bin Width /° : 10")  # Initial value
        self.bin_value_label.pack()
        
        temp_label = ttk.Label(settings_frame, text="Temperature / K:")
        temp_label.pack(pady=5)
        self.temp_entry = ttk.Entry(settings_frame)
        self.temp_entry.insert(0, "298")  # default to 298 K, why not!
        self.temp_entry.pack()

        fit_gaussian_label = ttk.Label(settings_frame, text="Fit Gaussian?") # Gauss fit menu
        fit_gaussian_label.pack(pady=5)
        
        self.fit_gaussian_var = tk.StringVar(settings_frame)
        self.fit_gaussian_var.set("Yes")  # Set an initial value
        
        fit_gaussian_options = ["Yes", "No"]
        fit_gaussian_option = ttk.OptionMenu(settings_frame, self.fit_gaussian_var, "Yes", *fit_gaussian_options)
        fit_gaussian_option.pack()

        save_data_button = ttk.Button(settings_frame, text="Save Data as CSV", command=self.save_data_as_csv)
        save_data_button.pack()
        
        save_data_button = ttk.Button(settings_frame, text="Save Hist. Data as CSV", command=self.save_hist_data_as_csv)
        save_data_button.pack()

        save_plot_button = ttk.Button(settings_frame, text="Save Plot as PNG", command=self.save_plot_as_png)
        save_plot_button.pack()
        
        generate_button = ttk.Button(settings_frame, text="Update Plot", command=self.generate_histograms_and_plots)
        generate_button.pack(pady=10)
        
        self.conformer_angle_stats = ttk.Label(settings_frame, text="")
        self.conformer_angle_stats.pack(pady=10)
            
        self.callback_handler.register_callback("data_updated", self.generate_histograms_and_plots)
            
    def update_bin_label(self, value):
        if int(len(self.shared_data.angle)) >=1: #check its got something in it.
            self.bin_value_label.config(text=f"Bins: {int(float(value))}")  # Convert to int to remove decimals

    def generate_histograms_and_plots(self):
        num_bins = self.bin_slider.get()                        # Get the selected number of bins from the slider
        temperature = float(self.temp_entry.get())              # Get user-input temperature
        
        if self.shared_data.angle == []:
            print('No conformer data present')
            return
            
        if self.shared_data.angle != []:
            self.returned_fig,self.shared_data.delta_energy,self.shared_data.probability,self.hist_x,self.hist_y,xmax,fwhm = conformer.calc_bend_hist(self.shared_data.angle, 
                                                    self.shared_data.energy, 
                                                    fit_gaussian= self.fit_gaussian_var.get() == "Yes", 
                                                    name=self.shared_data.settings_entries[1].get(), 
                                                    Temp=temperature, 
                                                    BinSteps=num_bins)

            for widget in self.canvas_histogram.winfo_children(): 
                widget.destroy() # Destroy the previous canvas widgets to remove the previous plots

            canvas_img = FigureCanvasTkAgg(self.returned_fig, master=self.canvas_histogram)
            canvas_img.get_tk_widget().pack()
            
            info_text = f"Mean Angle: {np.round(np.mean(self.shared_data.angle),2)}°\nMedian Angle: {np.round(np.median(self.shared_data.angle),2)}°\nFit X@Y(max): {xmax}°\nFit FWHM: {fwhm}"
            self.conformer_angle_stats.config(text=info_text)

    def save_data_as_csv(self):

        file_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV Files", "*.csv")])
        if file_path:
            with open(file_path, 'w', newline='') as csvfile:
                csvwriter = csv.writer(csvfile)
                csvwriter.writerow(['Angle / Deg', ' Energy / kj mol', 'dE / kj mol-1', 'Probability'])
                for angle, energy, delta_energy, probability in zip(self.shared_data.angle, self.shared_data.energy, self.shared_data.delta_energy, self.shared_data.probability):
                    csvwriter.writerow([angle, energy, delta_energy, probability])
        print('\nSaved data to ' + str(file_path))
        
    def save_hist_data_as_csv(self):
        file_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV Files", "*.csv")])
        if file_path:
            with open(file_path, 'w', newline='') as csvfile:
                csvwriter = csv.writer(csvfile)
                csvwriter.writerow(['Angle / Deg', 'Probability'])
                for angle, energy in zip(self.hist_x, self.hist_y):
                    csvwriter.writerow([angle, energy])
        print('\nSaved data to ' + str(file_path))
                
    def save_plot_as_png(self):
        file_path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG Image", "*.png")])
        if file_path:
            self.returned_fig.savefig(file_path)
            print('\nSaved figure to ' + str(file_path))
            