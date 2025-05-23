#### Functions for reading SEMACH files and flux calculations ####

#%pip install numpy==1.26.4
#%pip install os
#%pip install pandas
#%pip install scipy
#%pip install tkinter
# load packages ############################################################################################################################
import os
import pandas as pd
import numpy as np
from scipy import stats as st
import tkinter as tk

############################################################################################################################################################
# functions ################################################################################################################################################
############################################################################################################################################################


# read and format SEMACH files ############################################################################################################################
def load_SEMACHv3_file(path):
  meas_data=pd.read_csv(path,\
  sep =";",\
  index_col=False,\
  skiprows=range(1,2),\
  skip_blank_lines=True
  # escape_double = FALSE,\
  # trim_ws = TRUE,\
  #   show_col_types = FALSE))
  )

  meas_data.rename(columns={\
    # original var names
    "GPS (longitude)": "longitude",
    "GPS (latitude)": "latitude",
    "GPS (altitude_m)": "altitude",
    "GPS (satellites)": "gps_satellites",
    "GPS (hdop)":"gps_hdop",
    "GPS (fix)": "gps_signal",
    "Temperature intern (°C)": "temperature_in",
    "Humidity intern (%relH)": "humidity_in",
    "Pressure intern (mbar)": "pressure_in",
    "Temperature extern (°C)": "temperature_ex",
    "Humidity extern (%relH)": "humidity_ex",
    "Pressure extern (mbar)": "pressure_ex",
    "SCD30 CO2 (ppm)": "co2_scd30",
    "SCD30 Temperature (°C)": "temperature_in_scd30",
    "SDC30 Humidity (%relH)": "humidity_in_scd30",
    "Vaisala GMP252 CO2 (ppm)": "co2_gmp252",
    "SMT100 WaterContent (vol %)": "soilmoisture_smt100", 
    "SMT100 Temperature (°C)": "soiltemperature_smt100",     
    "SMT100 Permittivity dielectric coefficient": "soilpermittivity_smt100",
    # new or changed var names
    "Temperature intern (degC)": "temperature_in",
    "Temperature extern (degC)": "temperature_ex",
    "SCD30 Temperature (degC)": "temperature_in_scd30",
    "SMT100 Temperature (degC)": "soiltemperature_smt100", 
    "BatteryLevel (%)":"battery"
    },inplace=True)

  # format time
  meas_data["Timestamp"]=pd.to_datetime(meas_data.Timestamp,format="%d.%m.%Y %H:%M:%S")
  meas_data["meas_time"]=(meas_data.Timestamp-meas_data.Timestamp[0]).dt.total_seconds()


  return(meas_data)



# load multiple files from folder (or single file) and add metadata ############################################################################################################################
# not needed? i think... code for batch loading is integrated in SMEACHv3fluxcalc
def load_SEMACHv3_batch(path:str):
  
  if not path.endswith(".csv"):
    files=os.listdir(path=path)
  else:
    files=path
  
  output=pd.DataFrame()
  #print(files) #debug
  # for future reference / error catching: skip if 0 rows in data file (case when meas. is aborted within first 10sec)
  for file_i in files:
    
    #print(file_i) #debug
    measurement_id=file_i.removesuffix(".csv")
    measurement_data=load_SEMACHv3_file(os.path.join(path,file_i))
    measurement_data["id"]=str(measurement_id)

    output=pd.concat(objs=[output,measurement_data])
  return(output)



# Flux calculation ############################################################################################################################

def SEMACHv3_fluxcalc(path,
                      V_chamber=.01302, #m**3 volume
                      A_chamber=.0434, # m**2 base area (KG250 inner diameter)
                      height_offset=0., #delta, e.g. for base rings
                      cutoff_start=181., # use for eval from (incl.)
                      cutoff_end=3600., # use for eval to (incl.)
                      scd30_scaling=1, #linear calibration for scd30, slope 
                      scd30_offset=1, #linear calibration for scd30, intercept
                      ): 


    # true V_chamber is calculated (delta from hight offset is taken into account)
    V_chamber=V_chamber+A_chamber*height_offset

    output=pd.DataFrame()


    
    if not path.endswith(".csv"):
        # print(2) #debug
        files=os.listdir(path=path)
        
    else:
        # print(2) #debug
        # fixes single file read-in
        files=[os.path.basename(path)]
        path=os.path.dirname(path)

    
    

    for file_i in files:
        print(file_i) #debug
        # load 
        measurement_id=file_i.removesuffix(".csv")
        measurement_data=load_SEMACHv3_file(os.path.join(path,file_i))
        
        print("#### raw ##############################")
        print(measurement_data)
        print(cutoff_start)
        print(cutoff_end)    
        # check total length
        full_duration=measurement_data.meas_time.max()
        
        # check if cutoff end makes sense / if later than last measurement, take all thats there 
        # bugfix local cutoff_time variable
        cutoff_end_i=cutoff_end
        if full_duration<cutoff_end_i:
            cutoff_end_i=full_duration
        print(cutoff_end_i)
        # crop to eval timeframe 
        measurement_data=measurement_data[measurement_data["meas_time"]>=cutoff_start] 
        measurement_data=measurement_data[measurement_data["meas_time"]<=cutoff_end_i]


        # calc available eval duration
        eval_duration=cutoff_end_i-cutoff_start

        
        print(measurement_data) #debug
        



        # avg all support variables (temp, pressure etc.) over evaluation duration , if empty NaN
        output_i=measurement_data.drop(labels=[
                                        "meas_time",
                                        "co2_scd30",
                                        "co2_gmp252"
                                        ],axis=1).describe()[1:2]
        
        
        

        # metadata / fixed params for 
        output_i["id"]=measurement_id
        output_i["full_duration"]=full_duration,
        output_i["eval_duration"]=eval_duration,
        output_i["cutoff_start"]=cutoff_start,
        output_i["cutoff_end"]=cutoff_end_i,
        output_i["height_offset"]=height_offset,
        output_i["scd30_scaling"]=scd30_scaling,
        output_i["scd30_offset"]=scd30_offset
        output_i["V_chamber"]=V_chamber # actual volume not just bare chamber
        output_i["A_chamber"]=A_chamber
        output_i["comment"]=str() #just init, stays empty if no issues
        
        print(output_i)
        
        if not len(output_i): #if no values at all
            output_i["comment"]="No data. " # extra space in case additional errors are appended later
            # fill with nan for appending ... probably a better way of doing this
            #output_i["co2_scd30_corr_tscd"]=float("nan")
            #output_i["co2_scd30_corr"]=float("nan")
            #output_i["co2_gmp252_corr"]=float("nan")
            output_i["f_co2_scd30_corr_tscd"]=float("nan")
            output_i["f_co2_scd30_corr"]=float("nan")
            output_i["f_co2_gmp252"]=float("nan")
            output_i["slope_scd30_tscd"]=float("nan")
            output_i["r2_scd30_tscd"]=float("nan")
            output_i["stderr_scd30_tscd"]=float("nan")
            output_i["slope_scd30"]=float("nan")
            output_i["r2_scd30"]=float("nan")
            output_i["stderr_scd30"]=float("nan")
            output_i["slope_gmp252"]=float("nan")
            output_i["r2_gmp252"]=float("nan")
            output_i["stderr_gmp252"]=float("nan")
        
        else:

            # check data availability after cropping  ### note: something like taht at the beginning to check if any data (omits error when empty file is read)
            if len(measurement_data)<6: # less than ~ 1 min of usable data @ ~ .1 Hz sampling 
                #print(nrow(measurement_data))
            
                output_i["comment"]="Less than 1 min of usable data. No reliable flux calculation possible. "
                # fill with nan for appending
                #measurement_data["co2_scd30_corr_tscd"]=float("nan")
                #measurement_data["co2_scd30_corr"]=float("nan")
                #measurement_data["co2_gmp252_corr"]=float("nan")
                output_i["f_co2_scd30_corr_tscd"]=float("nan")
                output_i["f_co2_scd30_corr"]=float("nan")
                output_i["f_co2_gmp252"]=float("nan")
                output_i["slope_scd30_tscd"]=float("nan")
                output_i["r2_scd30_tscd"]=float("nan")
                output_i["stderr_scd30_tscd"]=float("nan")
                output_i["slope_scd30"]=float("nan")
                output_i["r2_scd30"]=float("nan")
                output_i["stderr_scd30"]=float("nan")
                output_i["slope_gmp252"]=float("nan")
                output_i["r2_gmp252"]=float("nan")
                output_i["stderr_gmp252"]=float("nan")

            else:
              

                # concentration correction and conversion from ppmv to µmol m**-3 
                measurement_data["co2_scd30_corr_tscd"]=(measurement_data.co2_scd30*scd30_scaling+scd30_offset) * \
                measurement_data.pressure_in*100 / ((273.15 + measurement_data.temperature_in_scd30) * 8.314)

                measurement_data["co2_scd30_corr"]=(measurement_data.co2_scd30*scd30_scaling+scd30_offset) * \
                measurement_data.pressure_in*100 / ((273.15 + measurement_data.temperature_in) * 8.314)

                measurement_data["co2_gmp252_corr"]=measurement_data.co2_gmp252 * \
                measurement_data.pressure_in*100 / ((273.15 + measurement_data.temperature_in) * 8.314)

                # slope calculation in µmol sec**-1 m**-3
                #slope_scd30_tscd, intercept_scd30_tscd, r_val_scd30_tscd, p_val_scd30_tscd, std_err_scd30_tscd 
                reg_scd30_tscd = st.linregress(measurement_data.meas_time,measurement_data.co2_scd30_corr_tscd)
                #slope_scd30, intercept_scd30, r_val_scd30, p_val_scd30, std_err_scd30 
                reg_scd30 = st.linregress(measurement_data.meas_time,measurement_data.co2_scd30_corr)
                #slope_gmp252, intercept_gmp252, r_val_gmp252, p_val_gmp252, std_err_gmp252 
                reg_gmp252 = st.linregress(measurement_data.meas_time,measurement_data.co2_gmp252_corr)

                
                output_i["slope_scd30_tscd"] = reg_scd30_tscd.slope
                output_i["r2_scd30_tscd"] = reg_scd30_tscd.rvalue
                output_i["stderr_scd30_tscd"] = reg_scd30_tscd.stderr

                output_i["slope_scd30"] = reg_scd30.slope
                output_i["r2_scd30"] = reg_scd30.rvalue
                output_i["stderr_scd30"] = reg_scd30.stderr

                output_i["slope_gmp252"] = reg_gmp252.slope
                output_i["r2_gmp252"] = reg_gmp252.rvalue
                output_i["stderr_gmp252"] = reg_gmp252.stderr

                if reg_scd30_tscd.rvalue < .9:
                    output_i["comment"]=output_i["comment"]+"Unreliable scd30_tscd slope - r2-flag."
                if reg_scd30_tscd.pvalue > .05:
                    output_i["comment"]=output_i["comment"]+"Unreliable scd30_tscd slope - p_value-flag."

                if reg_scd30.rvalue < .9:
                    output_i["comment"]=output_i["comment"]+"Unreliable scd30 slope - r2-flag."
                if reg_scd30.pvalue > .05:
                    output_i["comment"]=output_i["comment"]+"Unreliable scd30 slope - p_value-flag."

                if reg_gmp252.rvalue < .9:
                    output_i["comment"]=output_i["comment"]+"Unreliable gmp252 slope - r2-flag."
                if reg_gmp252.pvalue > .05:
                    output_i["comment"]=output_i["comment"]+"Unreliable gmp252 slope - p_value-flag."



                # flux calc in µmol sec**-1 m**-2
                output_i["f_co2_scd30_corr_tscd"]=reg_scd30_tscd.slope*V_chamber/A_chamber
                output_i["f_co2_scd30_corr"]=reg_scd30.slope*V_chamber/A_chamber
                output_i["f_co2_gmp252"]=reg_gmp252.slope*V_chamber/A_chamber

        # append to dataset
        output=pd.concat(objs=[output,output_i])
        output.index=output.id
        print(output_i)

           
    return output








##################################################################################################################################################
#### GUI #########################################################################################################################################
# ################################################################################################################################################



import os
import pandas as pd
import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk
class CsvProcessorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("SEMACHv3 app")

        # Parameters init
        self.V = tk.DoubleVar(value=.01302)
        self.A = tk.DoubleVar(value=.0434)
        self.start = tk.IntVar(value=181)
        self.end = tk.IntVar(value=3600)
        self.dh = tk.DoubleVar(value=0.)
        self.scd30_offset = tk.DoubleVar(value=1)
        self.scd30_scaling = tk.DoubleVar(value=1)

        # GUI Components
        self.create_widgets()

    def create_widgets(self):
        # Frame for parameter inputs
        frame = tk.Frame(self.root)
        frame.grid(row=0, column=0, sticky="nsew", pady=10)

        tk.Label(frame, text="V [m3]:").grid(row=0, column=0)
        tk.Entry(frame, textvariable=self.V).grid(row=0, column=1)

        tk.Label(frame, text="A [m2]:").grid(row=1, column=0)
        tk.Entry(frame, textvariable=self.A).grid(row=1, column=1)

        tk.Label(frame, text="Start [sec]:").grid(row=2, column=0)
        tk.Entry(frame, textvariable=self.start).grid(row=2, column=1)

        tk.Label(frame, text="End [sec]:").grid(row=3, column=0)
        tk.Entry(frame, textvariable=self.end).grid(row=3, column=1)

        tk.Label(frame, text="dh [m]:").grid(row=4, column=0)
        tk.Entry(frame, textvariable=self.dh).grid(row=4, column=1)

        tk.Label(frame, text="scd30 slope:").grid(row=5, column=0)
        tk.Entry(frame, textvariable=self.scd30_scaling).grid(row=5, column=1)

        tk.Label(frame, text="scd30 offset:").grid(row=6, column=0)
        tk.Entry(frame, textvariable=self.scd30_offset).grid(row=6, column=1)

        # File selection
        self.filepath = tk.StringVar()
        tk.Entry(frame, textvariable=self.filepath, width=100).grid(row=7, column=0, columnspan=2)
        tk.Button(frame, text="Select Folder", command=self.select_folder).grid(row=8, column=0, columnspan=2)
        tk.Button(frame, text="Select File", command=self.select_file).grid(row=8, column=1, columnspan=2)

        # Process button
        tk.Button(frame, text="Process", command=self.process_files).grid(row=9, column=0, columnspan=2)
        tk.Button(self.root, text="Save as CSV", command=self.save_result).grid(row=9, column=1, columnspan=2)

        # Result display (Treeview)
        self.result_tree = ttk.Treeview(self.root, height=200)
        self.result_tree.grid(row=12, column=0, columnspan=7, rowspan=1, sticky="nsew")

        # Make the window resizable and the layout flexible
        self.root.grid_rowconfigure(0, weight=0)  # Input section does not resize
        self.root.grid_rowconfigure(1, weight=1)  # Table row will take up remaining space

        self.root.grid_columnconfigure(0, weight=1)  # Ensure the first column is resizable
        self.root.grid_columnconfigure(1, weight=1)  # Ensure the second column is resizable

    def select_file(self):
        path = filedialog.askopenfilename(title="Select CSV file", filetypes=[("CSV files", "*.csv")])
        if path:
            self.filepath.set(path)

    def select_folder(self):
        path = filedialog.askdirectory(title="Select Folder containing CSV files")
        if path:
            self.filepath.set(path)

    def process_files(self):
        path = self.filepath.get()
        if not path:
            messagebox.showerror("Error", "Please select a file or folder.")
            return
        try:
            # Calculate the result DataFrame
            result_df = self.calculate_result(path)
            self.display_result(result_df)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def calculate_result(self, path):
        # Placeholder for calculation logic
        outdata = SEMACHv3_fluxcalc(path=path,
                                    V_chamber=self.V.get(),
                                    A_chamber=self.A.get(),
                                    cutoff_start=self.start.get(),
                                    cutoff_end=self.end.get(),
                                    height_offset=self.dh.get(),
                                    scd30_offset=self.scd30_offset.get(),
                                    scd30_scaling=self.scd30_scaling.get())
        return outdata

    def display_result(self, result_df):
        # Clear previous result
        def format_value(x):
            if isinstance(x, (float, int)):
                # add sufficient digits for accuracy and to avoid bindings in future processing
                return f"{x:.6f}"
            elif isinstance(x, pd.Timestamp) or isinstance(x, pd.DatetimeIndex):
                return x.strftime("%d-%m-%Y %H:%M:%S")
            return x

        result_df_formatted = result_df.map(format_value)

        # Clear the existing rows in the treeview
        for item in self.result_tree.get_children():
            self.result_tree.delete(item)

        # Remove previous scrollbars if they exist
        if hasattr(self, 'scrollbar_y'):
            self.scrollbar_y.grid_forget()
        if hasattr(self, 'scrollbar_x'):
            self.scrollbar_x.grid_forget()

        # Set columns based on the result dataframe
        self.result_tree["columns"] = list(result_df_formatted.columns)
        self.result_tree["show"] = "headings"

        # Set headers for each column
        for col in result_df_formatted.columns:
            self.result_tree.heading(col, text=col)  # Sets the column header as the column name

        # Create the vertical scrollbar
        self.scrollbar_y = ttk.Scrollbar(self.result_tree, orient=tk.VERTICAL, command=self.result_tree.yview)
        self.scrollbar_y.grid(row=0, column=2, sticky="ns")
        self.result_tree.configure(yscrollcommand=self.scrollbar_y.set)

        # Create the horizontal scrollbar
        self.scrollbar_x = ttk.Scrollbar(self.result_tree, orient=tk.HORIZONTAL, command=self.result_tree.xview)
        self.scrollbar_x.grid(row=2, column=0, sticky="ew")
        self.result_tree.configure(xscrollcommand=self.scrollbar_x.set)

        # Resize columns dynamically based on content
        for col in result_df_formatted.columns:
            max_len = max(result_df_formatted[col].apply(lambda x: len(str(x))))
            self.result_tree.column(col, width=max_len * 10, anchor="center")

        for index, row in result_df_formatted.iterrows():
            self.result_tree.insert("", "end", values=list(row))

    def save_result(self):
        save_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
        if save_path:
            # Get the current displayed DataFrame
            result_df = self.get_current_result()
            if result_df is not None:
                result_df.to_csv(save_path, index=False)

    def get_current_result(self):
        # Extract data from the Treeview for saving
        if not self.result_tree.get_children():
            return None
        
        columns = self.result_tree["columns"]
        data = []
        for item in self.result_tree.get_children():
            data.append(self.result_tree.item(item)["values"])
        return pd.DataFrame(data, columns=columns)

if __name__ == "__main__":
    root = tk.Tk()
    app = CsvProcessorApp(root)
    root.mainloop()


####################################
''' 
BUGS / ISSUES:
- Horizontal scrolling of prview window does not work
- Short measurements not processed... return NaN / NA entries


FIXED bugs (for future reference)
- flux calculation for GMP252 seems faulty. -> FIXED... did use gmp252 not gmp252_corr for slope calc.
- future warning DataFrame.applymap has been deprecated. Use DataFrame.map instead.  result_df_formatted=result_df.applymap(format_value) ...FIXED
- single file load-in: ERRNO 2: no such file... ...path/folder/filename.csv\\C -> fixed. single str path is read characterwise by for-loop. Also: path set to dirpath for loop  
- when processing: cannot use geometry manager inside which has already slaves managed by pack / pack vs. grid -> fixed using grid only. Looks ugly af now but functional
- preview window does not work at all -> fixed see above
- col names: Ach, Vch values as col names -> fixed... some idiot (me) put the variable as colname, but not with quotation marks
- cutoff time is not reset to desired for multifile read in -> fixed no local variable for cutoff time in loop. changed persistently, affecting subsequent measurements
'''
