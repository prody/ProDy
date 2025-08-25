# -*- coding: utf-8 -*-
__author__ = 'Carlos Ventura'
__email__ = ['carlosventura@stonybrook.edu']

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
import subprocess
import shutil
import glob
import sys
import ast
from pathlib import Path
import importlib

current_dir = Path(__file__).resolve().parent

package_dir = current_dir / "DruGUI-script"

if package_dir.is_dir() and str(package_dir) not in sys.path:
    sys.path.insert(0, str(package_dir)) 

try:
    druggability = importlib.import_module("druggability")
except ImportError as e:
    raise ImportError(
        f"Could not import 'druggability' from {package_dir}. "
        f"({e})"
    )

from prody import LOGGER
from numpy import array, ceil, histogramdd, arange
from prody.proteins.compare import matchAlign
from prody.proteins.pdbfile import parsePDB, writePDB
from prody.measure.measure import calcCenter
from prody.measure.transform import moveAtoms, wrapAtoms
from prody.trajectory.dcdfile import DCDFile, parseDCD
from prody.trajectory.psffile import parsePSF
from prody.trajectory.trajectory import Trajectory

from druggability.grid import OpenDX

__all__ = ['DruGUI']

class DruGUI:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("Druggability GUI Version 2.0")
        self.protein_pdb = tk.StringVar()
        self.protein_psf = tk.StringVar()
        self.probe01_type = tk.StringVar(value = "IPRO")
        self.probe02_type = tk.StringVar(value = "ACTT")
        self.probe03_type = tk.StringVar(value = "BENZ")
        self.probe04_type = tk.StringVar(value = "IMID")
        self.probe05_type = tk.StringVar(value = "IBUT")
        self.probe06_type = tk.StringVar(value = "ACAM")
        self.probe07_type = tk.StringVar(value = "IPAM")
        self.percent_probe01 = tk.IntVar(value=16)
        self.percent_probe02 = tk.IntVar(value=14)
        self.percent_probe03 = tk.IntVar(value=14)
        self.percent_probe04 = tk.IntVar(value=14)
        self.percent_probe05 = tk.IntVar(value=14)
        self.percent_probe06 = tk.IntVar(value=14)
        self.percent_probe07 = tk.IntVar(value=14)
        self.percent_total = tk.IntVar(value=100)
        self.solvent_padding = tk.IntVar(value=15)
        self.boundary_padding = tk.DoubleVar(value=2.4)
        self.neutralize = tk.BooleanVar(value=True)
        self.lipid = tk.BooleanVar(value=False)
        self.write_conf= tk.BooleanVar(value=True)
        self.n_sims = tk.IntVar(value=4)
        self.sim_length = tk.IntVar(value=40)
        self.outdir_location = tk.StringVar()
        self.output_prefix = tk.StringVar()
        self.par_files_list = [
        "probe.prm", "par_all36_cgenff.prm", "par_all27_prot_lipid_na.inp",
        "par_all36_lipid.prm", "par_all36_prot.prm", "par_all36_carb.prm",
        "par_all36_na.prm", "toppar_water_ions.str"
        ]
        self.par_files = tk.StringVar(value='\n'.join(self.par_files_list))
        self.tk_vmd_executable = tk.StringVar(value="/Applications/VMD1.9.4a57-arm64-Rev12.app/Contents/vmd/vmd_MACOSXARM64")
        self.selection = tk.StringVar(value="noh and protein")
        self.grid_spacing = tk.DoubleVar(value=0.5)
        self.contact_distance = tk.DoubleVar(value=4.0)
        self.align = tk.StringVar(value="calpha")
        self.outputdir_location = tk.StringVar()
        self.prefix = tk.StringVar()
        self.temperature = tk.DoubleVar(value=300.)
        self.dia_merge_radius =tk.DoubleVar(value=5.6)
        self.dia_n_frames = tk.IntVar(value=1)
        self.dia_n_probes = tk.IntVar(value=7)
        self.dia_delta_g = tk.IntVar(value=-1.0)
        self.dia_min_n_probes = tk.IntVar(value=6)
        self.dia_low_affinity = tk.IntVar(value=10)
        self.dia_max_charge = tk.IntVar(value=2)
        self.dia_n_solutions = tk.IntVar(value=3)
        self.dia_n_charged = tk.IntVar(value=3)
        self.probe_selection = tk.StringVar(value=[self.probe01_type.get(),self.probe02_type.get(),self.probe03_type.get(),self.probe04_type.get(),self.probe05_type.get(),self.probe06_type.get(),self.probe07_type.get()])
        self.PROBETOPPAR = {
        "PBDA": "probe2.top probe.prm",
        "CGenff": "top_all36_cgenff.rtf par_all36_cgenff.prm"
        }
        self.PROBETYPES = {
        "core": "Core probes",
        "polar": "Polar probes",
        "hydrophobe": "Hydrophobes",
        "negative": "Negatively charged",
        "positive": "Positively charged",
        "ring5": "5-membered rings",
        "ring6": "6-membered rings"
        }
        self.PROBEDATA = {}
        os.chdir('/Users/carlosventura/Desktop/prody_drugui/ProDy/prody/drugui/DruGUI-script')
        self.Druggability_path = os.getcwd()

        # Define the main dropdown menu options and view functions
        self.options = {
            "Prepare System": self.prepare_system_view,
            "Analysis Results": self.analysis_results_view
        }

        self.selected_option = tk.StringVar(value="Prepare System")

        dropdown = ttk.OptionMenu(
            self.root,
            self.selected_option,
            self.selected_option.get(),
            *self.options.keys(),
            command=self.show_view
        )
        dropdown.pack(side='top', fill='x', pady=10)

        self.main_frame = ttk.Frame(self.root)
        self.main_frame.pack(fill='both', expand=True)

        # Show default view
        self.show_view(self.selected_option.get())

        self.root.mainloop()

    def show_view(self, option):
        for widget in self.main_frame.winfo_children():
            widget.destroy()
        self.options[option](self.main_frame)

    def prepare_system_view(self, frame):
        mfa = tk.Frame(frame)
        mfa.grid(row=1, column=0, padx=10, pady=10, sticky='nsew')

        mfaif = tk.LabelFrame(mfa, text="Protein structure and coordinate files:", bd=2)
        mfaif.grid(row=0, column=0, padx=5, pady=5, sticky='ew')

        def show_psf_help():
            messagebox.showinfo(
                "HELP", 
                "Protein PSF file should contain all components of the system "
                "before the solvent and probe molecules are added. This may include "
                "structural/functional ions, cofactors, etc. As a kindly reminder, please also "
                "make sure that the protonation states of histidines, cysteines, or other relevant "
                "residues are set properly and, if any, sulfide bridging cysteines are patched correctly."
            )

        def browse_psf_file():
            tempfile = filedialog.askopenfilename(
                filetypes=[("PSF files", "*.psf"), ("All files", "*.*")]
            )
            if tempfile:
                self.protein_psf.set(tempfile)
                print(f"Selected PSF file: {self.protein_psf.get()}")
            else:
                messagebox.showinfo("No Selection", "No file was selected.")

        psf_help_button = tk.Button(mfaif, text="?", padx=0, pady=0, command=show_psf_help)
        psf_help_button.grid(row=0, column=0, sticky='w')

        psf_label = tk.Label(mfaif, text="PSF:")
        psf_label.grid(row=0, column=1, sticky='w')

        psf_entry = tk.Entry(mfaif, width=60, textvariable=self.protein_psf)
        psf_entry.grid(row=0, column=2, sticky='ew')

        psf_browse_button = tk.Button(mfaif, text="Browse", width=6, pady=1, command=browse_psf_file)
        psf_browse_button.grid(row=0, column=3, sticky='w')

        mfaif.grid_columnconfigure(2, weight=1)

        def show_protein_help():
            messagebox.showinfo(
                "HELP",
                "This file must contain coordinate data corresponding to the "
                "entities described in the PSF file."
                )

        protein_help_button = tk.Button(mfaif, text="?", padx=0, pady=0, command=show_protein_help)
        protein_help_button.grid(row=1, column=0, sticky='w')

        protein_label = tk.Label(mfaif, text="PDB:")
        protein_label.grid(row=1, column=1, sticky='w')

        pdb_entry = tk.Entry(mfaif, width=60, textvariable = self.protein_pdb)
        pdb_entry.grid(row=1, column=2, sticky='w')

        def browse_pdb_file():
            tempfile = filedialog.askopenfilename(
                    filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")]
            )
            if tempfile: 
                self.protein_pdb.set(tempfile)

        pdb_browse_button = tk.Button(mfaif, text="Browse", width=6, pady=1, command=browse_pdb_file)
        pdb_browse_button.grid(row=1, column=3, sticky='w')

        mfahi = tk.LabelFrame(mfa, text="Probe selection:", bd=2, pady=2)
        mfahi.grid(row=2, column=0, padx=5, pady=5, sticky='ew')


        def show_probe_help():
            messagebox.showinfo(
                "HELP", 
                "By default, the core probes will be selected for druggability simulations. "
                "Different probes can be selected. Look at http://prody.csb.pitt.edu/tutorials/drugui_tutorial/cgenff.html for other probe molecules."
            )

        probe_help_button = tk.Button(mfahi, text="?", padx=0, pady=0, command=show_probe_help)
        probe_help_button.grid(row=0, column=0, sticky='w')

        core_label = tk.Label(mfahi, text="Core Probes: ")
        core_label.grid(row=0, column=1, sticky='w')

        # Variable for core probe selection
        core_probe_var = tk.BooleanVar(value=True)
        different_probe_var = tk.BooleanVar(value=False)

        def update_probe_type():
            if core_probe_var.get(): 
                self.probe01_type
                self.probe02_type
                self.probe03_type
                self.probe04_type
                self.probe05_type
                self.probe06_type
                self.probe07_type
                different_probe_var.set(False)
            elif different_probe_var.get():
                self.probe01_type.set("")  # Clear the entry so user can type a different probe
                self.probe02_type.set("")
                self.probe03_type.set("")
                self.probe04_type.set("")
                self.probe05_type.set("")
                self.probe06_type.set("")
                self.probe07_type.set("")

            if core_probe_var.get():  
                self.probe01_type.set("IPRO") 
                self.probe02_type.set("ACTT")
                self.probe03_type.set("BENZ")
                self.probe04_type.set("IMID")
                self.probe05_type.set("IBUT")
                self.probe06_type.set("ACAM")
                self.probe07_type.set("IPAM")
                different_probe_var.set(False)  # Uncheck different probe

        core_probe_check = tk.Checkbutton(mfahi, variable=core_probe_var, command=update_probe_type)
        core_probe_check.grid(row=0, column=2, sticky='w')

        different_label = tk.Label(mfahi, text="Different Probes: ")
        different_label.grid(row=0, column=3, sticky='w')


        different_probe_check = tk.Checkbutton(mfahi, variable=different_probe_var, command=update_probe_type)
        different_probe_check.grid(row=0, column=4, sticky='w')

        mfapo = tk.LabelFrame(mfa, text="Probe composition:", bd=2, pady=2)
        mfapo.grid(row=4, column=0, padx=5, pady=5, sticky='ew')

        def change_percentage(probe_var, value):              
            current = int(probe_var.get() or 0)
            new_value = current + value
            probe_var.set(str(new_value))
            print(f"Current value after change: {probe_var.get()}")

        def update_total():
            """Updates the total percentage when individual probes are changed."""
            total = int(self.percent_probe01.get() or 0) + int(self.percent_probe02.get() or 0) + \
                        int(self.percent_probe03.get() or 0) + int(self.percent_probe04.get() or 0) + \
                        int(self.percent_probe05.get() or 0) + int(self.percent_probe06.get() or 0) + \
                        int(self.percent_probe07.get() or 0)
            self.percent_total.set(total)

        def on_probe_change(probe_var, *args):
            """Callback for when any percent_probe variable changes."""
            # If the value is empty, set it to "0" to avoid errors
            if probe_var.get() == "":
                probe_var.set("0")
            update_total()

        self.percent_probe01.trace("w", lambda *args: on_probe_change(self.percent_probe01, *args))
        self.percent_probe02.trace("w", lambda *args: on_probe_change(self.percent_probe02, *args))
        self.percent_probe03.trace("w", lambda *args: on_probe_change(self.percent_probe03, *args))
        self.percent_probe04.trace("w", lambda *args: on_probe_change(self.percent_probe04, *args))
        self.percent_probe05.trace("w", lambda *args: on_probe_change(self.percent_probe05, *args))
        self.percent_probe06.trace("w", lambda *args: on_probe_change(self.percent_probe06, *args))
        self.percent_probe07.trace("w", lambda *args: on_probe_change(self.percent_probe07, *args))

        #Probe01 

        probe01_label = tk.Label(mfapo, text="% Probe 01:        ")
        probe01_label.grid(row=0, column=0, sticky='w')

        probe01_entry = tk.Entry(mfapo, width=20, textvariable=self.probe01_type)
        probe01_entry.grid(row=0, column=1, sticky='w')

        percent_probe01_entry = tk.Entry(mfapo, width=3, textvariable=self.percent_probe01)
        percent_probe01_entry.grid(row=0, column=2, sticky='w')

        probe01_add10_button = tk.Button(mfapo, text="+10", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe01, 10))
        probe01_add10_button.grid(row=0, column=3, sticky='w')

        probe01_add5_button = tk.Button(mfapo, text="+5", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe01, 5))
        probe01_add5_button.grid(row=0, column=4, sticky='w')

        probe01_add1_button = tk.Button(mfapo, text="+1", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe01, 1))
        probe01_add1_button.grid(row=0, column=5, sticky='w')

        probe01_set0_button = tk.Button(mfapo, text="0", padx=0, pady=0, command=lambda: self.percent_probe01.set("0"))
        probe01_set0_button.grid(row=0, column=6, sticky='w')

        probe01_minus1_button = tk.Button(mfapo, text="-1", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe01, -1))
        probe01_minus1_button.grid(row=0, column=7, sticky='w')

        probe01_minus5_button = tk.Button(mfapo, text="-5", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe01, -5))
        probe01_minus5_button.grid(row=0, column=8, sticky='w')

        probe01_minus10_button = tk.Button(mfapo, text="-10", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe01, -10))
        probe01_minus10_button.grid(row=0, column=9, sticky='w')

        update_total()

        #Probe02

        probe02_label = tk.Label(mfapo, text="% Probe 02:        ")
        probe02_label.grid(row=1, column=0, sticky='w')

        probe02_entry = tk.Entry(mfapo, width=20, textvariable=self.probe02_type)
        probe02_entry.grid(row=1, column=1, sticky='w')

        percent_probe02_entry = tk.Entry(mfapo, width=3, textvariable=self.percent_probe02)
        percent_probe02_entry.grid(row=1, column=2, sticky='w')

        probe02_add10_button = tk.Button(mfapo, text="+10", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe02, 10))
        probe02_add10_button.grid(row=1, column=3, sticky='w')

        probe02_add5_button = tk.Button(mfapo, text="+5", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe02, 5))
        probe02_add5_button.grid(row=1, column=4, sticky='w')

        probe02_add1_button = tk.Button(mfapo, text="+1", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe02, 1))
        probe02_add1_button.grid(row=1, column=5, sticky='w')

        probe02_set0_button = tk.Button(mfapo, text="0", padx=0, pady=0, command=lambda: self.percent_probe02.set("0"))
        probe02_set0_button.grid(row=1, column=6, sticky='w')

        probe02_minus1_button = tk.Button(mfapo, text="-1", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe02, -1))
        probe02_minus1_button.grid(row=1, column=7, sticky='w')

        probe02_minus5_button = tk.Button(mfapo, text="-5", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe02, -5))
        probe02_minus5_button.grid(row=1, column=8, sticky='w')

        probe02_minus10_button = tk.Button(mfapo, text="-10", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe02, -10))
        probe02_minus10_button.grid(row=1, column=9, sticky='w')

        #Probe03

        probe03_label = tk.Label(mfapo, text="% Probe 03:        ")
        probe03_label.grid(row=2, column=0, sticky='w')

        probe03_entry = tk.Entry(mfapo, width=20, textvariable=self.probe03_type)
        probe03_entry.grid(row=2, column=1, sticky='w')

        percent_probe03_entry = tk.Entry(mfapo, width=3, textvariable=self.percent_probe03)
        percent_probe03_entry.grid(row=2, column=2, sticky='w')

        probe03_add10_button = tk.Button(mfapo, text="+10", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe03, 10))
        probe03_add10_button.grid(row=2, column=3, sticky='w')

        probe03_add5_button = tk.Button(mfapo, text="+5", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe03, 5))
        probe03_add5_button.grid(row=2, column=4, sticky='w')

        probe03_add1_button = tk.Button(mfapo, text="+1", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe03, 1))
        probe03_add1_button.grid(row=2, column=5, sticky='w')

        probe03_set0_button = tk.Button(mfapo, text="0", padx=0, pady=0, command=lambda: self.percent_probe03.set("0"))
        probe03_set0_button.grid(row=2, column=6, sticky='w')

        probe03_minus1_button = tk.Button(mfapo, text="-1", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe03, -1))
        probe03_minus1_button.grid(row=2, column=7, sticky='w')

        probe03_minus5_button = tk.Button(mfapo, text="-5", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe03, -5))
        probe03_minus5_button.grid(row=2, column=8, sticky='w')

        probe03_minus10_button = tk.Button(mfapo, text="-10", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe03, -10))
        probe03_minus10_button.grid(row=2, column=9, sticky='w')

        #Probe04

        probe04_label = tk.Label(mfapo, text="% Probe 04:        ")
        probe04_label.grid(row=3, column=0, sticky='w')

        probe04_entry = tk.Entry(mfapo, width=20, textvariable=self.probe04_type)
        probe04_entry.grid(row=3, column=1, sticky='w')

        percent_probe04_entry = tk.Entry(mfapo, width=3, textvariable=self.percent_probe04)
        percent_probe04_entry.grid(row=3, column=2, sticky='w')

        probe04_add10_button = tk.Button(mfapo, text="+10", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe04, 10))
        probe04_add10_button.grid(row=3, column=3, sticky='w')

        probe04_add5_button = tk.Button(mfapo, text="+5", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe04, 5))
        probe04_add5_button.grid(row=3, column=4, sticky='w')

        probe04_add1_button = tk.Button(mfapo, text="+1", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe04, 1))
        probe04_add1_button.grid(row=3, column=5, sticky='w')

        probe04_set0_button = tk.Button(mfapo, text="0", padx=0, pady=0, command=lambda: self.percent_probe04.set("0"))
        probe04_set0_button.grid(row=3, column=6, sticky='w')

        probe04_minus1_button = tk.Button(mfapo, text="-1", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe04, -1))
        probe04_minus1_button.grid(row=3, column=7, sticky='w')

        probe04_minus5_button = tk.Button(mfapo, text="-5", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe04, -5))
        probe04_minus5_button.grid(row=3, column=8, sticky='w')

        probe04_minus10_button = tk.Button(mfapo, text="-10", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe04, -10))
        probe04_minus10_button.grid(row=3, column=9, sticky='w')

        #Probe05

        probe05_label = tk.Label(mfapo, text="% Probe 05:        ")
        probe05_label.grid(row=4, column=0, sticky='w')

        probe05_entry = tk.Entry(mfapo, width=20, textvariable=self.probe05_type)
        probe05_entry.grid(row=4, column=1, sticky='w')

        percent_probe05_entry = tk.Entry(mfapo, width=3, textvariable=self.percent_probe05)
        percent_probe05_entry.grid(row=4, column=2, sticky='w')

        probe05_add10_button = tk.Button(mfapo, text="+10", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe05, 10))
        probe05_add10_button.grid(row=4, column=3, sticky='w')

        probe05_add5_button = tk.Button(mfapo, text="+5", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe05, 5))
        probe05_add5_button.grid(row=4, column=4, sticky='w')

        probe05_add1_button = tk.Button(mfapo, text="+1", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe05, 1))
        probe05_add1_button.grid(row=4, column=5, sticky='w')

        probe05_set0_button = tk.Button(mfapo, text="0", padx=0, pady=0, command=lambda: self.percent_probe05.set("0"))
        probe05_set0_button.grid(row=4, column=6, sticky='w')

        probe05_minus1_button = tk.Button(mfapo, text="-1", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe05, -1))
        probe05_minus1_button.grid(row=4, column=7, sticky='w')

        probe05_minus5_button = tk.Button(mfapo, text="-5", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe05, -5))
        probe05_minus5_button.grid(row=4, column=8, sticky='w')

        probe05_minus10_button = tk.Button(mfapo, text="-10", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe05, -10))
        probe05_minus10_button.grid(row=4, column=9, sticky='w')

        #Probe06

        probe06_label = tk.Label(mfapo, text="% Probe 06:        ")
        probe06_label.grid(row=5, column=0, sticky='w')

        probe06_entry = tk.Entry(mfapo, width=20, textvariable=self.probe06_type)
        probe06_entry.grid(row=5, column=1, sticky='w')

        percent_probe06_entry = tk.Entry(mfapo, width=3, textvariable=self.percent_probe06)
        percent_probe06_entry.grid(row=5, column=2, sticky='w')

        probe06_add10_button = tk.Button(mfapo, text="+10", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe06, 10))
        probe06_add10_button.grid(row=5, column=3, sticky='w')

        probe06_add5_button = tk.Button(mfapo, text="+5", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe06, 5))
        probe06_add5_button.grid(row=5, column=4, sticky='w')

        probe06_add1_button = tk.Button(mfapo, text="+1", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe06, 1))
        probe06_add1_button.grid(row=5, column=5, sticky='w')

        probe06_set0_button = tk.Button(mfapo, text="0", padx=0, pady=0, command=lambda: self.percent_probe06.set("0"))
        probe06_set0_button.grid(row=5, column=6, sticky='w')

        probe06_minus1_button = tk.Button(mfapo, text="-1", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe06, -1))
        probe06_minus1_button.grid(row=5, column=7, sticky='w')

        probe06_minus5_button = tk.Button(mfapo, text="-5", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe06, -5))
        probe06_minus5_button.grid(row=5, column=8, sticky='w')

        probe06_minus10_button = tk.Button(mfapo, text="-10", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe06, -10))
        probe06_minus10_button.grid(row=5, column=9, sticky='w')

        #Probe07

        probe07_label = tk.Label(mfapo, text="% Probe 07:        ")
        probe07_label.grid(row=6, column=0, sticky='w')

        probe07_entry = tk.Entry(mfapo, width=20, textvariable=self.probe07_type)
        probe07_entry.grid(row=6, column=1, sticky='w')

        percent_probe07_entry = tk.Entry(mfapo, width=3, textvariable=self.percent_probe07)
        percent_probe07_entry.grid(row=6, column=2, sticky='w')

        probe07_add10_button = tk.Button(mfapo, text="+10", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe07, 10))
        probe07_add10_button.grid(row=6, column=3, sticky='w')

        probe07_add5_button = tk.Button(mfapo, text="+5", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe07, 5))
        probe07_add5_button.grid(row=6, column=4, sticky='w')

        probe07_add1_button = tk.Button(mfapo, text="+1", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe07, 1))
        probe07_add1_button.grid(row=6, column=5, sticky='w')

        probe07_set0_button = tk.Button(mfapo, text="0", padx=0, pady=0, command=lambda: self.percent_probe07.set("0"))
        probe07_set0_button.grid(row=6, column=6, sticky='w')

        probe07_minus1_button = tk.Button(mfapo, text="-1", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe07, -1))
        probe07_minus1_button.grid(row=6, column=7, sticky='w')

        probe07_minus5_button = tk.Button(mfapo, text="-5", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe07, -5))
        probe07_minus5_button.grid(row=6, column=8, sticky='w')

        probe07_minus10_button = tk.Button(mfapo, text="-10", padx=0, pady=0, command=lambda: change_percentage(self.percent_probe07, -10))
        probe07_minus10_button.grid(row=6, column=9, sticky='w')   

        def total_help():
            messagebox.showinfo(
                "HELP", 
                "Percentages must sum up to 100 or 0. If percentages sum up to 0, "
                "system will only be solvated. If percentages are multiples of 5 or 10, system will "
                "contain lesser number of solvent/probe molecules."
            )

        total_probe_button = tk.Button(mfapo, text="?", padx=0, pady=0, command=total_help)
        total_probe_button.grid(row=7, column=0, sticky='w')

        total_label =tk.Label(mfapo, text="Total of probe percentages: ")
        total_label.grid(row=7, column=1, sticky='w')

        total_percent = tk.Entry(mfapo, width=3, textvariable=self.percent_total)
        total_percent.grid(row=7, column=2, stick='w') 

        mfasi = tk.LabelFrame(mfa, text="Solvation and ionization options:", bd=2, pady=2)
        mfasi.grid(row=5, column=0, padx=5, pady=5, sticky='ew')     

        def padding_help():
            messagebox.showinfo(
                "HELP", 
                "This is the half of the initial distance between the protein "
                "and its imaginary copies under periodic boundary conditions. For systems with "
                "probes, the resulting padding distance will be slightly larger, due to "  
                "constraint of preserving the ratio of 20 water molecules per probe molecule."
            )

        padding_help_button = tk.Button(mfasi, text="?", padx=0, pady=0, command=padding_help)
        padding_help_button.grid(row=0, column=0, sticky='w')

        padding_label = tk.Label(mfasi, text="Simulation box padding (Å): ")
        padding_label.grid(row=0, column=1, sticky="w")

        padding_entry = tk.Entry(mfasi, width=3, textvariable=self.solvent_padding)
        padding_entry.grid(row=0, column=2, sticky='w',)

        def neutralize_help():
            messagebox.showinfo(
                "HELP",
                "By default, counter ions will be added to neutralize a charged "
                "system. A charged system (if the protein is charged) may be obtained by unchecking this option."
            )

        neutralize_help_button = tk.Button(mfasi, text="?", padx=0, pady=0, command=neutralize_help)
        neutralize_help_button.grid(row=0, column=3, sticky='w', padx=8)

        neutralize_label = tk.Label(mfasi, text="Add counter ions: ")
        neutralize_label.grid(row=0, column=4, sticky="w")

        neutralize_check = tk.Checkbutton(mfasi, variable=self.neutralize)
        neutralize_check.grid(row=0, column=5, sticky='w', )

        def boundary_help():
            messagebox.showinfo(
                "HELP",
                "Minimum distance between water/probe and solute."
            )

        boundary_help_button = tk.Button(mfasi, text="?", padx=0, pady=0, command=boundary_help)
        boundary_help_button.grid(row=1, column=0, sticky='w')

        boundary_label = tk.Label(mfasi, text="Boundary (Å): ")
        boundary_label.grid(row=1, column=1, sticky='w')

        boundary_entry = tk.Entry(mfasi, width=3, textvariable=self.boundary_padding)
        boundary_entry.grid(row=1, column=2, sticky='w',)

        def lipid_help():
            messagebox.showinfo(
                "HELP",
                "If the protein is in a membrane, the system is solvated using ",
                "the xy plane of the lipid bilayer."
            )

        lipid_help_button = tk.Button(mfasi, text="?", padx=0, pady=0, command=lipid_help)
        lipid_help_button.grid(row=1, column=3, sticky='w', padx=8)

        lipid_label = tk.Label(mfasi, text="Enclosed in a lipid bilyar: ")
        lipid_label.grid(row=1, column=4, sticky='w')

        lipid_check = tk.Checkbutton(mfasi, variable=self.lipid)
        lipid_check.grid(row=1, column=5,sticky='w')

        mfaoo = tk.LabelFrame(mfa, text="Output options:", bd=2, pady=2)
        mfaoo.grid(row=6, column=0, padx=5, pady=5, sticky='ew') 

        def outdir_help():
            messagebox.showinfo(
                "HELP",
                "Output folder, default is current working directory."
            )

        outdir_help_button = tk.Button(mfaoo, text="?", padx=0, pady=0, command=outdir_help)
        outdir_help_button.grid(row=0, column=0, sticky='w')

        outdir_label = tk.Label(mfaoo, text="Output folder:",)
        outdir_label.grid(row=0, column=1, sticky='w')

        def browse_output_folder():
            folder_path = filedialog.askdirectory(
                title="Select Output Folder"  
            )
            if folder_path:  
                self.outdir_location.set(folder_path)

        outdir_entry = tk.Entry(mfaoo, width=40, textvariable=self.outdir_location)
        outdir_entry.grid(row=0, column=2, sticky='w')

        outdir_browse_button = tk.Button(mfaoo, text="Browse", width=6, pady=1, command=browse_output_folder)
        outdir_browse_button.grid(row=0, column=3, sticky='w')

        def prefix_help():
            messagebox.showinfo(
                "HELP",
                "All output files and folders will start with this prefix. "
                "A unique and descriptive prefix choice may allow running multiple simulations in the same folder."
            )

        prefix_help_button = tk.Button(mfaoo, text="?", padx=0, pady=0, command=prefix_help)
        prefix_help_button.grid(row=1, column=0, sticky='w')

        prefix_label = tk.Label(mfaoo, text="Output prefix:")
        prefix_label.grid(row=1, column=1, sticky='w')

        prefix_path = tk.Entry(mfaoo, width=10, textvariable=self.output_prefix)
        prefix_path.grid(row=1, column=2, sticky='w')

        def write_conf_help():
            messagebox.showinfo(
                "HELP",
                "Minimization, equilibration, and simulation configuration "
                "files, and necessary parameter files will be written. Simulation parameters "
                "cannot be edited by the user within this GUI."
            )

        write_conf_frame = tk.Frame(mfaoo)
        write_conf_frame.grid(row=1, column=3, sticky='w')

        write_conf_help_button = tk.Button(write_conf_frame, text="?", padx=0, pady=0, command=write_conf_help)
        write_conf_help_button.pack(side='left')

        write_conf_label = tk.Label(write_conf_frame, text="Write NAMD input:")
        write_conf_label.pack(side='left')

        write_conf_check = tk.Checkbutton(write_conf_frame, variable=self.write_conf)
        write_conf_check.pack(side='left')

        def nsim_help():
            messagebox.showinfo(
                "HELP",
                "If more than 1 is specified, multiple simulation input files "
                "will be generated. These simulations will differ by their random number "
                "generator seeds. This will result in different trajectories. Multiple simulations "
                "will share output of the same minmization run."
            )

        nsim_help_button = tk.Button(mfaoo, text="?", padx=0, pady=0, command=nsim_help)
        nsim_help_button.grid(row=2, column=0, sticky='w')

        nsim_label = tk.Label(mfaoo, text="Number of sims:")
        nsim_label.grid(row=2, column=1, sticky='w')

        nsim_entry = tk.Entry(mfaoo, width=3, textvariable=self.n_sims)
        nsim_entry.grid(row=2, column=2, sticky='w')

        write_conf_frame1 = tk.Frame(mfaoo)
        write_conf_frame1.grid(row=2, column=3, sticky='w')

        def sim_length_help():
            messagebox.showinfo(
                "HELP",
                "This is the length of the productive run in nanoseconds. "
                "Defaults for minimization (2000 steps) and equilibration (60 ps for water only "
                "systems, 900 ps for probe containing systems) is not included."
            )

        sim_length_help_button = tk.Button(write_conf_frame1, text="?", padx=0, pady=0, command=sim_length_help)
        sim_length_help_button.pack(side='left')

        sim_length_label = tk.Label(write_conf_frame1, text="Sim length (ns):")
        sim_length_label.pack(side='left')

        sim_length_entry = tk.Entry(write_conf_frame1, textvariable=self.sim_length, width = 4)
        sim_length_entry.pack(side='left')

        def par_help():
            messagebox.showinfo(
                "HELP",
                "If a system requires parameters in addition those defined in "
                "par_all36_lipid.prm, additional filenames may be specified here. "
                "Specified files will be copied to parameters folder. If a specified file "
                "does not contain CHARMM format parameters, NAMD runtime error will occur."
            )

        par_help_button = tk.Button(mfaoo, text="?", padx=0, pady=0, command=par_help)
        par_help_button.grid(row=3, column=0, sticky='w')

        par_label = tk.Label(mfaoo, text="Additional\nparameters:")
        par_label.grid(row=3, column=1, rowspan=2, sticky='wn')

        par_frame = tk.Frame(mfaoo)
        par_frame.grid(row=3, column=2, sticky="snew")

        listbox = tk.Listbox(par_frame, 
                            activestyle='dotbox', 
                            listvariable=self.par_files, 
                            selectmode='browse', 
                            width=40, 
                            height=3, 
                            setgrid=True)
        listbox.grid(row=3, column=3, sticky="nsew")

        scroll = tk.Scrollbar(par_frame, command=listbox.yview)
        scroll.grid(row=3, column=4, sticky="ns")

        listbox.config(yscrollcommand=scroll.set)

        self.root.grid_columnconfigure(0, weight=1)
        self.root.grid_rowconfigure(0, weight=1)
        par_frame.grid_columnconfigure(0, weight=1)
        par_frame.grid_rowconfigure(0, weight=1)

        def add_par_files():
            tempfiles = filedialog.askopenfilenames(
                title="Select CHARMM parameter files",
                filetypes=[("CHARMM parameter files", "*.prm *.inp *.str"), ("All files", "*.*")]
            )
            if tempfiles:
                added = False
                for tempfile in tempfiles:
                    if tempfile in self.par_files_list:
                        messagebox.showwarning("WARNING", f"{tempfile} has already been added to the list.")
                    else:
                        self.par_files_list.append(tempfile)
                        added = True
                if added:
                    par_files.set('\n'.join(self.par_files_list))

        test = tk.Frame(mfaoo)
        test.grid(row=3, column=3, sticky='w')

        par_add = tk.Button(test, text="Add", width=6, command=add_par_files, padx=0, pady=0)
        par_add.pack(side='left')

        def remove_par_files():
            selected_indices = listbox.curselection()
            if not selected_indices:
                return
            for i in reversed(selected_indices):
                del self.par_files_list[i]  # Remove from your python list
            par_files.set('\n'.join(self.par_files_list))


        par_delete = tk.Button(test, text="Remove", command=remove_par_files, width=6, padx=0, pady=0)
        par_delete.pack(side='left')

        #The vmd executable will be used to prepare the system for druggability simulations

        def vmd_executable_help():
            messagebox.showinfo(
                "HELP",
                "Specify the path to your computers VMD executable here. "
                "You can use the which function of ProDy. Example: from prody.utilities import which \n which('vmd_MACOSXARM64') and your path location will be printed"
            )

        vmd_executable_help_button = tk.Button(mfaoo, text="?", padx=0, pady=0, command=vmd_executable_help)
        vmd_executable_help_button.grid(row=4, column=0, sticky='w')

        vmd_executable_label = tk.Label(mfaoo, text="VMD Executable:")
        vmd_executable_label.grid(row=4, column=1, rowspan=2, sticky='wn')

        vmd_executabl_entry = tk.Entry(mfaoo, textvariable=self.tk_vmd_executable, width = 40)
        vmd_executabl_entry.grid(row=4, column=2, sticky='w')

        def Prepare_system():
            global par_files
            global percent_probe01 
            global percent_probe02 
            global percent_probe03
            global percent_probe04
            global percent_probe05
            global percent_probe06
            global percent_probe07
            global solvent_padding
            global PROBEDATA
            global PROBETYPES
            global neutralize
            global outdir_location
            neutral = self.neutralize.get()
            global n_sims
            Probes = True
            lipids = self.lipid.get()
            boundary=self.boundary_padding.get()
            padding=self.solvent_padding.get()
            constrain = "heavy"
            vmd = self.tk_vmd_executable.get()
            global write_conf

            drugui_data(self)

            if self.protein_pdb.get() == "" or  self.protein_psf.get() == "":
                messagebox.showerror("ERROR", "Both PSF and PDB files must be specified.")
            else : 
                inputpsf = self.protein_psf.get()
                inputpdb = self.protein_pdb.get()

            if self.percent_total.get() != 100 and self.percent_total.get() != 0:
                messagebox.showerror("ERROR", "Probe percentages must sum up to 100 or 0.\nCurrent total is "+ str(self.percent_total.get()))

            if self.percent_total.get() == 0:
                Probes = False

            os.chdir('/Users/carlosventura/Desktop/prody_drugui/ProDy/prody/drugui/DruGUI-script')

            if package_dir.is_dir():
                os.chdir(package_dir)
            else:
                raise FileNotFoundError(f"DruGUI-script not found at {package_dir}")

            Druggability_path = os.getcwd()

            probepsf = os.path.join(Druggability_path, "probe.psf")
            probepdb = os.path.join(Druggability_path, "probe.pdb")
            probetop = os.path.join(Druggability_path, "probe2.top")
            probeprm = os.path.join(Druggability_path, "probe.prm")
            cgenfftop = os.path.join(Druggability_path, "top_all36_cgenff.rtf")
            cgenffprm = os.path.join(Druggability_path, "par_all36_cgenff.prm")
            probebox = 62.3572
            probekey = "name OH2"

            if not (os.path.exists(probepsf) and os.path.exists(probepdb) and os.path.exists(probetop) and os.path.exists(probeprm) and os.path.exists(cgenfftop) and os.path.exists(cgenffprm)):
                messagebox.showerror(f"ERROR", f"One of the probe PSF, PDB, TOP, or PRM files is not found in {Druggability_path}")

            opts = {}
            parameterfiles =[]
            parm_files = self.par_files.get()

            parameterfiles = ast.literal_eval(parm_files)

            prefix = self.output_prefix.get()

            if prefix == "":
                messagebox.showerror("ERROR", "Prefix must not be left blank")

            if vmd =="":
                messagebox.showerror("ERROR", "The location of your VMD executable is needed to setup druggability simulations.")
                
            drugui_title = self.root.wm_title()
            intermediate = os.path.join(Druggability_path, "intermediate")
            log = open(f"{prefix}.log",'w')
            with open(f"{prefix}.log",'a') as log:
                log.write(f"{drugui_title}\n")
                log.write(f"Input PSF: {self.protein_psf.get()}.\n")
                log.write(f"Input PDB: {self.protein_pdb.get()}.\n")
                log.write(f"Output directory: {self.outdir_location.get()}.\n")
                log.write(f"Output prefix: {prefix}.\n")
                log.write(f"Intermediate = {intermediate}\n")

                if Probes == True:
                    opts = [self.probe01_type.get(), self.probe02_type.get(), self.probe03_type.get(), self.probe04_type.get(), self.probe05_type.get(), self.probe06_type.get(), self.probe07_type.get()]
                    opts_percentage = [self.percent_probe01.get(), self.percent_probe02.get(), self.percent_probe03.get(), self.percent_probe04.get(), self.percent_probe05.get(), self.percent_probe06.get(), self.percent_probe07.get()]
                    log.write(f"Probe compostion: \n")
                    opts_dict = dict(zip(opts, opts_percentage))

                    tcl_opts = "set opts [dict create " + " ".join(f"{key} \"{value}\"" for key, value in opts_dict.items()) + "]"

                    probe_analysis = f"""
                    set PROBEDATA {{{PROBEDATA}}}
                    set PROBETOPPAR {{{self.PROBETOPPAR}}}
                    set PROBETYPES {{{self.PROBETYPES}}}
                    {tcl_opts}
                    set logfile [open probe_analysis.log a]

                    """
                    probe_analysis += """
                    set percent_charge [open "percent_charge.txt" w]

                    dict for {key info} $PROBEDATA {
                            if {[dict exists $opts "$key"]} {
                                dict with info {
                                    set percentage [dict get $opts "$key"]
                                    if {![string is digit $percentage]} {
                                    error "\"$key $percentage\" is not valid, probe percentages must be positive integers."
                                    }
                                    dict unset opts "$key"
                                    if {$source == "CGenFF"} {
                                        set general 1
                                    }           
                                    dict set probePercent $key $percentage
                                    incr probeTotal $percentage
                                    set holder "$key $percentage% ($name; $charge e)"
                                    puts $percent_charge $holder
                                }
                            }
                        }
                    close $percent_charge
                    exit
                    """
                    with open('probe_analysis.tcl', 'w') as tcl_file:
                        tcl_file.write(probe_analysis)
                    subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'probe_analysis.tcl'])
                    with open('percent_charge.txt', 'r') as percent_charge_file:
                        probe_information = percent_charge_file.read().strip()
                    print(probe_information)
                    log.write(f"{probe_information} \n")

                    for probe in opts:
                        if probe not in PROBEDATA:
                            raise ValueError("ERROR", f"Unknown probe '{probe}' not found in PROBEDATA.")

                    for probe_percentage in opts_percentage:
                        if probe_percentage < 0:
                            raise ValueError("ERROR", "Probe percentages must be a positive number")

                    percent_total = sum(opts_percentage)

                    if percent_total != 100:
                        raise ValueError("ERROR", "Probe percentages must sum up to 100.\nCurrent total is "+ str(percent_total))

                else :
                    log.write("Probe molecules are not added to the system.")

                if Probes == True:
                    n = 0
                    probePercent = {}
                    while n<7:
                        probe = opts[n]
                        percentage = opts_percentage[n]
                        probePercent[probe] = []
                        probePercent[probe].append(int(percentage))
                        n+=1
                if neutral == True:
                    log.write("System is neutralized by adding counter ions.\n")
                else :
                    log.write("System is not neutralized by adding counter ions.\n")
                log.write(f"Minimum solvent padding is set to {padding} Å.\n")
                log.write(f"Minimum solvent boundary is to {boundary} Å.\n")
                if lipids == False:
                    log.write(f"System does not contain lipid bilayer.\n")
                else :
                    lipids == True
                    log.write(f"System does contain lipid bilayer.\n")

                log.write(f"NAMD configuration files for {self.n_sims.get()} independent {self.sim_length.get()} ns simuation(s) are written.\n")

                #set padding for the system and probes
                padding_x = padding
                padding_y = padding
                padding_z = padding
                padx = 0
                pady = 0
                padz = 0

                if Probes == True:
                    padx = 5
                    pady = 5
                    padz = 5 
                else :
                    solvent_pad = padding

                if lipids == True:
                    padding_x = 0
                    padding_y = 0
                    padding_z = padding
                    padx = -3
                    pady = -3
                    padz = 5

                solvate_options = f"{inputpsf} {inputpdb} -o {intermediate}"

                if Probes == True:
                    solvate_options += f" -b {boundary} -x {padding_x + padx} -y {padding_y + pady} -z {padding_z + padz}"
                    solvate_options += f" +x {padding_x + padx} +y {padding_y + pady} +z {padding_z + padz}"
                    solvate_options += f' -spdb "{probepdb}" -spsf "{probepsf}" -stop "{probetop}" -ks "{probekey}" -ws {probebox}\n'
                   
                    log.write(f"Command solvate: {solvate_options}")
                   
                    #Solvating the system with VMD
                    tcl_script = f"""
                    package require solvate
                    solvate {solvate_options}
                    set all [atomselect top "protein"]
                    set proteincharge [vecsum [$all get charge]]
                    set pc_output [open "pc_output.txt" w]
                    puts $pc_output $proteincharge
                    close $pc_output
                    exit
                    """
                    with open('solvate_script.tcl', 'w') as tcl_file:
                        tcl_file.write(tcl_script)
                    subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'solvate_script.tcl'])
                    with open('pc_output.txt', 'r') as output_file:
                        proteincharge = output_file.read().strip()
                else :
                    if lipids == True:

                        solvate_options += f" -b {boundary} -x {padding_x + padx} -y {padding_y + pady} -z {padding_z + padz}"
                        solvate_options += f" +x {padding_x + padx} +y {padding_y + pady} +z {padding_z + padz}\n"

                        log.write(f"Command solvate: {solvate_options}")
                        #Solvating the system with VMD
                        tcl_script = f"""
                        package require solvate
                        solvate {solvate_options}
                        set all [atomselect top "protein"]
                        set proteincharge [vecsum [$all get charge]]
                        set pc_output [open "pc_output.txt" w]
                        puts $pc_output $proteincharge
                        close $pc_output
                        exit
                        """
                        with open('solvate_script.tcl', 'w') as tcl_file:
                            tcl_file.write(tcl_script)
                        subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'solvate_script.tcl'])
                        with open('pc_output.txt', 'r') as output_file:
                            proteincharge = float(output_file.read().strip())
                    else:
                        solvate_options += f" -b {boundary} -t {solvent_pad}\n"
                        log.write(f"Command solvate: {solvate_options}")
                        tcl_script = f"""
                        package require solvate
                        solvate {solvate_options}
                        set all [atomselect top "protein"]
                        set proteincharge [vecsum [$all get charge]]
                        set pc_output [open "pc_output.txt" w]
                        puts $pc_output $proteincharge
                        close $pc_output
                        exit
                        """
                        with open('solvate_script.tcl', 'w') as tcl_file:
                            tcl_file.write(tcl_script)
                        subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'solvate_script.tcl'])
                        with open('pc_output.txt', 'r') as output_file:
                            proteincharge = float(output_file.read().strip())

                #Neutraling the system with VMD
                if neutral == True:
                    if Probes == False:
                        log.write(f"Ionization: System has a total charge of {proteincharge} electrons.\n")
                        if proteincharge > 0:
                            nna = 0
                            ncl = round(proteincharge)
                            log.write(f"Ionization: {ncl} chloride ions will be added.\n")
                        else :
                            ncl = 0
                            nna = round(proteincharge)
                            log.write(f"Ionization: {nna} sodium ions will be added.\n")
                        ionization_script = f"""
                        package require autoionize
                        autoionize -psf {intermediate}.psf -pdb {intermediate}.pdb -o {prefix} -from 5 -between 5 -ncl {ncl} -nna {nna} -seg ION
                        exit
                        """
                        with open('ionization_script.tcl', 'w') as tcl_file:
                            tcl_file.write(ionization_script)
                        subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'ionization_script.tcl'])

                if Probes == True:
                #===============================================================================
                # START - ADJUST RELATIVE NUMBER of WATER and PROBE molecules
                # minimum and maximum coordinates for PROTEIN is calculated.
                # if a molecule other than a PROTEIN is of interest, change selection in the next line.
                    protein = 'atomselect top "not water and not segid ION \\"WT.*\\""'
                    minmaxP = f"measure minmax $protein"
                    minP = "lindex $minmaxP 0"
                    maxP = "lindex $minmaxP 1"
                    minPx = "lindex $minP 0"
                    minPy = "lindex $minP 1"
                    minPz = "lindex $minP 2"
                    maxPx = "lindex $maxP 0"
                    maxPy = "lindex $maxP 1"
                    maxPz = "lindex $maxP 2"

                    padx = int(padding_x)
                    pady = int(padding_y)
                    padz = int(padding_z)

                    if lipids == False:
                        selWstring = f"water and name OH2 and x > [expr $minPx-{padx}] and y > [expr $minPy-{pady}] and z > [expr $minPz-{padz}] and x < [expr $maxPx+{padx}] and y < [expr $maxPy+{pady}] and z < [expr $maxPz+{padz}]"
                
                    else :
                        sel = 'set sel [atomselect top "lipid"]'
                        minmaxL = "measure minmax $sel"
                        minLz = "expr [lindex [lindex $minmaxL 0] 2] + 10"
                        maxLz = "expr [lindex [lindex $minmaxL 0] 2] - 10"
                        selWstring = f"water and name OH2 and x > [expr $minPx - {padx}] and y > [expr $minPy - {pady}] and x < [expr $maxPx+{padx}] and y < [expr $maxPy+{pady}] and ((z < [expr $maxPz+{padz}] and z > {maxLz}) or (z < {minLz} and z > [expr $minPz-{padz}]))"
                    # select waters in the box of requested size

                    if lipids == False:
                        selWater = f'atomselect top "{selWstring}"'
                        nWater = "$selWater num"
                        water_script = f"""
                        mol load psf {intermediate}.psf
                        mol load pdb {intermediate}.pdb
                        set protein [{protein}]
                        set minmaxP [{minmaxP}]
                        set minP [{minP}]
                        set maxP [{maxP}]
                        set minPx [{minPx}]
                        set minPy [{minPy}]
                        set minPz [{minPz}]
                        set maxPx [{maxPx}]
                        set maxPy [{maxPy}]
                        set maxPz [{maxPz}]
                        set selWater [{selWater}]
                        set nWater [{nWater}]
                        set indicesWater [$selWater get index]
                        set output_file [open "nwater_output.txt" w]
                        puts $output_file $nWater
                        close $output_file
                        exit
                        """
                    else :
                        selWater = f'atomselect top "{selWstring}"'
                        nWater = "$selWater num"
                        water_script = f"""
                        mol load psf {intermediate}.psf
                        mol load pdb {intermediate}.pdb
                        set protein [{protein}]
                        set minmaxP [{minmaxP}]
                        set minP [{minP}]
                        set maxP [{maxP}]
                        set minPx [{minPx}]
                        set minPy [{minPy}]
                        set minPz [{minPz}]
                        set maxPx [{maxPx}]
                        set maxPy [{maxPy}]
                        set maxPz [{maxPz}]
                        {sel}
                        set minmaxL [{minmaxL}]
                        $sel delete
                        set minLz [{minLz}]
                        set maxLz [{maxLz}]
                        set selWater [{selWater}]
                        set nWater [{nWater}]
                        set indicesWater [$selWater get index]
                        set output_file [open "nwater_output.txt" w]
                        puts $output_file $nWater
                        close $output_file
                        exit
                        """
                    with open('water_script.tcl', 'w') as tcl_file:
                        tcl_file.write(water_script)
                    subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'water_script.tcl'])
                    with open('nwater_output.txt', 'r') as output_file:
                        nWater = int(output_file.read().strip())
                    log.write(f"Number of waters: {nWater}\n")

                    # THE RATIO OF WATER TO PROBE MOLECULE (INTEGERS)
                    # - 20 is an ideal ratio. It has worked fine in the test cases.
                    # - Less than 20 leads to underestimates.
                    if lipids == False:
                        selPstring = f"resname IPRO and name OH2 and x > [expr $minPx-{padx}] and y > [expr $minPy-{pady}] and x < [expr $maxPx+{padx}] and y < [expr $maxPy+{padz}] and z < [expr $maxPz+{padz}] and z > [expr $minPz-{padz}]"
                    else :
                        selPstring = f"resname IPRO and name OH2 and x > [expr $minPx-{padx}] and y > [expr $minPy-{pady}] and x < [expr $maxPx+{padx}] and y < [expr $maxPy+{padz}] and ((z < [expr $maxPz+{padz}] and z > $maxLz) or (z < $minLz and z > [expr $minPz-{padz}]))"
                    selPstring
                    water2probeRatio = 20
    
                    def calc_gcd(p,q):
                    #calculate the greater common divisor
                        while q != 0:
                            p, q = q, p % q
                        return abs(p)
                    gcd = 0
                    for key, value in probePercent.items():
                        for v in value:
                            gcd = calc_gcd(gcd, v)
                    modWater = (water2probeRatio * 100) / gcd
                    log.write(f"Number of waters must be multiples of {modWater}\n")
                    if nWater % modWater == 0:
                        howManyMoreWater = 0
                    else :
                        howManyMoreWater = int(modWater - (nWater % modWater))

                    log.write(f"Change in number of waters: {howManyMoreWater}\n")

                    if howManyMoreWater:
                        pad = 0.1 
                        if lipids == False:
                            addWater = f'[atomselect top "water and name OH2 and exwithin {pad} of index $indicesWater"]'
                        else :
                            addWater = f'[atomselect top "water and name OH2 and exwithin {pad} of index $indicesWater and (z > $maxLz or z < $minLz)"]'

                        if lipids == False:
                            add_water_script = f"""
                            package require psfgen
                            mol load psf {intermediate}.psf
                            mol load pdb {intermediate}.pdb
                            set protein [{protein}]
                            set minmaxP [{minmaxP}]
                            set minP [{minP}]
                            set maxP [{maxP}]
                            set minPx [{minPx}]
                            set minPy [{minPy}]
                            set minPz [{minPz}]
                            set maxPx [{maxPx}]
                            set maxPy [{maxPy}]
                            set maxPz [{maxPz}]
                            set selWater [{selWater}]
                            set nWater {nWater}
                            set indicesWater [$selWater get index]
                            set addWater {addWater}
                            set pad 0.1
                            set padding {padding}
                            set howManyMoreWater {howManyMoreWater}
                            """
                            add_water_script += """
                            while {[$addWater num] < $howManyMoreWater && $pad < [expr {$padding + 5}] } {
                                $addWater delete
                                set pad [expr {$pad + 0.1}]
                                set addWater [atomselect top "water and name OH2 and exwithin $pad of index $indicesWater"]
                            } 
                            set indicesWater "$indicesWater [lrange [$addWater get index] 0 [expr $howManyMoreWater - 1]]"
                            $addWater delete
                            set numWater [llength $indicesWater]
                            set output_file [open "nwater_output.txt" w]
                            puts $output_file $numWater
                            close $output_file
                            set indiceswater_file [open "indiceswater_ouput.txt" w]
                            puts $indiceswater_file $indicesWater
                            close $indiceswater_file
                            exit
                            """
                        else :
                            add_water_script = f"""
                            package require psfgen
                            mol load psf {intermediate}.psf
                            mol load pdb {intermediate}.pdb
                            set protein [{protein}]
                            set minmaxP [{minmaxP}]
                            set minP [{minP}]
                            set maxP [{maxP}]
                            set minPx [{minPx}]
                            set minPy [{minPy}]
                            set minPz [{minPz}]
                            set maxPx [{maxPx}]
                            set maxPy [{maxPy}]
                            set maxPz [{maxPz}]
                            {sel}
                            set minmaxL [{minmaxL}]
                            $sel delete
                            set minLz [{minLz}]
                            set maxLz [{maxLz}]
                            set selWater [{selWater}]
                            set nWater {nWater}
                            set indicesWater [$selWater get index]
                            set addWater {addWater}
                            set pad 0.1
                            set padding {padding}
                            set howManyMoreWater {howManyMoreWater}
                            """
                            add_water_script += """
                            while {[$addWater num] < $howManyMoreWater && $pad < [expr {$padding + 5}] } {
                                $addWater delete
                                set pad [expr {$pad + 0.1}]
                                set addWater [atomselect top "water and name OH2 and exwithin $pad of index $indicesWater and (z > $maxLz or z < $minLz)"]
                            } 
                            set indicesWater "$indicesWater [lrange [$addWater get index] 0 [expr $howManyMoreWater - 1]]"
                            $addWater delete
                            set numWater [llength $indicesWater]
                            set output_file [open "nwater_output.txt" w]
                            puts $output_file $numWater
                            close $output_file
                            set indiceswater_file [open "indiceswater_ouput.txt" w]
                            puts $indiceswater_file $indicesWater
                            close $indiceswater_file
                            exit
                            """
                        with open('add_water_script.tcl', 'w') as tcl_file:
                            tcl_file.write(add_water_script)
                        subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'add_water_script.tcl'])
                        with open('nwater_output.txt', 'r') as output_file:
                            numWater = int(output_file.read().strip())
                        with open('indiceswater_ouput.txt', 'r') as waterindice_file:
                            indicesWater = waterindice_file.read().strip()
                        log.write(f"Final number of waters: {numWater}\n")

                    #Select Probes
                    if lipids == False:
                        selProbe = f'[atomselect top "{selPstring}"]'
                        nProbe = "$selProbe num"

                        probe_script = f"""
                        mol load psf {intermediate}.psf
                        mol load pdb {intermediate}.pdb
                        set protein [{protein}]
                        set minmaxP [{minmaxP}]
                        set minP [{minP}]
                        set maxP [{maxP}]
                        set minPx [{minPx}]
                        set minPy [{minPy}]
                        set minPz [{minPz}]
                        set maxPx [{maxPx}]
                        set maxPy [{maxPy}]
                        set maxPz [{maxPz}]
                        set selProbe {selProbe}
                        set nProbe [{nProbe}]
                        set indicesProbe [$selProbe get index]
                        set output_file [open "nprobe_output.txt" w]
                        puts $output_file $nProbe
                        close $output_file
                        set indices_file [open "indices_output.txt" w]
                        puts $indices_file $indicesProbe
                        close $indices_file
                        set Probelength [llength $indicesProbe]
                        exit
                        """
                    else :
                        selProbe = f'[atomselect top "{selPstring}"]'
                        nProbe = "$selProbe num"

                        probe_script = f"""
                        mol load psf {intermediate}.psf
                        mol load pdb {intermediate}.pdb
                        set protein [{protein}]
                        set minmaxP [{minmaxP}]
                        set minP [{minP}]
                        set maxP [{maxP}]
                        set minPx [{minPx}]
                        set minPy [{minPy}]
                        set minPz [{minPz}]
                        set maxPx [{maxPx}]
                        set maxPy [{maxPy}]
                        set maxPz [{maxPz}]
                        {sel}
                        set minmaxL [{minmaxL}]
                        $sel delete
                        set minLz [{minLz}]
                        set maxLz [{maxLz}]
                        set selProbe {selProbe}
                        set nProbe [{nProbe}]
                        set indicesProbe [$selProbe get index]
                        set output_file [open "nprobe_output.txt" w]
                        puts $output_file $nProbe
                        close $output_file
                        set indices_file [open "indices_output.txt" w]
                        puts $indices_file $indicesProbe
                        close $indices_file
                        set Probelength [llength $indicesProbe]
                        exit
                        """

                    with open('probe_sel_script.tcl', 'w') as tcl_file:
                        tcl_file.write(probe_script)
                    subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'probe_sel_script.tcl'])
                    with open('indices_output.txt', 'r') as output_file1:
                        indicesProbe = output_file1.read().strip()
                    with open('nprobe_output.txt', 'r') as output_file:
                        nProbe = int(output_file.read().strip())
                    log.write(f"Number of probes: {nProbe}\n")

                    #Selection of final probe number
                    howManyMoreProbe = int((numWater / water2probeRatio) - nProbe)
                    log.write(f"Change in number of probes: {howManyMoreProbe}\n")

                    if howManyMoreProbe > 0:
                        pad = 5.0
                        if lipids == False:
                            addProbe = f'[atomselect top "resname IPRO and name OH2 and exwithin {pad} of index {nProbe}"]'
                        else :
                            addProbe = f'[atomselect top "resname IPRO and name OH2 and exwithin {pad} of index {nProbe} and (z > $maxLz or z < $minLz)"]'

                        if lipids == False:
                            add_probe_script = f"""
                            mol load psf {intermediate}.psf
                            mol load pdb {intermediate}.pdb
                            set protein [{protein}]
                            set minmaxP [{minmaxP}]
                            set minP [{minP}]
                            set maxP [{maxP}]
                            set minPx [{minPx}]
                            set minPy [{minPy}]
                            set minPz [{minPz}]
                            set maxPx [{maxPx}]
                            set maxPy [{maxPy}]
                            set maxPz [{maxPz}]
                            set selProbe {selProbe}
                            set nProbe [$selProbe num]
                            set indicesProbe [$selProbe get index]
                            set padding {padding}
                            set howManyMoreProbe {howManyMoreProbe}
                            set addProbe {addProbe}
                            """
                            add_probe_script += """
                            if {$howManyMoreProbe > 0} {
                                set pad 5.0
                                while {[$addProbe num] < $howManyMoreProbe && $pad < [expr {$padding + 5}] } {
                                    $addProbe delete
                                    set pad [expr {$pad + 0.25}]
                                    set addProbe [atomselect top "resname IPRO and name OH2 and exwithin $pad of index $indicesProbe"]
                                }
                                set indicesProbe "$indicesProbe [lrange [$addProbe get index] 0 [expr $howManyMoreProbe - 1]]"
                                $addProbe delete
                            } elseif {$howManyMoreProbe < 0} {
                                set indicesProbe [lrange $indicesProbe 0 end+$howManyMoreProbe]
                            }
                            set numProbe [llength $indicesProbe]
                            set indicesProbe [lsort $indicesProbe]
                            set output_file [open "nprobe_output.txt" w]
                            puts $output_file $numProbe
                            close $output_file
                            set indices_file [open "indices_files.txt" w]
                            puts $indices_file $indicesProbe
                            close $indices_file
                            exit
                            """
                        else :
                            add_probe_script = f"""
                            mol load psf {intermediate}.psf
                            mol load pdb {intermediate}.pdb
                            set protein [{protein}]
                            set minmaxP [{minmaxP}]
                            set minP [{minP}]
                            set maxP [{maxP}]
                            set minPx [{minPx}]
                            set minPy [{minPy}]
                            set minPz [{minPz}]
                            set maxPx [{maxPx}]
                            set maxPy [{maxPy}]
                            set maxPz [{maxPz}]
                            {sel}
                            set minmaxL [{minmaxL}]
                            $sel delete
                            set minLz [{minLz}]
                            set maxLz [{maxLz}]
                            set selProbe {selProbe}
                            set nProbe [$selProbe num]
                            set indicesProbe [$selProbe get index]
                            set padding {padding}
                            set howManyMoreProbe {howManyMoreProbe}
                            set addProbe {addProbe}
                            """
                            add_probe_script += """
                            if {$howManyMoreProbe > 0} {
                                set pad 5.0
                                while {[$addProbe num] < $howManyMoreProbe && $pad < [expr {$padding + 5}] } {
                                    $addProbe delete
                                    set pad [expr {$pad + 0.25}]
                                    set addProbe [atomselect top "resname IPRO and name OH2 and exwithin $pad of index $indicesProbe"]
                                }
                                set indicesProbe "$indicesProbe [lrange [$addProbe get index] 0 [expr $howManyMoreProbe - 1]]"
                                $addProbe delete
                            } elseif {$howManyMoreProbe < 0} {
                                set indicesProbe [lrange $indicesProbe 0 end+$howManyMoreProbe]
                            }
                            set numProbe [llength $indicesProbe]
                            set indicesProbe [lsort $indicesProbe]
                            set output_file [open "nprobe_output.txt" w]
                            puts $output_file $numProbe
                            close $output_file
                            set indices_file [open "indices_files.txt" w]
                            puts $indices_file $indicesProbe
                            close $indices_file
                            exit
                            """
                        with open('add_probe_script.tcl', 'w') as tcl_file:
                            tcl_file.write(add_probe_script)
                        subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'add_probe_script.tcl'])
                        with open('nprobe_output.txt', 'r') as output_file:
                            numProbe = int(output_file.read().strip())
                        with open('indices_files.txt', 'r') as output_file1:
                            indicesProbe = output_file1.read().strip()
                        log.write(f"Final number of probes: {numProbe}\n")
                        log.write(f"System contains {numWater} water and {numProbe} probe molecules\n")

                    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                    # END - ADJUST RELATIVE NUMBER of WATER and PROBE molecules

                    # WRITE PDB files for SOLVENT and IONS
                    # PSFGEN


                    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
                    # START - PROBE RENUMBER and MUTATE
                    # Renumber probe molecule resids starting from 1
                    # This is useful when mutating probe molecules
                    log.write(f"System contains the following probes: \n")
                    probecharge = 0
                    tcl_probepercent = "set probePercent [dict create " + " ".join(f"{key} \"{value}\"" for key, value in opts_dict.items()) + "]"

                    probe_charge = f"""
                    set probeCount [dict create]
                    set howmanylist [list]
                    set probeidlist [list]
                    set aliaslist [list]
                    set nProbe {numProbe}
                    set PROBEDATA {{{PROBEDATA}}}
                    set probecharge {probecharge}
                    {tcl_probepercent}
                    """
                    probe_charge +="""
                    set probe_numcharge [open "probe_numcharge.txt" w]
                    dict for {key value} $probePercent {
                            set howmanyPROB [::tcl::mathfunc::int [expr $nProbe * $value / 100.0]]
                            dict set probeCount $key $howmanyPROB
                            set holder "$howmanyPROB [dict get $PROBEDATA $key name] molecules ($key)"
                            puts $probe_numcharge $holder
                            lappend probeidlist $key
                            lappend howmanylist $howmanyPROB
                            lappend aliaslist [dict get $PROBEDATA $key alias]
                            set charge [::tcl::mathfunc::int [dict get $PROBEDATA $key charge]]
                            if {$charge} {
                                incr probecharge [expr $howmanyPROB * $charge]
                            }
                        }
                    close $probe_numcharge
                    set output [open "probe_charge.txt" w]
                    puts $output $probecharge
                    close $output

                    set howmanyoutput [open "howmanylist.txt" w]
                    puts $howmanyoutput $howmanylist
                    close $howmanyoutput

                    set probeidoutput [open "probeidlist.txt" w]
                    puts $probeidoutput $probeidlist
                    close $probeidoutput

                    set aliasoutput [open "aliaslist.txt" w]
                    puts $aliasoutput $aliaslist
                    close $aliasoutput

                    exit
                    """
                    with open('probe_charge.tcl', 'w') as tcl_file:
                        tcl_file.write(probe_charge)
                    subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'probe_charge.tcl'])
                    with open('probe_numcharge.txt', 'r') as probe_numcharge_file:
                        probe_amount = probe_numcharge_file.read().strip()
                    log.write(f"{probe_amount} \n")
                    with open('probe_charge.txt', 'r') as probe_charge_file:
                        probecharge = float(probe_charge_file.read().strip())
                    with open('howmanylist.txt', 'r') as howmanylist_file:
                        howmanylist = howmanylist_file.read().strip()
                    with open('probeidlist.txt', 'r') as probeidlist_file:
                        probeidlist = probeidlist_file.read().strip()
                    with open('aliaslist.txt', 'r') as aliaslist_file:
                        aliaslist = aliaslist_file.read().strip()

                    totalcharge = (probecharge + float(proteincharge))
                    log.write(f"The system has a total charge of {totalcharge} electrons.\n")

                    n_ions = 0

                    if totalcharge > 0:
                        n_ions = round(totalcharge)
                        log.write(f"Ionization: {n_ions} chloride ions will be added.\n")
                    else :
                        n_ions = (round(totalcharge) * -1)
                        log.write(f"Ionization: {n_ions} sodium ions will be added.\n")

                    probe_script = f"""
                    package require psfgen
                    set intermediate {intermediate}
                    set prefix {prefix}
                    mol load psf intermediate.psf
                    mol load pdb intermediate.pdb
                    psfcontext [psfcontext create]
                    topology {cgenfftop}
                    topology {probetop}
                    readpsf {inputpsf}
                    set sel [atomselect top "not segid ION \\"WT.*\\""]
                    $sel writepdb $intermediate.pdb 
                    $sel delete
                    coordpdb $intermediate.pdb
                    set residProbe 1
                    set indicesProbe {{{indicesProbe}}}
                    set PROBEDATA [dict create {PROBEDATA}]
                    {tcl_probepercent}
                    set nProbe {numProbe}
                    set probeidlist {{{probeidlist}}}
                    set howmanylist {{{howmanylist}}}
                    set aliaslist {{{aliaslist}}}
                    set totalcharge {totalcharge}
                    set n_ions {n_ions}
                    set indicesWater {{{indicesWater}}}
                    set neutral {{{neutral}}}
                    """
                    probe_script += """
                    foreach indexProbe $indicesProbe {
                        set sel [atomselect top "same residue as index $indexProbe"]
                        $sel set resid $residProbe
                        $sel set chain " X"
                        incr residProbe
                        $sel delete
                    }
                    if {![dict exist $probePercent 'IPRO'] || [dict get $probePercent 'IPRO'] < 100} {
                        set residProbe 1
                        while {$residProbe <= $nProbe} {
                            set whichProbe [::tcl::mathfunc::int [expr rand() * [llength $probeidlist]]]
                            if {[lindex $howmanylist $whichProbe] > 0} {
                                set sel [atomselect top "chain  X and resid $residProbe"]
                                $sel set resname [lindex $probeidlist $whichProbe]
                                $sel delete
                                foreach {old new} [lindex $aliaslist $whichProbe] {
                                set sel [atomselect top "chain  X and resid $residProbe and name $old"]
                                $sel set name $new
                                $sel delete
                                }
                                incr residProbe
                                lset howmanylist $whichProbe [expr [lindex $howmanylist $whichProbe] -1]
                            }
                        }
                        set selstr [list]
                        dict for {key value} $probePercent {
                            set info [dict get $PROBEDATA $key]
                            dict with info {
                            lappend selstr "(resname $key and name $atomnames)"
                            }
                        }
                        set selstr [join $selstr " or "]
                        set sel [atomselect top "(same residue as index $indicesProbe) and ($selstr)"]
                        } else {
                            set sel [atomselect top "same residue as index $indicesProbe"]
                        }
                    $sel writepdb $intermediate.pdb
                    $sel delete
                    set residProbe 1
                    segment PROB { pdb $intermediate.pdb }
                    coordpdb $intermediate.pdb PROB    

                    if {$totalcharge > 0} {
                        set n_ions [expr round($totalcharge)]
                        set ion_name "CLA"
                        set ion_resname "CLA"
                        set ncl $n_ions
                        set nna 0
                    } else  {
                        set n_ions [expr -1 * round($totalcharge)]
                        set ion_name "SOD"
                        set ion_resname "SOD"
                        set ncl 0
                        set nna $n_ions
                    }

                    set sel [atomselect top "segid \\"WT.*\\""]
                    set segidWTs [lsort -unique [$sel get segid]]
                    $sel delete
                    foreach segidWT $segidWTs {
                        set sel [atomselect top "segid $segidWT and index $indicesWater"]
                        # While at it, renumber water molecule resids starting from 1 for each segment
                        set residWater 1
                        foreach indexWater [$sel get index] {
                            set sel [atomselect top "same residue as index $indexWater"]
                            $sel set resid $residWater
                            incr residWater
                            $sel delete
                        }
                        set sel [atomselect top "segid $segidWT and (same residue as index $indicesWater)"]
                        $sel writepdb $intermediate.pdb
                        segment $segidWT {pdb $intermediate.pdb}
                        coordpdb $intermediate.pdb $segidWT
                        $sel delete
                    }

                    guesscoord

                    writepsf $intermediate.psf
                    writepdb $intermediate.pdb

                    if {$neutral}{
                        package require autoionize
                        autoionize -psf $intermediate.psf -pdb $intermediate.pdb -o $prefix -from 5 -between 5 -ncl $ncl -nna $nna -seg ION
                    }

                    exit
                    """
                    with open('probe_script.tcl', 'w') as tcl_file:
                        tcl_file.write(probe_script)
                    subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'probe_script.tcl'])

                log.write(f"Output: Structural data is written into {prefix}.psf file. \n")
                log.write(f"Output: Structural data is written into {prefix}.pdb file. \n")

                if constrain:
                    constrain_script = f"""
                    mol load psf {prefix}.psf
                    mol load pdb {prefix}.pdb
                    set prefix {prefix}

                    set all [atomselect top "all"]
                    $all set beta 0
                    $all set occupancy 0
                    # protein heavy atoms BETA 1
                    $all delete
                    set constrain {constrain}
                    """
                    constrain_script += """
                    if {$constrain == "heavy"} {
                    set protein [atomselect top "noh and not water and not segid PROB ION \\"WT.*\\""]
                    set protein_num [$protein num]
                    }
                    set output_file [open "protein_num_constrain.txt" w]
                    puts $output_file $protein_num
                    close $output_file
                    $protein set beta 1
                    $protein delete
                    # alpha carbons OCCUPANCY 1
                    set protein [atomselect top "protein and name CA and not segid PROB ION \\"WT.*\\""]
                    $protein set occupancy 1
                    set geomcent [measure center $protein]
                    $protein delete
                    set all [atomselect top "all"]
                    $all writepdb $prefix.pdb
                    $all delete
                    set geomcent_output [open "geomcent_output.txt" w]
                    puts $geomcent_output $geomcent
                    close $geomcent_output
                    exit
                    """
                    with open('constrain_script.tcl', 'w') as tcl_file:
                        tcl_file.write(constrain_script)
                    subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'constrain_script.tcl'])
                    with open('protein_num_constrain.txt', 'r') as output_file:
                        protein_num = int(output_file.read().strip())
                    with open('geomcent_output.txt', 'r') as geo_output_file:
                        geomcent = geo_output_file.read().strip()
                    log.write(f"Constraints: {protein_num} heavy atoms are constrained in equilibration\n")

                #XSC File generation
                xsc_script = f"""
                mol load psf {prefix}.psf
                mol load pdb {prefix}.pdb

                set selWater [atomselect top "noh water"]
                set minmaxW [measure minmax $selWater]
                set minW [lindex $minmaxW 0]
                set maxW [lindex $minmaxW 1]
                set minWx [lindex $minW 0]
                set minWy [lindex $minW 1]
                set minWz [lindex $minW 2]
                set maxWx [lindex $maxW 0]
                set maxWy [lindex $maxW 1]
                set maxWz [lindex $maxW 2]
                set geomcent {{{geomcent}}}

                set desired_density 0.62

                set total_mass [vecsum [[atomselect top "all"] get mass]]
                set dimScale [::tcl::mathfunc::pow [expr $total_mass / $desired_density / ($maxWx - $minWx) / ($maxWy - $minWy) / ($maxWz - $minWz)] [expr 1.0 / 3.0]]
                set xLength [expr ($maxWx - $minWx)*$dimScale]
                set yLength [expr ($maxWy - $minWy)*$dimScale]
                set zLength [expr ($maxWz - $minWz)*$dimScale]

                set xsc_file [open "{prefix}.xsc" w]
                puts $xsc_file "# NAMD extended system configuration output file"
                puts $xsc_file "#\$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z"
                puts $xsc_file "0 $xLength  0  0 0 $yLength  0 0  0 $zLength [join $geomcent]" 
                close $xsc_file

                set atomtotal [[atomselect top all] num]
                set num_atom [open "numatom_ouput.txt" w]
                puts $num_atom $atomtotal
                close $num_atom 

                set mass_total [open "masstotal_output.txt" w]
                puts $mass_total $total_mass
                close $mass_total

                set density [expr $total_mass / $xLength / $yLength / $zLength]
                set dens [open "density_output.txt" w]
                puts $dens $density
                close $dens
                exit
                """

                log.write(f"Output: Extended system coordinates are written into {prefix}.xsc file. \n")

                with open('xsc_script.tcl', 'w') as tcl_file:
                    tcl_file.write(xsc_script)
                subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'xsc_script.tcl'])
                with open('numatom_ouput.txt', 'r') as output_file:
                    numAtom = int(output_file.read().strip())
                log.write(f"Statistics: System contains {numAtom} atoms\n")
                with open('masstotal_output.txt', 'r') as output_file:
                    massTotal = float(output_file.read().strip())
                log.write(f"Statistics: Mass of the system is {massTotal} amu\n")
                with open('density_output.txt', 'r') as output_file:
                    density = float(output_file.read().strip())
                log.write(f"Statistics: Density of the system is {density} amu/A^3\n")

                parfolder = self.outdir_location.get() + f'/{prefix}/parameters'
                par_lines = []
                par_lines.append("paraTypeCharmm  on")

                par_files = [] 

                for par_file in parameterfiles:
                    par_files = os.path.basename(par_file)
                    par_files_location = os.path.join(Druggability_path, par_files)
                    par_files_path = os.path.join(parfolder,par_files)

                    if not os.path.exists(parfolder):
                        os.makedirs(parfolder)

                    if not os.path.exists(par_files_path):
                        shutil.copy(par_files_location, parfolder)
                    par_lines.append(f"parameters      ../parameters/{par_files}")

                par_lines = "\n".join(par_lines)
                log.write("Simulation: Parameter files are copied into /parameter folder.\n")
                log.close()

                if self.write_conf:
                    nonbonded = []
                    nonbonded.append("cutoff          10.0")
                    nonbonded.append("switching       on")
                    nonbonded.append("switchdist      8.0")  
                    nonbonded.append("pairlistdist    12.0")
                    nonbonded.append("margin          1.0")
                    nonbonded.append("exclude         scaled1-4")
                    nonbonded = "\n".join(nonbonded)

                    minfix = "min"
                    min = open(f"min.conf",'w')
                    min.write(f"coordinates     ./{prefix}.pdb\n")
                    min.write(f"structure       ./{prefix}.psf\n")
                    min.write(f"{par_lines}\n")
                    min.write(f"outputname      {minfix}\n")
                    min.write(f"binaryoutput    no\n")
                    min.write(f"restartname     {minfix}\n")
                    min.write(f"restartfreq     10000\n")
                    min.write(f"timestep        1.0\n")
                    min.write(f"{nonbonded}\n")
                    min.write(f"temperature     0\n")
                    min.write(f"constraints     on\n")
                    min.write(f"consref         ./{prefix}.pdb\n")
                    min.write(f"conskfile       ./{prefix}.pdb\n")
                    min.write(f"conskcol        B\n")
                    min.write(f"constraintScaling  1.0\n")
                    min.write(f"PME             yes\n")
                    min.write(f"PMEGridSpacing  1.0\n")
                    min.write(f"extendedSystem  ./{prefix}.xsc\n")
                    min.write(f"wrapWater       on\n")
                    min.write(f"wrapAll         on\n")
                    if lipids == True:
                        min.write(f"minimize        20000\n")
                    else :
                        min.write(f"minimize        10000\n")
                    min.close()

                    eq1 = open(f"eq1.conf",'w')
                    eq1.write(f"coordinates     ./min.coor\n")
                    eq1.write(f"structure       ./{prefix}.psf\n")
                    eq1.write(f"{par_lines}\n")
                    eq1.write(f"outputname      eq1\n")
                    eq1.write(f"binaryoutput    no\n")
                    eq1.write(f"restartname     eq1\n")
                    eq1.write(f"restartfreq     2000\n")
                    eq1.write(f"binaryrestart   no\n")
                    eq1.write(f"DCDfreq         2000\n")
                    eq1.write(f"DCDfile         eq1.dcd\n")
                    eq1.write(f"outputEnergies  2000\n")
                    if lipids == True:
                        eq1.write(f"timestep        1.0\n")
                    else :
                        eq1.write(f"timestep        2.0\n")
                    eq1.write(f"fullElectFrequency 2\n")
                    eq1.write(f"nonbondedFreq   1\n")
                    eq1.write(f"rigidBonds      all\n")
                    eq1.write(f"{nonbonded}\n")
                    eq1.write(f"temperature     0\n")
                    eq1.write(f"constraints     on\n")
                    eq1.write(f"consref         ./min.coor\n")
                    eq1.write(f"conskfile       ./{prefix}.pdb\n")
                    eq1.write(f"conskcol        B\n")
                    eq1.write(f"constraintScaling  0.5\n")
                    eq1.write(f"PME             yes\n")
                    eq1.write(f"PMEGridSpacing  1.0\n")
                    eq1.write(f"langevin        on\n")
                    eq1.write(f"langevinTemp    100\n")
                    eq1.write(f"langevinDamping 5\n")
                    eq1.write(f"langevinHydrogen off\n")
                    eq1.write(f"useGroupPressure yes\n")
                    eq1.write(f"useFlexibleCell  no\n")
                    eq1.write(f"useConstantArea  no\n")
                    eq1.write(f"useConstantRatio no\n")
                    eq1.write(f"langevinPiston   on\n")
                    eq1.write(f"langevinPistonTarget  1.01325\n")
                    eq1.write(f"langevinPistonPeriod  100.0\n")
                    eq1.write(f"langevinPistonDecay   50.0\n")
                    eq1.write(f"langevinPistonTemp    100\n")
                    eq1.write(f"extendedSystem    min.xsc\n")
                    eq1.write(f"wrapWater       on\n")
                    eq1.write(f"wrapAll         on\n")
                    eq1.write(f"reinitvels      100\n")
                    eq1.write("for {set T 100} {$T < 300} {incr T 10} {\n")
                    eq1.write("    langevinTemp      $T;\n")
                    eq1.write("    run              1000;\n")
                    eq1.write("}\n")
                    eq1.write(f"langevinTemp    300\n")
                    eq1.write(f"run           40000;")
                    eq1.close()

                    if Probes == True:
                        eq2 = open(f"eq2.conf",'w')
                        eq2.write(f"coordinates     eq1.coor\n")
                        eq2.write(f"structure       ./{prefix}.psf\n")
                        eq2.write(f"{par_lines}\n")
                        eq2.write(f"outputname      eq2\n")
                        eq2.write(f"binaryoutput    no\n")
                        eq2.write(f"restartname     eq2\n")
                        eq2.write(f"restartfreq     2000\n")
                        eq2.write(f"binaryrestart   no\n")
                        eq2.write(f"DCDfreq         2000\n")
                        eq2.write(f"DCDfile         eq2.dcd\n")
                        eq2.write(f"outputEnergies  2000\n")
                        if lipids == True:
                            eq2.write(f"timestep        1.0\n")
                        else :
                            eq2.write(f"timestep        2.0\n")
                        eq2.write(f"fullElectFrequency 2\n")
                        eq2.write(f"nonbondedFreq   1\n")
                        eq2.write(f"rigidBonds      all\n")
                        eq2.write(f"{nonbonded}\n")
                        eq2.write(f"velocities      eq1.vel\n")
                        eq2.write(f"constraints     on\n")
                        eq2.write(f"consref         ./min.coor\n")
                        eq2.write(f"conskfile       ./{prefix}.pdb\n")
                        eq2.write(f"conskcol        B\n")
                        eq2.write(f"constraintScaling  1.0\n")
                        eq2.write(f"PME             yes\n")
                        eq2.write(f"PMEGridSpacing  1.0\n")
                        eq2.write(f"langevin        on\n")
                        eq2.write(f"langevinDamping 5\n")
                        eq2.write(f"langevinHydrogen off\n")
                        eq2.write(f"extendedSystem  eq1.xsc\n")
                        eq2.write(f"wrapWater       on\n")
                        eq2.write(f"wrapAll         on\n")
                        eq2.write("for {set T 300} {$T < 600} {incr T  10} {\n")
                        eq2.write("    langevinTemp      $T;\n")
                        eq2.write("    run              1000;\n")
                        eq2.write("}\n")
                        eq2.write(f"langevinTemp    600\n")
                        eq2.write(f"run           300000;")
                        eq2.write("for {set T 570} {$T >= 300} {incr T -30} {\n")
                        eq2.write(f"    langevinTemp     $T;\n")
                        eq2.write(f"	   run             1000;\n")
                        eq2.write("}\n")
                        eq2.close()

                        eq3 = open(f"eq3.conf",'w')
                        eq3.write(f"coordinates     eq2.coor\n")
                        eq3.write(f"structure       ./{prefix}.psf\n")
                        eq3.write(f"{par_lines}\n")
                        eq3.write(f"outputname      eq3\n")
                        eq3.write(f"binaryoutput    no\n")
                        eq3.write(f"restartname     eq3\n")
                        eq3.write(f"restartfreq     2000\n")
                        eq3.write(f"binaryrestart   no\n")
                        eq3.write(f"DCDfreq         2000\n")
                        eq3.write(f"DCDfile         eq3.dcd\n")
                        eq3.write(f"outputEnergies  2000\n")
                        if lipids == True:
                            eq3.write(f"timestep        1.0\n")
                        else :
                            eq3.write(f"timestep        2.0\n")
                        eq3.write(f"fullElectFrequency 2\n")
                        eq3.write(f"nonbondedFreq   1\n")
                        eq3.write(f"rigidBonds      all\n")
                        eq3.write(f"{nonbonded}\n")
                        eq3.write(f"velocities      eq2.vel\n")
                        eq3.write(f"PME             yes\n")
                        eq3.write(f"PMEGridSpacing  1.0\n")
                        eq3.write(f"langevin        on\n")
                        eq3.write(f"langevinTemp    300\n")
                        eq3.write(f"langevinDamping 5\n")
                        eq3.write(f"langevinHydrogen off\n")
                        eq3.write(f"useGroupPressure        yes\n")
                        eq3.write(f"useFlexibleCell         no\n")
                        eq3.write(f"useConstantArea         no\n")
                        eq3.write(f"useConstantRatio        no\n")
                        eq3.write(f"langevinPiston          on\n")
                        eq3.write(f"langevinPistonTarget    1.01325\n")  
                        eq3.write(f"langevinPistonPeriod    100.0\n")     
                        eq3.write(f"langevinPistonDecay     50.0\n")  
                        eq3.write(f"langevinPistonTemp      300.0\n")                                               
                        eq3.write(f"extendedSystem  eq2.xsc\n")
                        eq3.write(f"wrapWater       on\n")
                        eq3.write(f"wrapAll         on\n")
                        eq3.write(f"run                  300000")
                        eq3.close()

                        sim_step = int(self.sim_length.get()) * 500000

                        sim = open(f"sim.conf",'w')
                        sim.write(f"coordinates     eq3.coor\n")
                        sim.write(f"structure       ./{prefix}.psf\n")
                        sim.write(f"{par_lines}\n")
                        sim.write(f"outputname      sim\n")
                        sim.write(f"binaryoutput    no\n")
                        sim.write(f"restartname     sim\n")
                        sim.write(f"restartfreq     2000\n")
                        sim.write(f"binaryrestart   no\n")
                        sim.write(f"DCDfreq         2000\n")
                        sim.write(f"DCDfile         sim.dcd\n")
                        sim.write(f"outputEnergies  2000\n")
                        if lipids == True:
                            sim.write(f"timestep        2.0\n")
                        else :
                            sim.write(f"timestep        2.0\n")
                        sim.write(f"fullElectFrequency 2\n")
                        sim.write(f"nonbondedFreq   1\n")
                        sim.write(f"rigidBonds      all\n")
                        sim.write(f"{nonbonded}\n")
                        sim.write(f"velocities      eq3.vel\n")
                        sim.write(f"PME             yes\n")
                        sim.write(f"PMEGridSpacing  1.0\n")
                        sim.write(f"langevin        on\n")
                        sim.write(f"langevinTemp    300\n")
                        sim.write(f"langevinDamping 5\n")
                        sim.write(f"langevinHydrogen off\n")
                        sim.write(f"useGroupPressure        yes\n")
                        sim.write(f"useFlexibleCell         no\n")
                        sim.write(f"useConstantArea         no\n")
                        sim.write(f"useConstantRatio        no\n")
                        sim.write(f"langevinPiston          on\n")
                        sim.write(f"langevinPistonTarget    1.01325\n")  
                        sim.write(f"langevinPistonPeriod    100.0\n")     
                        sim.write(f"langevinPistonDecay     50.0\n")  
                        sim.write(f"langevinPistonTemp      300.0\n")                                               
                        sim.write(f"extendedSystem  eq3.xsc\n")
                        sim.write(f"wrapWater       on\n")
                        sim.write(f"wrapAll         on\n")
                        sim.write(f"run                  {sim_step}")
                        sim.close()

                    sh_file = open(f"{prefix}.sh",'w')
                    sh_file.write("namd2 +p8 min.conf > min.log\n")
                    sh_file.write("namd2 +p8 eq1.conf > eq1.log\n")
                    if Probes == True:
                        sh_file.write("namd2 +p8 eq2.conf > eq2.log\n")
                        sh_file.write("namd2 +p8 eq3.conf > eq3.log\n")
                    sh_file.write("namd2 +p8 sim.conf > sim.log\n")
                    sh_file.close()

                if Probes == True:
                    conf_files = ["min.conf", "eq1.conf", "eq2.conf", "eq3.conf", "sim.conf",f"{prefix}.sh", f"{prefix}.psf", f"{prefix}.pdb",f"{prefix}.xsc", f"{prefix}.log"]
                else :
                    conf_files = ["min.conf", "eq1.conf", "sim.conf",f"{prefix}.sh", f"{prefix}.psf", f"{prefix}.pdb",f"{prefix}.xsc", f"{prefix}.log"]

                ztrj_files = ["0chlcmpall.tcl", "1a.sh"]

                output_location = os.path.join(self.outdir_location.get(), prefix)
                con_folder = os.path.join(output_location, 'simulation_run')
                total_num_sims = self.n_sims.get() + 1
                sim_num = 1

                while sim_num < total_num_sims:
                    final_folder = f"{con_folder}{sim_num}"

                    os.makedirs(final_folder, exist_ok = True)

                    for files in conf_files:
                        shutil.copy(files, final_folder)

                    ztrj_folder = f'{final_folder}/ztrj'
                    os.makedirs(ztrj_folder, exist_ok = True)

                    for files in ztrj_files:
                        shutil.copy(files, ztrj_folder)  

                    sim_num += 1

                remove_files = ['*.pdb', '*.psf', '*.log', '*.xsc', '*.tcl', '*.txt', '*.conf', '*.sh']
                important_files = ['probe.pdb', 'probe.psf', '0chlcmpall.tcl', '1a.sh']
                for items in remove_files:
                    files = glob.glob(items)
                    for file in files:
                        if file not in important_files:
                            shutil.os.remove(file)

        prepare_button = tk.Button(mfaoo, text="Prepare System ", bd=3, command=Prepare_system)
        prepare_button.grid(row=5, column=2)
    
    def analysis_results_view(self, frame):
        mfb = tk.Frame(frame)
        mfb.grid(row=0, column=0, padx=10, pady=10)

        mfbif = tk.LabelFrame(mfb, text = "System structure, probes, and trajecotry files:", bd =2)
        mfbif.grid(row=1, column=0, padx=5, pady=5, sticky='ew')  

        def show_psf_help():
            messagebox.showinfo(
                "HELP", 
                "The PSF file of the simulated system."
            )

        # Function to browse for PSF file
        def browse_psf_file():
            tempfile = filedialog.askopenfilename(
                filetypes=[("PSF files", "*.psf"), ("All files", "*.*")]
            )
            if tempfile:  # If a file is selected
                self.protein_psf.set(tempfile)
                print(f"Selected PSF file: {self.protein_psf.get()}")  # Print the selected file for debugging
            else:
                messagebox.showinfo("No Selection", "No file was selected.")

        # Variable to store the PSF file path
        self.protein_psf = tk.StringVar()

        # Create the PSF help button
        psf_help_button = tk.Button(mfbif, text="?", padx=0, pady=0, command=show_psf_help)
        psf_help_button.grid(row=1, column=0, sticky='w') 

        psf_label = tk.Label(mfbif, text="PSF:")
        psf_label.grid(row=1, column=1, sticky='w')

        # Create and place the entry for the PSF file path
        psf_entry = tk.Entry(mfbif, width=40, textvariable=self.protein_psf)
        psf_entry.grid(row=1, column=2, sticky='ew')

        # Create the Browse button for PSF file selection
        psf_browse_button = tk.Button(mfbif, text="Browse", width=6, pady=1, command=browse_psf_file)
        psf_browse_button.grid(row=1, column=3, sticky='w')

        mfbif.grid_columnconfigure(2, weight=1)

        def show_protein_help():
            messagebox.showinfo(
                "HELP",
                "The PDB file that contains the initial coordinates of the simulation."
                )
    
        self.protein_pdb = tk.StringVar()
        protein_help_button = tk.Button(mfbif, text="?", padx=0, pady=0, command=show_protein_help)
        protein_help_button.grid(row=2, column=0, sticky='w')
    
        protein_label = tk.Label(mfbif, text="PDB:")
        protein_label.grid(row=2, column=1, sticky='w')

        pdb_entry = tk.Entry(mfbif, width=40, textvariable = self.protein_pdb)
        pdb_entry.grid(row=2, column=2, sticky='w')

        def browse_pdb_file():
            tempfile = filedialog.askopenfilename(
                filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")]
            )
            if tempfile:  # If a file is selected
                self.protein_pdb.set(tempfile)

        pdb_browse_button = tk.Button(mfbif, text="Browse", width=6, pady=1, command=browse_pdb_file)
        pdb_browse_button.grid(row=2, column=3, sticky='w')


        def probe_selection_help():
            messagebox.showinfo(
                "HELP", 
                "Selected probes used for druggability simulation. "
            )

        probe_sele_help_button = tk.Button(mfbif, text="?", padx=0, pady=0, command=probe_selection_help)
        probe_sele_help_button.grid(row=3, column=0, sticky='w')

        probe_sele_label = tk.Label(mfbif, text="Probe Selection: ")
        probe_sele_label.grid(row=3, column=1, sticky='w')

        probe_sele_listbox = tk.Entry(mfbif, width=40, textvariable = self.probe_selection)
        probe_sele_listbox.grid(row=3, column=2, sticky='w')
    
        def selection_help():
            messagebox.showinfo(
                "HELP", 
                "Selection will be used to align the protein. If protein has flexible loops or termini, they may be excluded "
                "from superimposition using this selection box. When Show button is clicked, protein as ribbon and selected "
                "atoms as spheres will be shown."
            )

        selection_help_button = tk.Button(mfbif, text="?", padx=0, pady=0, command=selection_help)
        selection_help_button.grid(row=4, column=0, sticky='w')

        selection_label = tk.Label(mfbif, text="Selection:")
        selection_label.grid(row=4, column=1, sticky='w')

        selection_entry = tk.Entry(mfbif, width=40, textvariable= self.selection)
        selection_entry.grid(row=4, column=2, sticky='w')

        def dcd_help():
            messagebox.showinfo(
                "HELP", 
                "Multiple DCD files from the same simulation or from different simulations of the same system can be specified. DCD files will be read and processed in the order they are specified here."
            )
        dcd_help_button = tk.Button(mfbif, text="?", padx=0, pady=0, command=dcd_help)
        dcd_help_button.grid(row=5, column=0, sticky='w')

        dcd_label = tk.Label(mfbif, text="DCDs:")
        dcd_label.grid(row=5, column=1, sticky='w')
    
        par_frame = tk.Frame(mfbif)
        par_frame.grid(row=5, column=2, sticky="snew")

        dcd_files = tk.StringVar()

        dcd_listbox = tk.Listbox(par_frame, 
                            activestyle='dotbox', 
                            listvariable=dcd_files, 
                            selectmode='browse', 
                            width=40, 
                            height=3, 
                            setgrid=True)
        dcd_listbox.grid(row=5, column=3, sticky="nsew")
    
        dcd_scroll = tk.Scrollbar(par_frame, command=dcd_listbox.yview)
        dcd_scroll.grid(row=5, column=4, sticky="ns")
    
        dcd_listbox.config(yscrollcommand=dcd_scroll.set)

        self.root.grid_columnconfigure(0, weight=1)
        self.root.grid_rowconfigure(0, weight=1)
        par_frame.grid_columnconfigure(0, weight=1)
        par_frame.grid_rowconfigure(0, weight=1)

        dcd_frame = tk.Frame(mfbif)
        dcd_frame.grid(row=5, column=3, sticky='w')

        system_dcds = []

        def add_dcd_files():
            tempfiles = filedialog.askopenfilenames(
                title="Select DCD files",
                filetypes=[("DCD files", "*.dcd"), ("All files", "*.*")]
            )
            if tempfiles:
                added = False
                for tempfile in tempfiles:
                    if tempfile in system_dcds:
                        messagebox.showwarning("WARNING", f"{tempfile} has already been added to the list.")
                    else:
                        system_dcds.append(tempfile)
                        added = True
                if added:
                # Update the StringVar by joining the list with newline so Listbox shows one file per line
                    dcd_files.set('\n'.join(system_dcds))

        dcd_add = tk.Button(dcd_frame, text="Add", width=6, command=add_dcd_files, padx=0, pady=0)
        dcd_add.pack(side='top')

        def remove_files():
            selected_indices = dcd_listbox.curselection()

            if not selected_indices:
                return  
            for i in reversed(selected_indices):
                del system_dcds[i]
            dcd_files.set('\n'.join(system_dcds))

        dcd_delete = tk.Button(dcd_frame, text="Remove", command=remove_files, width=6, padx=0, pady=0)
        dcd_delete.pack(side='bottom')       

        mfbto = tk.LabelFrame(mfb, text= "Grid Calculation options and parameters:", bd =2)
        mfbto.grid(row=2, column=0, padx=5, pady=5, sticky='ew') 

        def grid_help():
            messagebox.showinfo(
                "HELP", 
                "The size of a grid element, along X, Y, and Z dimensions. 0.5 works best for druggability index calculations."
            )

        grid_help_button = tk.Button(mfbto, text="?", padx=0, pady=0, command=grid_help)
        grid_help_button.grid(row=1, column=0, sticky='w')

        grid_spacing_label = tk.Label(mfbto, text="Grid resolution (A): ")
        grid_spacing_label.grid(row=1, column=1,  sticky='w')

        grid_spacing_entry = tk.Entry(mfbto, width=3, textvariable=self.grid_spacing)
        grid_spacing_entry.grid(row=1, column=2, sticky='w', padx=(0, 0)) #this works

        def condist_help():
            messagebox.showinfo(
            "HELP",
            "Only probe molecules having at least one atom within the contact distance of the protein atoms will be counted in grid calculations."
            )

        condist_help_button = tk.Button(mfbto, text="?", padx=0, pady=0, command=condist_help)
        condist_help_button.grid(row=1, column=3, sticky='w')

        condist_label = tk.Label(mfbto, text="Contact distance (A): ")
        condist_label.grid(row=1, column=4, sticky='w')

        condist_entry = tk.Entry(mfbto, width=3, textvariable=self.contact_distance)
        condist_entry.grid(row=1, column= 5, sticky='w')

        def align_help ():
           messagebox.showinfo(
                "HELP",
                "For druggability calculations, molecules will be wrapped that fall out of alignment of align snapshots."
            )

        align_help_button = tk.Button(mfbto, text="?", padx=0, pady=0, command=align_help)
        align_help_button.grid(row=1, column=6, sticky='w')

        align_label = tk.Label(mfbto, text="Align:")
        align_label.grid(row=1, column=7, sticky='w')

        align_entry = tk.Entry(mfbto, width=10, textvariable=self.align)
        align_entry.grid(row=1, column=8, sticky='w')

        mfbdo = tk.LabelFrame(mfb, text= "Output options:", bd =2)
        mfbdo.grid(row=3, column=0, padx=5, pady=5, sticky='ew')

        def outputdir_help():
            messagebox.showinfo(
                "HELP",
                "Output folder, default is current working directory."
            )

        outputdir_help_button = tk.Button(mfbdo, text="?", padx=0, pady=0, command=outputdir_help)
        outputdir_help_button.grid(row=0, column=0, sticky='w')

        outputdir_label = tk.Label(mfbdo, text="Output folder:     ",)
        outputdir_label.grid(row=0, column=1, sticky='w')

        def browse_output_folder():
            folder_path = filedialog.askdirectory(
                title="Select Output Folder"  # Title for the directory dialog
            )
            if folder_path:  # If a folder is selected
                self.outputdir_location.set(folder_path)

        outputdir_entry = tk.Entry(mfbdo, width=40, textvariable=self.outputdir_location)
        outputdir_entry.grid(row=0, column=2, sticky='w')

        outputdir_browse_button = tk.Button(mfbdo, text="Browse", command=browse_output_folder)
        outputdir_browse_button.grid(row=0, column=3, sticky='ew')

        def prefix_help():
            messagebox.showinfo(
                "HELP",
                "All output files and folders will start with this prefix. "
                "A unique and descriptive prefix choice may allow running multiple simulations in the same folder."
            )

        prefix_help_button = tk.Button(mfbdo, text="?", padx=0, pady=0, command=prefix_help)
        prefix_help_button.grid(row=1, column=0, sticky='w')

        prefix_label = tk.Label(mfbdo, text="Output prefix:")
        prefix_label.grid(row=1, column=1, sticky='w')

        prefix_path = tk.Entry(mfbdo, width=15, textvariable=self.prefix)
        prefix_path.grid(row=1, column=2, sticky='w')

        mfbfo = tk.LabelFrame(mfb, text= "Druggability Analysis:", bd =2)
        mfbfo.grid(row=4, column=0, padx=5, pady=5, sticky='ew') 

        def temp_help():
            messagebox.showinfo(
                "Help",
                "Temperature of the system in the productive simulation (NAMD configuration files generated by this GUI sets temperature to 300 K)."
            )

        temp_help_button = tk.Button(mfbfo, text="?", padx=0, pady=0, command=temp_help)
        temp_help_button.grid(row=0, column=0, sticky='w',padx=(40, 0))

        temp_label = tk.Label(mfbfo, text="Temperature (K):")
        temp_label.grid(row=0, column=1, sticky='w')

        temp_entry = tk.Entry(mfbfo, width=3, textvariable=self.temperature)
        temp_entry.grid(row=0, column= 2, sticky='w', padx=(0, 60))

        def radius_help():
            messagebox.showinfo(
                "HELP",
                "Probe merge radius in angstroms. Twice the size of largest effective probe radius gives reasonable solutions."
            )

        radius_help_button = tk.Button(mfbfo, text="?", padx=0, pady=0, command=radius_help)
        radius_help_button.grid(row=0, column=3, sticky='w')

        radius_label = tk.Label(mfbfo, text="Probe merge radius (A):")
        radius_label.grid(row=0, column=4, sticky='w')

        radius_entry = tk.Entry(mfbfo, width=3, textvariable=self.dia_merge_radius)
        radius_entry.grid(row=0, column=5, sticky='w')

        def frame_num_help():
            messagebox.showinfo(
                "HELP",
                "Number of frames used in determining the grid data. Volmap by default uses one frame to calculate frame averagte grid data."
            )

        frame_num_help_button = tk.Button(mfbfo, text="?", padx=0, pady=0, command=frame_num_help)
        frame_num_help_button.grid(row=1, column=0, sticky='w',padx=(40, 0))


        frame_label = tk.Label(mfbfo, text="Number of frames:")
        frame_label.grid(row=1, column=1, sticky="w")

        frame_entry = tk.Entry(mfbfo, width=3, textvariable=self.dia_n_frames)
        frame_entry.grid(row=1, column=2, sticky='w', padx=(0, 50))

        def numprobes_help():
            messagebox.showinfo(
                "HELP",
                "Number of probe binding hotspots to merge to make a drug-size solution."
            )

        numprobes_help_button = tk.Button(mfbfo, text="?", padx=0, pady=0, command=numprobes_help)
        numprobes_help_button.grid(row=1, column=3, sticky='w')

        numprobes_label = tk.Label(mfbfo, text="Number of hotspots to merge:")
        numprobes_label.grid(row=1, column=4, sticky='w')

        numprobes_entry = tk.Entry(mfbfo, width=3, textvariable=self.dia_n_probes)
        numprobes_entry.grid(row=1, column=5, sticky='w')

        def deltag_help():
            messagebox.showinfo(
                "HELP",
                "Probe binding free energy to determine binding hotspots."
            )

        deltag_help_button = tk.Button(mfbfo, text="?", padx=0, pady=0, command=deltag_help)
        deltag_help_button.grid(row=2, column=0, sticky='w',padx=(40, 0))

        deltag_label = tk.Label(mfbfo, text="Hotspot dG (kcal/mol):")
        deltag_label.grid(row=2, column=1, sticky='w')

        deltag_entry = tk.Entry(mfbfo, width=3, textvariable=self.dia_delta_g)
        deltag_entry.grid(row=2, column=2, sticky='w', padx=(0, 50))

        def minnpro_help():
            messagebox.showinfo(
                "HELP",
                "Minimum number of hotspots in an acceptable drug-size solution."
            )

        minnpro_help_button = tk.Button(mfbfo, text="?", padx=0, pady=0, command=minnpro_help)
        minnpro_help_button.grid(row=2, column=3, sticky='w')

        minnpro_label = tk.Label(mfbfo, text="Minimum number of hotspots:")
        minnpro_label.grid(row=2, column=4, sticky='w')

        minnpro_entry = tk.Entry(mfbfo, width=3, textvariable=self.dia_min_n_probes)
        minnpro_entry.grid(row=2, column=5, sticky='w')

        def affinity_help():
            messagebox.showinfo(
                "HELP",
                "Lowest affinity to report a solution in micromolar units."
            )

        affinity_help_button = tk.Button(mfbfo, text="?", padx=0, pady=0, command=affinity_help)
        affinity_help_button.grid(row=3, column=0, sticky='w',padx=(40, 0))

        affinity_label = tk.Label(mfbfo, text="Lowest affinity (uM):")
        affinity_label.grid(row=3, column=1, sticky="w")

        affinity_entry = tk.Entry(mfbfo, width=3, textvariable=self.dia_low_affinity)
        affinity_entry.grid(row=3, column=2, sticky='w', padx=(0, 50))

        def charge_help():
            messagebox.showinfo(
                "HELP",
                "Maximum absolute charge to accept solutions."
            )

        charge_help_button = tk.Button(mfbfo, text="?", padx=0, pady=0, command=charge_help)
        charge_help_button.grid(row=3, column=3, sticky='w')

        charge_label = tk.Label(mfbfo, text="Maximum absolute charge (e):")
        charge_label.grid(row=3, column=4, sticky='w')

        charge_entry = tk.Entry(mfbfo, width=3, textvariable=self.dia_max_charge)
        charge_entry.grid(row=3, column=5, sticky='w')

        def solution_help():
            messagebox.showinfo(
                "HELP",
                "Number of solutions to report in each distinct potential binding site."
            )

        solution_help_button = tk.Button(mfbfo, text="?", padx=0, pady=0, command=solution_help)
        solution_help_button.grid(row=4, column=0, sticky='w',padx=(40, 0))

        solution_label = tk.Label(mfbfo, text="Number of solutions:")
        solution_label.grid(row=4, column=1, sticky='w')

        solution_entry = tk.Entry(mfbfo, width=3, textvariable=self.dia_n_solutions)
        solution_entry.grid(row=4, column=2, sticky='w', padx=(0, 50))

        def ncharged_help():
            messagebox.showinfo(
                "HELP",
                "Maximum number of charged hotspots in a solution."
            )

        ncharged_help_button = tk.Button(mfbfo, text="?", padx=0, pady=0, command=ncharged_help)
        ncharged_help_button.grid(row=4, column=3, sticky='w')

        ncharged_label = tk.Label(mfbfo, text="Number of charged hotspots: ")
        ncharged_label.grid(row=4, column=4, sticky='w')

        ncharged_entry = tk.Entry(mfbfo, width=3, textvariable=self.dia_n_charged)
        ncharged_entry.grid(row=4, column=5, sticky='w')

        mfbgo = tk.LabelFrame(mfb, bd =2)
        mfbgo.grid(row=5, column=0, padx=5, pady=5, sticky='ew') 

        def analyze_system():

            verbose = 'info'

            def buildGrids(prefix, pdb_psf_dcds, probes, align=self.align.get(), protein=self.selection.get(), contacti=self.contact_distance.get(), resolution=self.grid_spacing.get(), savedcd=False):
                
                if len(pdb_psf_dcds) > 1:
                    reference = parsePDB(pdb_psf_dcds[0][0]).select(align).copy()
                else:
                    reference = None

                UNITCELL = []
                DCDOUT = {}
                nframe = 0

                for p in probes:
                    DCDOUT[p] = DCDFile(prefix + '_' + p + '.dcd', 'w')

                from prody import startLogfile
                startLogfile(prefix + '_grid.log')

                pdb, psf = pdb_psf_dcds[0][:2]
                dcds = pdb_psf_dcds[0][2:]

                pdb = parsePDB(pdb, AtomGroup=parsePSF(psf), long_resname = True)
                palign = pdb.select(align)
                if reference is not None:
                    matchAlign(palign, reference)
                pcenter = calcCenter(palign)

                # make sure all probe names select some residues
                from prody import plog
                probe_selstr = 'noh and resname'
                for p in probes:
                    sel = pdb.select('noh and resname ' + p)
                    if sel is None:
                        continue
                        #raise ValueError('probe ' + p + ' is not found in the system')
                    hv = sel.getHierView()
                    n = hv.numResidues()
                    res = next(hv.iterResidues())
                    writePDB(prefix + '_' + p + '.pdb', res, long_resname = True)
                    plog(str(n) + ' copies of ' + p + ' is found.')
                    probe_selstr += ' ' + p

                # start trajectories for reading
                plog('Opening DCD files for reading.')
                dcd = Trajectory(dcds[0]) 
                for fn in dcds[1:]:
                    dcd.addFile(fn)
                plog(str(len(dcd)) + ' frames from ' + str(len(dcds)) +  
                    ' file(s) will be evaluated.')

                dcd.link(pdb)
                # make alignment selection, calling `frame.superpose` will align frames 
                # based on this selection
                dcd.setAtoms(palign)

                # make a probe selection
                PRBSEL = pdb.select('noh and ' + probe_selstr)
                PRBIDX = PRBSEL.getIndices()
                # make a copy of probe selection, it will be faster to make contact search
                # in this copy after copying the coordinates
                PROBES = PRBSEL.copy().toAtomGroup()

                pcontact = pdb.select('protein')
                writePDB(prefix + '_protein_heavyatoms.pdb', pcontact)

                LOGGER.progress('Evaluating frames:', len(dcd))
                if savedcd:
                    dcdout = DCDFile(pdb.getTitle() + '_aligned_wrapped.dcd', 'w')

                for i, frame in enumerate(dcd):
                    #if i % 20 != 0:
                    #    continue
                    # first center the frame
                    moveAtoms(palign, to=pcenter, ag=True)

                    # wrap probes that are out of the box
                    unitcell = frame.getUnitcell()[:3]
                    UNITCELL.append(unitcell)
                    coords = pdb._coords[0]
                    coords[PRBIDX] = wrapAtoms(coords[PRBIDX], unitcell, pcenter)  

                    # superpose frame coordinates onto selected atoms
                    frame.superpose()

                    if savedcd:
                        dcdout.write(coords)

                    # identify probes in contact and write them in DCD file
                    PROBES.setCoords(PRBSEL.getCoords())
                    cont = PROBES.select(f'same residue as within {contacti} of pcontact', 
                                         pcontact=pcontact)

                    if cont:
                        for res in cont.getHierView().iterResidues():
                            DCDOUT[res.getResname()].write(res._getCoords()) 
                    nframe += 1
                    LOGGER.update(i)

                dcd.close()

                for p in probes: DCDOUT.pop(p).close()

                UNITCELL = array(UNITCELL).max(0)

                offset = pcenter - UNITCELL / 2.
                length = UNITCELL
                n_bins = ceil(length / resolution).astype(int)
                bins = [arange(n) * resolution + offset[i] for i, n in enumerate(n_bins)]
                probe_grids = []
                for p in probes:    
                    fn = prefix + '_' + p
                    if not os.path.getsize(fn + '.dcd'):
                        plog('No ' + p + ' found in contact with the protein.')
                        continue
                    e = parseDCD(fn + '.dcd')
                    c = calcCenter(e._getCoordsets())
                    garray = histogramdd(c, bins)
                    grid = OpenDX()
                    grid.filename = fn + '.dx'
                    grid.name = fn
                    grid.spacing = array([resolution, resolution, resolution])
                    grid._origin = offset
                    grid.offset = offset
                    grid.array = garray[0] / nframe
                    grid.shape = grid.array.shape
                    grid._comments = ['# comment', 'object']
                    grid.write(fn + '.dx')
                    probe_grids.append((p, fn + '.dx'))

                from prody import closeLogfile
                closeLogfile(prefix + '_grid.log')
                return probe_grids

            def calcDruggability(prefix, probe_grids, **kwargs):


                dia = druggability.DIA(prefix, workdir=self.outputdir_location.get(), verbose=verbose)
                # check parameters for their correctness
                dia.set_parameters(temperature=kwargs.get('temperature', self.temperature.get())) # K (productive simulation temperature)
                dia.set_parameters(delta_g=kwargs.get('delta_g', self.dia_delta_g.get())) # kcal/mol (probe binding hotspots with lower values will be evaluated)
                dia.set_parameters(n_probes=kwargs.get('n_probes', self.dia_n_probes.get())) # (number of probes to be merged to determine achievable affinity of a potential site)
                dia.set_parameters(min_n_probes=kwargs.get('min_n_probes', self.dia_min_n_probes.get())) # (minimum number of probes to be merged for an acceptable soltuion)
                dia.set_parameters(merge_radius=kwargs.get('merge_radius', self.dia_merge_radius.get())) # A (distance within which two probes will be merged)
                dia.set_parameters(low_affinity=kwargs.get('low_affinity', self.dia_low_affinity.get())) # microMolar (potential sites with affinity better than this value will be reported)
                dia.set_parameters(n_solutions=kwargs.get('n_solutions', self.dia_n_solutions.get())) # (number of drug-size solutions to report for each potential binding site)
                dia.set_parameters(max_charge=kwargs.get('max_charge', self.dia_max_charge.get())) # (maximum absolute total charge, where total charge is occupancy weighted sum of charges on probes)
                dia.set_parameters(n_charged=kwargs.get('n_charged', self.dia_n_charged.get())) # (maximum number of charged hotspots in a solution)
                dia.set_parameters(n_frames=kwargs.get('n_frames', self.dia_n_frames.get())) # number of frames (if volmap was used (.dx), 1 is correct)

                for probe_type, grid_file in probe_grids: dia.add_probe(probe_type, grid_file)

                # Do the calculations
                dia.perform_analysis()
                dia.pickle()
                # Evaluate a ligand. Be sure that the ligand bound structure is superimposed
                # onto PROTEIN_heavyatoms.pdb
                ligand = kwargs.get('ligand', None)
                if ligand:
                    dia.evaluate_ligand(ligand)

            #def evalLigandSite(prefix, ligand, radius=1.5, delta_g=-0.5):

                #dia = pickler(os.path.join(prefix, prefix + '.dso.gz'))
                #dia.evaluate_ligand(ligand, radius=radius, delta_g=delta_g)

            
            # output file and folder names will start with the following
            prefix_n = self.prefix.get()
            os.chdir(self.outputdir_location.get())

            dcds_raw = dcd_files.get()

            try:
                dcds_list = ast.literal_eval(dcds_raw)
                if not isinstance(dcds_list, (list, tuple)):
                    raise ValueError("Parsed object is not a list or tuple")
            except Exception:
                dcds_list = dcds_raw.split()

            probes = list(ast.literal_eval(self.probe_selection.get()))
            pdb_psf_dcds = [[
            self.protein_pdb.get(),
            self.protein_psf.get(),
            *dcds_list
            ]]

            # build grids
            buildGrids(prefix_n, pdb_psf_dcds, probes, savedcd= True)

            # DRUGGABILITY
            parameters = {}
            parameters['temperature'] = self.temperature.get()    # K (productive simulation temperature)
            parameters['delta_g'] = self.dia_delta_g.get()        # kcal/mol (probe binding hotspots with lower values will be evaluated)
            parameters['n_probes'] = self.dia_n_probes.get()         # (number of probes to be merged to determine achievable affinity of a potential site)
            parameters['min_n_probes'] = self.dia_min_n_probes.get()     # (minimum number of probes to be merged for an acceptable soltuion)
            parameters['merge_radius'] = self.dia_merge_radius.get()    # A (distance within which two probes will be merged)
            parameters['low_affinity'] = self.dia_low_affinity.get()     # microMolar (potential sites with affinity better than this value will be reported)
            parameters['n_solutions'] = self.dia_n_solutions.get()       # (number of drug-size solutions to report for each potential binding site)
            parameters['max_charge'] = self.dia_max_charge.get()        # (maximum absolute total charge, where total charge is occupancy weighted sum of charges on probes)
            parameters['n_charged'] = self.dia_n_charged.get()         # (maximum number of charged hotspots in a solution)
            parameters['n_frames'] = self.dia_n_frames.get()          # number of frames (if volmap was used (.dx), 1 is correct)

            # Probe grid files are automatically determined based on prefix
            from glob import glob
            probe_grids = [(os.path.splitext(fn)[0].split('_')[-1], fn) 
                            for fn in glob(prefix_n + '_*.dx')]  

            calcDruggability(prefix_n, probe_grids, **parameters)

            # LIGAND SITE
            # Evaluate a ligand. Be sure that the ligand bound structure is superimposed
            # onto PROTEIN_heavyatoms.pdb
            #evalLigandSite(prefix, 'ligand.pdb', radius=1.5, delta_g=-0.5)
                

        analyze_system_button = tk.Button(mfbgo, text="Access Druggability", bd=3, command=analyze_system)
        analyze_system_button.grid(row=0, column=0, padx=(225,0))



def drugui_data(self):


     #global PACKAGEPATH
    global PROBEDATA
    global PROBETYPES
    global PROBETOPPAR


    probedata = f"""
    set PROBEDATA [dict create]
    set PROBETYPES [dict create core "Core probes"]
    dict set PROBETYPES polar "Polar probes"
    dict set PROBETYPES hydrophobe "Hydrophobes"
    dict set PROBETYPES negative "Negatively charged"
    dict set PROBETYPES positive "Positively charged"
    dict set PROBETYPES ring5 "5-membered rings"
    dict set PROBETYPES ring6 "6-membered rings"
    set PROBETOPPAR [dict create PBDA "probe2.top probe.prm"]
    dict set PROBETOPPAR CGenFF "top_all36_cgenff.rtf par_all36_cgenff.prm"
    set PACKAGEPATH {self.Druggability_path}
    """
    probedata += """
    foreach {key} [dict keys $PROBETYPES] {
    dict set PROBEDATA $key ""
    }

    set inf [open [file join $PACKAGEPATH "probeV2.dat"] r]
    foreach line [split [read -nonewline $inf] \\n] {
        if {[llength $line] < 3 || [string equal -length 1 $line "#"]} {continue}
        set resi [lindex $line 0]
        set key [lindex $line 1]
        if {![dict exists $PROBEDATA $resi]} {
            dict set PROBEDATA $resi default 0
            dict set PROBEDATA $resi alias ""
            dict set PROBEDATA $resi atomnames ""
            dict set PROBEDATA $resi charge 0
            dict set PROBEDATA $resi source ""
        }
        if {$key == "default" && [lindex $line 2] > 0} {
            dict lappend PROBEDATA defaults $resi
        }
        if {$key == "type" && [lindex $line 2] > 0} {
            dict lappend PROBEDATA [lindex $line 2] $resi
        }
        dict set PROBEDATA $resi $key [lrange $line 2 end]
        if {$key == "alias" && [expr [llength [dict get $PROBEDATA $resi $key]] % 2] != 0} {
            error "Problem in aliases of $resi!"
        }
     }
    close $inf

    set ipronames "C2 H21 C1 H11 H12 H13 C3 C4 C5 C6 H31 H32 H33 OH2 HO2"
    dict for {src toppar} $PROBETOPPAR {
        set inf [open [file join $PACKAGEPATH [lindex $toppar 0]] r]
        puts "$inf"
        set resi "____"
        set prev "____"
        foreach line [split [read -nonewline $inf] \\n] {
            if {[dict exists $PROBEDATA $resi] && [string range $line 0 3] == "ATOM"} {
                dict set PROBEDATA $resi atomnames "[dict get $PROBEDATA $resi atomnames] [lindex $line 1]"
            } elseif {[string range $line 0 3] == "RESI"} {
              set prev $resi
              set resi [string trim [lindex $line 1]]
              if {[dict exists $PROBEDATA $resi]} {
                dict set PROBEDATA $resi charge [lindex $line 2]
                dict set PROBEDATA $resi source $src
              }
              if {[dict exists $PROBEDATA $prev]} {
                set thisnames [dict get $PROBEDATA $prev atomnames]
                set remove 0
                set aliases [dict get $PROBEDATA $prev alias]
                foreach {name alias} $aliases {
                    if {[lsearch $ipronames $name] == -1} {
                        error "Invalid $prev alias, $name is not found in IPRO!"
                        #set remove 1
                    } elseif {[lsearch $thisnames $alias] == -1} {
                      error "Invalid $prev alias, $alias is not found in $prev!"
                    #set remove 1
                    }
                    #if {$remove} {
                    #  set PROBEDATA [dict remove $PROBEDATA $prev]
                    #  break
                    #}
                }
            }
        }
     }
    close $inf
    }

    set probe_data [open "probedata.txt" w]
    puts $probe_data $PROBEDATA
    close $probe_data

    set probetypes [open "probetypes.txt" w]
    puts $probetypes $PROBETYPES
    close $probetypes

    set probetoppar [open "probetoppar.txt" w]
    puts $probetoppar $PROBETOPPAR
    close $probetoppar
    exit
    """
    with open('probedata.tcl', 'w') as tcl_file:
        tcl_file.write(probedata)
    vmd = self.tk_vmd_executable.get()
    subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'probedata.tcl'])
    with open('probedata.txt', 'r') as output_file:
            PROBEDATA = output_file.read().strip()

