# import libraries
from textwrap import fill

import numpy as np
import pandas as pd
from fontTools.cffLib import width
from pandastable import Table, config
import random
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from Bio import PDB
from Bio.PDB import PPBuilder
import math
from tkinter import filedialog
from tkinter.filedialog import askdirectory
from tkinter import *
import os

# sets the code classes and their methods

class HAMP:
    # class HAMP that analyses the protein file

    def _set_degrees(self):
        # Load the structure
        parser = PDB.PDBParser()
        structure = parser.get_structure("prot", self.path)
        # Build polypeptides to handle chain breaks
        ppb = PPBuilder()
        data = {'ResNum': [], 'Residue': [], 'PHI': [], 'PSI': []}
        for model in structure:
            for chain in model:
                polypeptides = ppb.build_peptides(chain)
                for poly in polypeptides:
                    # get_phi_psi_list() returns a list of (phi, psi) tuples in radians
                    phi_psi_list = poly.get_phi_psi_list()

                    for res, (phi, psi) in zip(poly, phi_psi_list):
                        # Convert radians to degrees if needed
                        phi_deg = math.degrees(phi) if phi else 0
                        psi_deg = math.degrees(psi) if psi else 0
                        data['ResNum'].append(res.id[1])
                        data['Residue'].append(res.resname)
                        data['PHI'].append(phi_deg)
                        data['PSI'].append(psi_deg)
                        # print(f"Residue: {res.id[1]} | {res.resname} | {round(phi_deg, 1)} | {round(psi_deg, 1)}")
        self.__df1 = pd.DataFrame(data)

    def _set_scales(self):
        # set the scales used in the method
        scale = {
            'Residue': ['PHE', 'MET', 'ILE', 'LEU', 'TRP', 'VAL', 'CYS', 'TYR', 'ALA', 'HIS', 'GLY', 'THR', 'PRO',
                        'ARG','SER', 'GLN', 'GLU', 'ASN', 'ASP', 'LYS'],
            'ResCode': ['F', 'M', 'I', 'L', 'W', 'V', 'C', 'Y', 'A', 'H', 'G', 'T', 'P', 'R', 'S', 'Q', 'E', 'N', 'D',
                        'K'],
            'ParkerScale': [-9.20, -4.20, -8.00, -9.20, -10.00, -3.70, 1.40, -1.90, 2.10, 2.10, 5.70, 5.20, 2.10, 4.20,
                            6.50, 6.00, 7.80, 7.00, 10.00, 5.70],
            'YoungScale': [5.12, 4.91, 4.88, 4.65, 4.36, 4.17, 4.00, 3.24, 2.82, 2.75, 2.34, 2.30, 2.22, 2.18, 2.07,
                           1.98, 1.94, 1.90, 1.81, 1.50],
            'ChargeScale': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, -1, 1],
            'LinScale': [5, 4, 3, 3, 11, 3, 1, 10, 2, 19, 1, 11, 7, 20, 6, 14, 18, 6, 9, 44]}
        scales = pd.DataFrame(scale)
        self.__df1 = self.__df1.merge(scales, on='Residue', how='left')

    def _set_delta_degrees(self):
        self.__df1['DeltaPHI'] = 0.0
        self.__df1['DeltaPSI'] = 0.0
        for row in self.__df1.index.tolist():
            if row == 0:
                #the first residue's delta vectors are set as 0(0.0, 0.0)
                continue
            else:
                # calcs delta angles on every other lines, general case
                # DeltaPHI
                if self.__df1.at[row - 1, 'DeltaPHI'] - (
                        ((self.__df1.at[row - 1, 'PSI'] - 0.1) / abs(self.__df1.at[row - 1, 'PSI'] - 0.1)) *
                        self.__df1.at[row, 'PHI']) > 180:
                    self.__df1.at[row, 'DeltaPHI'] = self.__df1.at[row - 1, 'DeltaPHI'] - (
                            ((self.__df1.at[row - 1, 'PSI'] - 0.1) / abs(self.__df1.at[row - 1, 'PSI'] - 0.1)) *
                            self.__df1.at[row, 'PHI']) - 360
                elif self.__df1.at[row - 1, 'DeltaPHI'] - (
                        ((self.__df1.at[row - 1, 'PSI'] - 0.1) / abs(self.__df1.at[row - 1, 'PSI'] - 0.1)) *
                        self.__df1.at[row, 'PHI']) < -180:
                    self.__df1.at[row, 'DeltaPHI'] = self.__df1.at[row - 1, 'DeltaPHI'] - (
                            ((self.__df1.at[row - 1, 'PSI'] - 0.1) / abs(self.__df1.at[row - 1, 'PSI'] - 0.1)) *
                            self.__df1.at[row, 'PHI']) + 360
                else:
                    self.__df1.at[row, 'DeltaPHI'] = self.__df1.at[row - 1, 'DeltaPHI'] - (
                            ((self.__df1.at[row - 1, 'PSI'] - 0.1) / abs(self.__df1.at[row - 1, 'PSI'] - 0.1)) *
                            self.__df1.at[row, 'PHI'])

                # DeltaPSI
                if self.__df1.at[row - 1, 'DeltaPSI'] - (
                        ((self.__df1.at[row - 1, 'PHI'] - 0.1) / abs(self.__df1.at[row - 1, 'PHI'] - 0.1)) *
                        self.__df1.at[row - 1, 'PSI']) > 180:
                    self.__df1.at[row, 'DeltaPSI'] = self.__df1.at[row - 1, 'DeltaPSI'] - (
                            ((self.__df1.at[row - 1, 'PHI'] - 0.1) / abs(self.__df1.at[row - 1, 'PHI'] - 0.1)) *
                            self.__df1.at[row - 1, 'PSI']) - 360
                elif self.__df1.at[row - 1, 'DeltaPSI'] - (
                        ((self.__df1.at[row - 1, 'PHI'] - 0.1) / abs(self.__df1.at[row - 1, 'PHI'] - 0.1)) *
                        self.__df1.at[row - 1, 'PSI']) < -180:
                    self.__df1.at[row, 'DeltaPSI'] = self.__df1.at[row - 1, 'DeltaPSI'] - (
                            ((self.__df1.at[row - 1, 'PHI'] - 0.1) / abs(self.__df1.at[row - 1, 'PHI'] - 0.1)) *
                            self.__df1.at[row - 1, 'PSI']) + 360
                else:
                    self.__df1.at[row, 'DeltaPSI'] = self.__df1.at[row - 1, 'DeltaPSI'] - (
                            ((self.__df1.at[row - 1, 'PHI'] - 0.1) / abs(self.__df1.at[row - 1, 'PHI'] - 0.1)) *
                            self.__df1.at[row - 1, 'PSI'])

    def _set_vectors(self):
        # calc the vector's components for Young, Park and Charge scales
        self.__df1['YPHI'] = 0.0
        self.__df1['YPSI'] = 0.0
        self.__df1['PPHI'] = 0.0
        self.__df1['PPSI'] = 0.0
        for row in self.__df1.index.tolist():
            # calcs Young vector's components
            if row == 0:
                # the components on the first residue are set as (0.0, 0.0)
                pass
            else:
                self.__df1.at[row, 'YPHI'] = np.sign(self.__df1.at[row, 'YoungScale']) * self.__df1.at[
                row, 'DeltaPHI'] * np.sqrt(
                (self.__df1.at[row, 'YoungScale'] ** 2) / (
                            (self.__df1.at[row, 'DeltaPHI'] ** 2) + (self.__df1.at[row, 'DeltaPSI'] ** 2)))
                self.__df1.at[row, 'YPSI'] = np.sign(self.__df1.at[row, 'YoungScale']) * self.__df1.at[
                row, 'DeltaPSI'] * np.sqrt(
                (self.__df1.at[row, 'YoungScale'] ** 2) / (
                            (self.__df1.at[row, 'DeltaPHI'] ** 2) + (self.__df1.at[row, 'DeltaPSI'] ** 2)))
                # calcs Parker vector's components
                self.__df1.at[row, 'PPHI'] = np.sign(self.__df1.at[row, 'ParkerScale']) * self.__df1.at[
                row, 'DeltaPHI'] * np.sqrt(
                (self.__df1.at[row, 'ParkerScale'] ** 2) / (
                            (self.__df1.at[row, 'DeltaPHI'] ** 2) + (self.__df1.at[row, 'DeltaPSI'] ** 2)))
                self.__df1.at[row, 'PPSI'] = np.sign(self.__df1.at[row, 'ParkerScale']) * self.__df1.at[
                row, 'DeltaPSI'] * np.sqrt(
                (self.__df1.at[row, 'ParkerScale'] ** 2) / (
                            (self.__df1.at[row, 'DeltaPHI'] ** 2) + (self.__df1.at[row, 'DeltaPSI'] ** 2)))

    def _set_quadrant_points(self):
        # calc the quadrant points for every residue in the peptide chain based on their delta degrees
        self.__df1['QP'] = 0
        for row in self.__df1.index.tolist():
            if self.__df1.at[row, 'DeltaPHI'] >= 0 and self.__df1.at[row, 'DeltaPSI'] >= 0:
                self.__df1.at[row, 'QP'] = 1
            elif self.__df1.at[row, 'DeltaPHI'] < 0 and self.__df1.at[row, 'DeltaPSI'] >= 0:
                self.__df1.at[row, 'QP'] = 2
            elif self.__df1.at[row, 'DeltaPHI'] < 0 and self.__df1.at[row, 'DeltaPSI'] < 0:
                self.__df1.at[row, 'QP'] = 3
            else:
                self.__df1.at[row, 'QP'] = 4

    def _set_random_fragments(self):
        # set the random fragments for the sequence
        data = {'PI': [], 'PF': [], 'PM': []}
        for i in range(int((43 / 2) * ((self.__df1['ResNum'].max() - 42) + (self.__df1['ResNum'].max() - 7)))):
            data['PI'].append(random.randint(1, self.__df1['ResNum'].max() - 7))
            if self.__df1['ResNum'].max() - data['PI'][i] >= 42:
                data['PF'].append(data['PI'][i] + random.randint(7, 42))
            else:
                data['PF'].append(data['PI'][i] + random.randint(7, self.__df1['ResNum'].max() - data['PI'][i]))
            data['PM'].append(int((data['PI'][i] + data['PF'][i]) / 2))
        self.__df2 = pd.DataFrame(data)

    def _set_fragment_vectors(self):
        # calcs the resultant vectors for each random fragment generated
        self.__df2['VY'] = 0
        self.__df2['VP'] = 0
        self.__df2['K'] = ''
        for row in self.__df2.index.tolist():
            # calcs the Young vector
            sYPHI = self.__df1.loc[(self.__df1['ResNum'] >= self.__df2.at[row, 'PI']) & (
                    self.__df1['ResNum'] <= self.__df2.at[row, 'PF'])]['YPHI'].sum()
            sYPSI = self.__df1.loc[(self.__df1['ResNum'] >= self.__df2.at[row, 'PI']) & (
                    self.__df1['ResNum'] <= self.__df2.at[row, 'PF'])]['YPSI'].sum()
            if sYPHI >= 0 and sYPSI >= 0:
                self.__df2.at[row, 'VY'] = 1
            elif sYPHI < 0 and sYPSI >= 0:
                self.__df2.at[row, 'VY'] = 2
            elif sYPHI < 0 and sYPSI < 0:
                self.__df2.at[row, 'VY'] = 3
            else:
                self.__df2.at[row, 'VY'] = 4

            # calcs the Park vector
            sPPHI = self.__df1.loc[(self.__df1['ResNum'] >= self.__df2.at[row, 'PI']) & (
                    self.__df1['ResNum'] <= self.__df2.at[row, 'PF'])]['PPHI'].sum()
            sPPSI = self.__df1.loc[(self.__df1['ResNum'] >= self.__df2.at[row, 'PI']) & (
                    self.__df1['ResNum'] <= self.__df2.at[row, 'PF'])]['PPSI'].sum()
            if sPPHI >= 0 and sPPSI >= 0:
                self.__df2.at[row, 'VP'] = 1
            elif sPPHI < 0 and sPPSI >= 0:
                self.__df2.at[row, 'VP'] = 2
            elif sPPHI < 0 and sPPSI < 0:
                self.__df2.at[row, 'VP'] = 3
            else:
                self.__df2.at[row, 'VP'] = 4

            # sets the key label
            self.__df2.at[row, 'K'] = f"{self.__df2.at[row, 'VP']}{self.__df2.at[row, 'VY']}"

    def _set_hit(self):
        # defines if the fragment meet the hit conditions
        self.__df2['Hit'] = 0.0000

        for row in self.__df2.index.tolist():

            # check conditions for the charged face
            c1 = self.__df1[(self.__df1['ResNum'] >= self.__df2.at[row, 'PI']) & (
                        self.__df1['ResNum'] <= self.__df2.at[row, 'PF']) & (
                                        self.__df1['QP'] == self.__df2.at[row, 'VP'])]['ChargeScale'].sum() < 1
            c2 = sum(self.__df1.loc[(self.__df1['ResNum'] >= self.__df2.at[row, 'PI']) & (self.__df1['ResNum'] <= self.__df2.at[row, 'PF'])]['QP'] == self.__df2.at[row, 'VP']) < int(
                (self.__df2.at[row, 'PF'] - self.__df2.at[row, 'PI'] + 1) / 3)

            # check conditions for the apolar face
            c3 = self.__df1[(self.__df1['ResNum'] >= self.__df2.at[row, 'PI']) & (
                        self.__df1['ResNum'] <= self.__df2.at[row, 'PF']) & (
                                        self.__df1['QP'] == self.__df2.at[row, 'VY'])]['ChargeScale'].sum() < 0
            c4 = sum(self.__df1.loc[(self.__df1['ResNum'] >= self.__df2.at[row, 'PI']) & (self.__df1['ResNum'] <= self.__df2.at[row, 'PF'])]['QP'] == self.__df2.at[row, 'VY']) < int(
                (self.__df2.at[row, 'PF'] - self.__df2.at[row, 'PI'] + 1) / 4)

            # check if the faces dont coincide
            c5 = self.__df2.at[row, 'VY'] == self.__df2.at[row, 'VP']

            # check the net charge of the fragment
            c6 = self.__df1[(self.__df1['ResNum'] >= self.__df2.at[row, 'PI']) & (
                        self.__df1['ResNum'] <= self.__df2.at[row, 'PF'])]['ChargeScale'].sum() <= 0

            if (c1 or c2 or c3 or c4 or c5 or c6):
                self.__df2.at[row, 'Hit'] = 0.0000
            else:
                self.__df2.at[row, 'Hit'] = 0.9999

    def _calc_prob(self):
        # calcs the probability of hit for each fragment
        self.__df2['P'] = 0.0
        for row in self.__df2.index.tolist():
            if self.__df2.at[row, 'Hit'] == 0.9999:
                self.__df2.at[row, 'P'] = self.__df2[
                    (self.__df2['PI'] < self.__df2.at[row, 'PM']) & (self.__df2['PF'] > self.__df2.at[row, 'PM']) & (
                                self.__df2['K'] == self.__df2.at[row, 'K'])][
                    'Hit'].mean()

    def _set_results(self):
        # sets the results dataframe with the selected fragments
        self.__df3 = self.__df1[['ResNum']]
        self.__df3['Probability'] = 0.0
        self.__df3['Start'] = 0
        self.__df3['End'] = 0
        self.__df3['Sequence'] = ''
        self.__df3['Charged surface'] = ''
        self.__df3['Charge (e)'] = 0.0
        self.__df3['Apolar surface'] = ''
        self.__df3['Area (nm²)'] = 0
        self.__df3['VP'] = 0
        self.__df3['VY'] = 0
        for row in self.__df3.index.tolist():
            Pmax = self.__df2[(self.__df2['PI'] <= self.__df3.at[row, 'ResNum']) & (self.__df2['PF'] >= self.__df3.at[row, 'ResNum'])]['P'].max()

            self.__df3.at[row, 'Probability'] = Pmax

            index = self.__df2.index[(self.__df2['PI'] <= self.__df3.at[row, 'ResNum']) & (self.__df2['PF'] >= self.__df3.at[row, 'ResNum']) & (self.__df2['P'] == Pmax)][0]

            self.__df3.at[row, 'Start'] = self.__df2.at[index, 'PI']

            self.__df3.at[row, 'End'] = self.__df2.at[index, 'PF']

            self.__df3.at[row, 'Sequence'] = self.__df1[(self.__df1['ResNum'] >= self.__df3.at[row, 'Start']) & (
                        self.__df1['ResNum'] <= self.__df3.at[row, 'End'])][
                'ResCode'].str.cat()

            self.__df3.at[row, 'VP'] = self.__df2.at[index, 'VP']

            self.__df3.at[row, 'VY'] = self.__df2.at[index, 'VY']

            df = self.__df1.filter(items=['ResNum', 'ResCode', 'ChargeScale', 'LinScale', 'QP'])
            df = df[(df['ResNum'] >= self.__df3.at[row, 'Start']) & (df['ResNum'] <= self.__df3.at[row, 'End'])]
            df.reset_index(drop=True, inplace=True)

            for i in df.index.tolist():
                if df.at[i, 'QP'] == self.__df3.at[row, 'VP']:
                    self.__df3.at[row, 'Charged surface'] += df.at[i, 'ResCode']
                    self.__df3.at[row, 'Charge (e)'] += df.at[i, 'ChargeScale']
                    self.__df3.at[row, 'Apolar surface'] += '_'
                elif df.at[i, 'QP'] == self.__df3.at[row, 'VY']:
                    self.__df3.at[row, 'Apolar surface'] += df.at[i, 'ResCode']
                    self.__df3.at[row, 'Area (nm²)'] += df.at[i, 'LinScale']
                    self.__df3.at[row, 'Charged surface'] += '_'
                else:
                    self.__df3.at[row, 'Charged surface'] += '_'
                    self.__df3.at[row, 'Apolar surface'] += '_'

        self.__df3.drop(['VP', 'VY'], axis=1, inplace=True)

    def __init__(self, file_path):
        # initializes the class by receiving the file path and calling the construction methods
        self.path = file_path

        # create and fill the first dataframe with peptide chain's information
        self._set_degrees()
        self._set_scales()
        self._set_delta_degrees()
        self._set_vectors()
        self._set_quadrant_points()

        # create and fill the second dataframe with random fragments and "hit" information
        self._set_random_fragments()
        self._set_fragment_vectors()
        self._set_hit()
        self._calc_prob()

        # summarize the second dataset by the residue number
        self._set_results()

    def probability(self):
        # return the score along the protein sequence
        df4 = self.__df3.filter(items=['ResNum', 'Probability'])

        return df4

    def fragments(self, cut):
        # returns the selected fragments based on the cut
        df5 = self.__df3.drop('ResNum', axis=1, inplace=False)
        df5 = df5[df5['Probability'] >= cut]
        df5.drop_duplicates(inplace=True)
        df5.sort_values(by='Probability', ascending=False, inplace=True)
        df5.reset_index(drop=True, inplace=True)

        return (df5)

class Application:
    # builds the Graphical User Interface (GUI) for the user interaction with the HAMP class

    def import_file(self):
        #command to import the PDB file for reading
        global path
        path = filedialog.askopenfilename(title="Open",initialdir="/", filetypes=[("PDB", "*.pdb")])

        self.file_path = path

        if path != '':
            self.file_path_label2.config(text=path.split("/")[-1])

        #build the HAMP class from the PDB file
        self.prot = HAMP(self.file_path)

    def analyze(self):
        # command to get the cut value and call the HAMP class methods for build the graphs and tables
        cut = self.cut_entry.get()
        self.cut = float(cut)

        # build pandastable table from .fragments method
        pt = Table(self.frame_middle, dataframe=self.prot.fragments(self.cut), editable=False, width=1290, height=110, showstatusbar=True)
        pt.cellbackgr = 'honeydew'
        pt.font = 'calibri'
        pt.fontsize = 12
        pt.redraw()
        config.load_options()
        options = {'align': 'center'}
        config.apply_options(options, pt)
        pt.autoResizeColumns()
        pt.show()

        # build the canvas for the plot and sets it
        fig = Figure(figsize=(14, 3.7),
                     dpi=100)

        y = self.prot.probability()

        fig.set_facecolor('honeydew')

        plt.clf()
        plt.style.use('seaborn-v0_8-deep')
        plot1 = fig.add_subplot(111)
        plot1.plot(y['ResNum'], y['Probability'], linewidth=1)
        plot1.set_title('CAMP Probability')
        plot1.set_xlabel("Residue Number")
        plot1.set_ylabel("Probability")
        plot1.grid(alpha=0.2)
        plot1.axhline(y=self.cut, color='grey', linestyle='--')

        canvas = FigureCanvasTkAgg(fig, master=self.frame_bottom)

        canvas.draw()

        canvas.get_tk_widget().grid(row=0, column=0, sticky='n')

    def download(self):
        #command to download the results
        Tk().withdraw()

        # Open the "Select Folder" dialog box
        folder_path = askdirectory(title='Select Folder')
        file_path1 = os.path.join(folder_path, f'predicted_CAMPs_{(self.file_path.split("/")[-1]).split('.')[0]}.csv')

        self.prot.fragments(self.cut).to_csv(file_path1, index=False)

    def __init__(self, root):
        #set main variables and the root
        self.file_path = ''
        self.cut = 0.0
        self.root = root
        self.root.title('HAMP')
        self.root.geometry("800x600")

        #build main frames and the canvas window
        self.main_frame = Frame(self.root)
        self.main_frame.pack(fill='both', expand=True)

        self.frame_top = Frame(self.main_frame, bg='darkseagreen')
        self.frame_middle = Frame(self.main_frame, bg='honeydew')
        self.frame_bottom = Frame(self.main_frame, bg='honeydew')

        self.frame_top.pack(fill='x')
        self.frame_middle.pack()
        self.frame_bottom.pack()

        #sets the widgets of the top frame
        self.title_label = Label(master = self.frame_top, text="HAMP", font=('calibri', 24, 'bold'), bg='darkseagreen', fg='ghost white')
        self.cut_label = Label(self.frame_top, text='Threshold:', font=('calibri', 14, 'bold'), bg='darkseagreen', fg='ghost white')
        self.cut_entry = Entry(self.frame_top, font=('calibri', 12), bg='ghost white', justify='center', width=8)
        self.file_path_label1 = Label(self.frame_top, text='File:', font=('calibri', 14, 'bold'), bg='darkseagreen', fg='ghost white')
        self.file_path_label2 = Label(self.frame_top, text='*.pdb', font=('calibri', 12), bg='darkseagreen', fg='ghost white', justify='center')
        self.import_button = Button(master = self.frame_top, text = 'Import', bg='honeydew2', font=('calibri', 12), command = self.import_file)
        self.analyze_button = Button(master = self.frame_top, text='Analyze', bg='honeydew2', font=('calibri', 12), command = self.analyze)
        self.download_button = Button(master=self.frame_top, text='Download', bg='honeydew2', font=('calibri', 12), command=self.download)
        self.peptide_label = Label(self.frame_top, text='Peptides:', font=('calibri', 14, 'bold'), bg='darkseagreen', fg='ghost white')

        self.spacer_label0 = Label(self.frame_top, text='', bg='darkseagreen')
        self.spacer_label1 = Label(self.frame_top, text='', bg='darkseagreen')
        self.spacer_label2 = Label(self.frame_top, text='       ', bg='darkseagreen')
        self.spacer_label3 = Label(self.frame_top, text='   ', bg='darkseagreen')
        self.spacer_label4 = Label(self.frame_top, text='   ', bg='darkseagreen')

        self.cut_entry.insert(0, "0.5")

        #organizes the widgets by grid coordenates
        self.spacer_label0.grid(row=0, column=0)
        self.title_label.grid(row = 0, column = 1)
        self.spacer_label1.grid(row = 1, column = 1)
        self.cut_label.grid(row = 2, column = 1)
        self.cut_entry.grid(row = 2, column = 2)
        self.file_path_label1.grid(row=3, column=1)
        self.file_path_label2.grid(row = 3, column = 2)
        self.spacer_label2.grid(row = 3, column = 3)
        self.import_button.grid(row = 3, column = 4)
        self.spacer_label3.grid(row=3, column=5)
        self.analyze_button.grid(row = 3, column = 6)
        self.spacer_label4.grid(row=3, column=7)
        self.download_button.grid(row = 3, column = 8)
        self.peptide_label.grid(row = 7, column = 1)

if __name__ == '__main__':
    #runs the GUI
    #root = Tk()
    #Application(root)
    #root.mainloop()
    HM = HAMP('AF-B5W6Z0-F1-model_v6.pdb')








