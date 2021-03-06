import numpy as np
import re
import pandas as pd
import pylab as plt
from ase.io import read

from splitDOS import SplitDOS

""" 
rapiDOS v.0.5.1 
-- Feb. - 19 - 2018 -- 
Jose A. Garrido Torres & Michal Bajdich. 
jagt@stanford.edu  &  bajdich@slac.stanford.edu 
    
-- Mar. 2020 --   
Edited by Kirsten Winther
"""


class RapiDOS(SplitDOS):
    """ 

    """

    def __init__(self):
        """Get general information of the system"""
        self.atoms = read('CONTCAR')  # Open CONTCAR file
        self.ispin = self.get_spin()
        self.kblock = self.get_kblock()

        # Call SplitDOS
        super().__init__()
        self.write_dos0()
        if self.ispin == 2:
            self.write_spin()
        else:
            self.write_nospin()

    def get_spin(self):

        incar_file = open("INCAR", "r")
        ispin = 1  # Non spin polarised calculations.
        for line in incar_file:
            if re.match("(.*)ISPIN(.*)2", line):
                ispin = 2  # For spin polarised calculations.
        return ispin

    def get_kblock(self):
        """Check scale"""
        # KBLOCK scales DOSCAR (See VASP Manual).
        for line in open('OUTCAR'):
            if line.find('KBLOCK') != -1:
                ckblock = line.split()
                index = [i+2 for i,
                         ck in enumerate(ckblock) if ck == 'KBLOCK'][0]
                kblock = float(ckblock[index])
                return kblock

    def get_total_dos(self):
        """ Get np array and columns for total DOS"""

        # Load file total DOS from DOS0.
        total_dos = np.loadtxt('DOS0', skiprows=0)
        # Scale (See kblock VASP).
        total_dos[:, 1:] = total_dos[:, 1:] * self.kblock
        self.total_dos = total_dos

        if self.ispin == 2:
            self.dos_columns = ['Energy (E-Ef)',
                                'Total DOS Spin Up',
                                'Total DOS Spin Down',
                                'Integrated DOS Spin Up',
                                'Integrated DOS Spin Down']
        elif self.ispin != 2:
            self.dos_columns = ['Energy (E-Ef)',
                                'Total DOS Spin Up',
                                'Integrated DOS']

        # total_dos_df = pd.DataFrame(total_dos,
        #                            columns=total_dos_columns)
        # total_dos_df.to_csv('TotalDOS.csv')

        return (self.dos_columns, self.total_dos)

    def get_pdos(self):
        """Build arrays for PDOS"""
        if self.ispin != 2:
            self.pdos_columns = ['Energy (E-Ef)',
                                 's_up',
                                 'py_up', 'pz_up', 'px_up',
                                 'dxy_up', 'dyz_up', 'dz2_up', 'dxz_up', 'dx2_up',
                                 'f1_up', 'f2_up', 'f3_up',
                                 'f4_up', 'f5_up', 'f6_up', 'f7_up'][:self.ncols]

        elif self.ispin == 2:
            self.pdos_columns = ['Energy (E-Ef)',
                                 's_up', 's_down',
                                 'py_up', 'py_down',
                                 'pz_up', 'pz_down',
                                 'px_up', 'px_down',
                                 'dxy_up', 'dxy_down',
                                 'dyz_up', 'dyz_down',
                                 'dz2_up', 'dz2_down',
                                 'dxz_up', 'dxz_down',
                                 'dx2_up', 'dx2_down',
                                 'f1_up', 'f1_down',
                                 'f2_up', 'f2_down',
                                 'f3_up', 'f3_down',
                                 'f4_up', 'f4_down',
                                 'f5_up', 'f5_down',
                                 'f6_up', 'f6_down',
                                 'f7_up', 'f7_down'][:self.ncols]

        atomic_symbols = self.atoms.get_chemical_symbols()
        pdos_data = {}
        for i in range(len(self.atoms)):
            pdos_data[i] = np.loadtxt('DOS'+str(i+1), skiprows=0)
        return self.pdos_columns, pdos_data

    def get_bandgap(self):
        # Load file total DOS from DOS0.
        total_dos = np.loadtxt('DOS0', skiprows=0)
        arg_fermi = np.where(total_dos[:, 0] > 0)[0][0]
        arg_fermi_upper = arg_fermi
        arg_fermi_lower = arg_fermi
        band_tolerance = np.max(
            total_dos[:, 1] + total_dos[:, 2]) / 100  # 5e-3

        converged = False
        while not converged:
            print(arg_fermi_lower, arg_fermi_upper)
            occup_lower = (total_dos[:, 1] + total_dos[:, 2])[arg_fermi_lower]
            occup_upper = (total_dos[:, 1] + total_dos[:, 2])[arg_fermi_upper]

            if occup_lower > band_tolerance and occup_upper > band_tolerance:
                e_lower = total_dos[arg_fermi_lower][0]
                e_upper = total_dos[arg_fermi_upper][0]
                band_gap = e_lower - e_upper
                converged = True

            if occup_lower < band_tolerance:
                arg_fermi_lower -= 1
            if occup_upper < band_tolerance:
                arg_fermi_upper += 1

        print("Approx. band gap: ", np.abs(band_gap), "eV")

        band_gap_columns = ['Lower Band Gap', 'Upper Band Gap', 'Band Gap']
        band_gap_df = pd.DataFrame([[e_lower, e_upper, band_gap]],
                                   columns=band_gap_columns)

        return band_gap_df

    def get_bandcenter(self, energy, dos):
        center = np.sum(energy*dos) / np.sum(dos)

        return center

    def get_fillfactor(self, energy, dos):

        full = np.sum(dos)

        filled = np.sum(dos[energy <= 0])

        return filled/full

    def plot_dos(self, xlim=[-10, 5], title="", ext='eps'):
        dos_columns, total_dos = self.get_total_dos()
        band_gap = self.get_bandgap()
        band_gap_lower = band_gap['Lower Band Gap'][0]
        band_gap_upper = band_gap['Upper Band Gap'][0]
        print('Approx. Band Gap:', np.round(
            np.abs(band_gap['Band Gap'][0]), 3), "eV")

        ix = dos_columns.index('Energy (E-Ef)')
        iy = dos_columns.index('Total DOS Spin Up')

        energy = total_dos[:, ix]
        dos = total_dos[:, iy]

        fig = plt.figure(figsize=(10.0, 6.0))  # Create figure.

        plt.plot(energy,
                 dos,
                 color='r',
                 label='spin up')

        plt.fill_between(energy, 0, dos,
                         facecolor='lightpink',
                         interpolate=True)  # Fill between spin up and down.

        dos_sum = dos.copy()
        if self.ispin == 2:
            iy = dos_columns.index('Total DOS Spin Down')
            dos_down = -total_dos[:, iy]
            dos_sum -= dos_down
            plt.plot(energy,
                     dos_down,
                     color='b',
                     label='spin down')
            # Fill between spin up and down.
            plt.fill_between(energy,
                             0,
                             dos_down,
                             facecolor='lightblue', interpolate=True)

        plt.plot(energy,
                 dos_sum,
                 color='k',
                 label='total')

        # Plot vertical line in Fermi
        plt.axvline(x=[0.0], color='k', linestyle='--', linewidth=1.2)
        plt.axvspan(band_gap_lower, band_gap_upper, color='silver',
                    label='band gap={}'.format(band_gap['Band Gap'][0].round(2)))
        plt.axhline(y=[0.0], color='k', linestyle='--', linewidth=1.2)

        # plt.axvline(x=[center], color='g', linestyle='--', linewidth=1.2,
        #            label='band center={}'.format(np.round(center, 2)))
        plt.legend()
        plt.xlabel('(E - E$_F$) [eV]')  # x axis label.
        plt.ylabel('DOS [states eV$^{-1}$]')  # x axis label.
        plt.xlim(xlim)  # Plot limits.
        plt.title(title)
        plt.show()
        fig.savefig('DOS{}.{}'.format(title, ext))

    def plot_pdos(self, elements=None,
                  title='PDOS',
                  ext='eps',
                  xlim=[-10, 5]):
        """
        Plot projected density of states
        """

        dos_columns, total_dos = self.get_total_dos()

        ix = dos_columns.index('Energy (E-Ef)')
        iy = dos_columns.index('Total DOS Spin Up')

        energy = total_dos[:, ix]
        idx1 = np.min(np.argwhere(energy >= xlim[0]))
        idx2 = np.max(np.argwhere(energy <= xlim[1]))
        total_dos = total_dos[idx1:idx2]
        energy = total_dos[:, ix]
        dos = total_dos[:, iy]

        if self.ispin == 2:
            iyd = dos_columns.index('Total DOS Spin Down')
            dos_down = -total_dos[:, iyd]

        pdos_columns, pdos_data = self.get_pdos()
        for k in pdos_data.keys():
            pdos_data[k] = pdos_data[k][idx1:idx2]
        symbols = self.atoms.get_chemical_symbols()
        all_elements = list(set(symbols))
        if not elements:
            elements = {}
            for e in all_elements:
                elements.update({e: "all"})

        color_dict = {'Oall': "red", 'H': "grey", 'N': "blue", 'S': 'yellow'}
        colors = ['seagreen', 'orange', 'skyblue', 'pink',
                  'magenta', 'indigo', 'coral']
        j = 0
        for e, orb in elements.items():
            for o in orb:
                if e + o not in color_dict:
                    color_dict.update({e+o: colors[j]})
                    j += 1
        print(color_dict)
        # Plots
        fig = plt.figure(figsize=(10.0, 6.0))  # Create figure.

        plt.plot(energy,
                 dos,
                 color='black',
                 label='Total DOS')  # Plot DOS spin up.

        if self.ispin == 2:
            plt.plot(energy,
                     dos_down,
                     label='',
                     color='black')  # Plot DOS spin down.

        for e, orbitals in elements.items():
            if e[-1].isdigit():
                idx = [i for i in range(len(e)) if e[i].isdigit()][0]

                atom_indices = [int(e[idx:])]
                print(atom_indices)
            else:
                atom_indices = [i for i, s in enumerate(symbols) if s == e]

            if not isinstance(orbitals, list):
                orbitals = [orbitals]

            for orbital_plot in orbitals:
                if orbital_plot == 'all':
                    orb_indices_up = [i for i, p in enumerate(pdos_columns)
                                      if p[-2:] == 'up']
                    orb_indices_down = [i for i, p in enumerate(pdos_columns)
                                        if p[-4:] == 'down']
                else:
                    if orbital_plot in ['s', 'p', 'd', 'f', 't2g', 'eg']:
                        orbitals = self.expand_orbitals(orbital_plot)
                    else:
                        orbitals = [orbital_plot]

                    orbitals_up = [o + '_up' for o in orbitals]
                    orbitals_down = [o + '_down' for o in orbitals]
                    orb_indices_up = [io for io, o in enumerate(pdos_columns) if
                                      o in orbitals_up]
                    orb_indices_down = [io for io, o in enumerate(pdos_columns) if
                                        o in orbitals_down]

                pdos_up = 0
                pdos_down = 0
                for i in atom_indices:
                    pdos_up += np.sum(pdos_data[i][:, orb_indices_up], axis=1)
                    pdos_down += np.sum(pdos_data[i]
                                        [:, orb_indices_down], axis=1)

                color = color_dict[e + orbital_plot]
                extra_label = ''

                s_labels = ['up', 'down']
                labels = [None, None]
                sgn = 1
                for i, pdos_i in enumerate([pdos_up, pdos_down]):
                    s_label = s_labels[i]
                    if len(pdos_i) == 0:
                        continue
                    if i == 0:
                        label = '{}({})'.format(e, orbital_plot)
                    else:
                        label = None
                    plt.plot(energy,
                             pdos_i,
                             color=color,
                             label=label)

        plt.axvline(x=[0.0], color='k', linestyle='--', linewidth=1.2)

        plt.xlabel('(E - E$_F$) [eV]')  # x axis label.
        plt.ylabel('(P)DOS [states eV$^{-1}$]')  # x axis label.
        plt.legend()
        plt.title(title)
        plt.show()

        fig.savefig('PDOS.eps')  # {}.{}'.format(title, ext))

    def expand_orbitals(self, name):
        if name in ['s', 'p', 'd', 'f']:
            all_orbitals = list(set([orb.split('_')[0] for orb in
                                     self.pdos_columns if orb[0] == name]))
        elif name == 't2g':
            all_orbitals = ['dxy', 'dyz', 'dxz']
        elif name == 'eg':
            all_orbitals = ['dz2', 'dx2']

        return all_orbitals
