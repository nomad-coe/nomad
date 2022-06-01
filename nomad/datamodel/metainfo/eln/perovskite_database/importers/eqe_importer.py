# -*- coding: utf-8 -*-
#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

###################################################################################################
# Evaluation of EQE measurement data + Urbach tail to determine the radiative open-circuit voltage.
# by: Lisa Krückemeier & Dane W. deQuilettes
# (supplemented version of 'SCRIPT_Vocrad_EQE_fit_Urbachtail.m' by Lisa Krückemeier)

# translated to Python 2.x and 3.x by Christian Wolff
###################################################################################################
##
###################################################################################################
from scipy import integrate
import numpy as np
import pandas as pd
import os
from scipy.signal import savgol_filter


def read_eqe(filename, header_lines=None):

    # header_lines = None
    if header_lines is None:
        header_lines = 0

    temperature = 300  # in [°K]

    # some optionals:
    ###########################################################################
    # !!!! VARIABLE: ADJUST URBACH ENERGY IF NECESSARY
    E_Urbach = 0.013                        # [eV] Urbach energy (usually between 13.5 to 16 meV for metal-halide perovskites)
    ###########################################################################
    # Decide if you want to choose the EQE value for Urbach Tail attachment manually or automatically
    manually = 'yes'                          # choose 'yes' or 'no' if you want to chose the transition point manually 'yes' if a manual correction is necessary or automatically 'no'. 'no' attaches Urbach tail midway between data on logscale.
    EQE_level = 0.01                         # EQE value at which the Urbach Tail should be attached (only if manually 'yes').

    coords = []

    if header_lines == 0:  # in case you have a header
        try:
            df = pd.read_csv(filename, header=None, sep='\t',)
        except IndexError:
            df = pd.read_csv(filename, header=None)

    else:
        try:
            df = pd.read_csv(filename, header=int(header_lines - 1), sep='\t')

        except IndexError:
            df = pd.read_csv(filename, header=int(header_lines - 1))

    # df = df.apply(pd.to_numeric, errors='coerce')
    # df = df.dropna()

    if 'Calculated' in list(df.columns):
        x = df[df.columns[0]]
        y = df['Calculated'].values
    else:
        x = df[df.columns[0]]
        y = df[df.columns[1]]

    x = np.array(x)
    y = np.array(y)

    if any(x > 10):  # check if energy (eV) or wavelength (nm)
        x = 1240 / x

    if any(y > 10):  # check if EQE is given in (%), if so it's translated to abs. numbers
        y = y / 100

    if x[1] - x[2] > 0:  # bring both arrays into correct order (i.e. w.r.t eV increasing) if one started with e.g. wavelength in increasing order e.g. 300nm, 305nm,...
        x.sort()
        y = np.flip(y)

    xold = x
    x = np.linspace(min(x), max(x), 500, endpoint=True)
    y = np.interp(x, xold, y)
    # print(y)
    coords = []

    # plt.waitforbuttonpress()
    if manually == 'no':
        pass

    # here you push the last part of the script
    q = 1.602176462e-19  # % [As], elementary charge
    h_Js = 6.62606876e-34  # % [Js], Planck's constant
    # h_eVs = 4.135667662e-15  # % [eVs], Planck's constant
    k = 1.38064852e-23  # % [(m^2)kg(s^-2)(K^-1)], Boltzmann constant
    T = temperature
    VT = (k * T) / q  # % [V], 25.8mV thermal voltage at 300K
    c = 299792458  # % [m/s], speed of light c_0
    # vor = ((h_Js * c) / q) / (1e-9)  # % prefactor for converting between energy and wavelength in (eVnm)

    x_interp = np.linspace(min(x), max(x), 1000, endpoint=True)
    y_interp = np.interp(x_interp, x, y)
    y_interp = savgol_filter(y_interp, 51, 4)

    if manually == 'yes':

        # TODO find an automatic way to auto find the optimun fitting range
        y_fit_min = min(y) * 5
        y_fit_max = y_fit_min * 8

        x_value = np.interp(y_fit_min, y, x)
        x_value_2 = np.interp(y_fit_max, y, x)
        # print(y_fit_min, y_fit_max, x_value, x_value_2)
        i_x1 = np.abs(x - x_value).argmin()
        i_x2 = np.abs(x - x_value_2).argmin()
        x1 = x[min(i_x1, i_x2)]
        x2 = x[max(i_x1, i_x2)]

        # y1 = y[min(i_x1, i_x2)]
        y2 = y[max(i_x1, i_x2)]

        coords = coords[-2:]
        yfit = y[min(i_x1, i_x2):max(i_x1, i_x2)]
        xfit = x[min(i_x1, i_x2):max(i_x1, i_x2)] - x1
        yfit2 = np.log(yfit)
        log_slope = np.mean(np.diff(yfit2) / np.diff(xfit))
        E_Urbach = 1 / log_slope
        x_extrap = np.linspace(-1, 0, 500, endpoint=False) + x2
        y_extrap = y2 * np.exp((x_extrap - x2) * log_slope)

        x_interp = np.linspace(x2, max(x), 1000, endpoint=True)
        y_interp = np.interp(x_interp, x[max(i_x1, i_x2):], y[max(i_x1, i_x2):])

        x_interp = x_interp[y_interp >= EQE_level]
        y_interp = y_interp[y_interp >= EQE_level]
        x_extrap = np.linspace(-1, 0, 500, endpoint=False) + min(x_interp)
        y_extrap = y_interp[0] * np.exp((x_extrap - min(x_interp)) / E_Urbach)
    else:
        i_x1 = np.abs(x - coords[0][0]).argmin()
        i_x2 = np.abs(x - coords[1][0]).argmin()
        x1 = x[min(i_x1, i_x2)]
        x2 = x[max(i_x1, i_x2)]
        # y1 = y[min(i_x1, i_x2)]
        y2 = y[max(i_x1, i_x2)]

        coords = coords[-2:]
        yfit = y[min(i_x1, i_x2):max(i_x1, i_x2)]
        xfit = x[min(i_x1, i_x2):max(i_x1, i_x2)] - x1
        yfit2 = np.log(yfit)
        log_slope = np.mean(np.diff(yfit2) / np.diff(xfit))
        E_Urbach = 1 / log_slope
        x_extrap = np.linspace(-1, 0, 500, endpoint=False) + x2
        y_extrap = y2 * np.exp((x_extrap - x2) * log_slope)

        x_interp = np.linspace(x2, max(x), 1000, endpoint=True)
        y_interp = np.interp(x_interp, x[max(i_x1, i_x2):], y[max(i_x1, i_x2):])

    x_new = np.append(x_extrap, x_interp)
    y_new = np.append(y_extrap, y_interp)

    phi_BB = (2 * 3.14159265 * q**3 * (x_new)**2) / (h_Js**3 * c**2 * (np.exp(x_new / VT) - 1))
    EL = phi_BB * y_new
    # j0rad = []
    # for i in range(len(EL)):

    j0rad = integrate.cumtrapz(EL, x_new)
    j0rad = j0rad * q

    dir_path = os.path.dirname(os.path.realpath(__file__))
    filename2 = os.path.join(dir_path, 'AM15G.dat.txt')

    df3 = pd.read_csv(filename2, header=None)

    # print(df3)
    energy_AM15 = df3[df3.columns[1]]
    energy_AM15 = np.array(energy_AM15)
    spectrum_AM15 = df3[df3.columns[2]]
    spectrum_AM15 = np.array(spectrum_AM15)

    spectrum_AM15G_interp = np.interp(x_new, energy_AM15, spectrum_AM15)
    JSC_calc = integrate.cumtrapz(y_new * spectrum_AM15G_interp, x_new)
    JSC_int = max(JSC_calc * q * 1e4)

    VOC_rad_calc = VT * np.log(JSC_int / max(j0rad))

    dEQE_interp = np.diff(y_new) / np.diff(np.flip(-x_new))
    e_g = x_new[dEQE_interp.argmax()]
    # print('band gap is max(d/dE (EQE) = %s eV' %(E_G)

    eqe_dict = {}
    eqe_dict['bandgap_eqe'] = round(e_g, 3)
    eqe_dict['integrated_Jsc'] = round(JSC_int, 3)
    eqe_dict['eqe_array'] = y_new
    eqe_dict['photon_energy_array'] = x_new

    if E_Urbach <= 0.026:
        eqe_dict['integrated_J_0_rad'] = float('{:0.3e}'.format(max(j0rad)))
        eqe_dict['voc_rad'] = round(VOC_rad_calc, 3)
        eqe_dict['urbach_e'] = round(E_Urbach, 4)

    else:
        pass

    # print(eqe_dict['urbach_e'])
    # print(eqe_dict['bandgap_eqe'])

    return eqe_dict


filename = '/home/pepe_marquez/NOMAD/nomad/nomad/datamodel/metainfo/eln/perovskite_database/importers/15 of iii batch.txt'
# filename = '/home/pepe_marquez/NOMAD/nomad/nomad/datamodel/metainfo/eln/perovskite_database/importers/EQE_Liu_ACSEnergyLett_19_recipeB.dat'
# read_eqe(filename, header_lines=0)
read_eqe(filename, 9)
