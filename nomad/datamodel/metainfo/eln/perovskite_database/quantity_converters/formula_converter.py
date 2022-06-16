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

import re
from pymatgen.core import Composition


cation_dict = {
    '(TMA)': 'N(CH3)3',
    '(Anyl)': 'C7H10N',
    '(CH3ND3)': 'CH3ND3',
    '(mF1PEA)': 'C16H22N2F2',
    '(pF1PEA)': 'C16H22N2F2',
    '(oFPEA)': 'C16H22N2F2',
    '(PA)': 'CH3(CH2)4NH3',
    '(NEA)': 'C12H14N',
    '(GABA)': '(HOOC(CH2)3NH3)2',
    '(mFPEA)': 'C16H22N2F2',
    'FA': 'CH5N2',
    '(4AMPY)': 'C6H10N2',
    '(PDMA)': 'C8H12N2',
    '(FPEA)': 'C8H11FN',
    '(PMA)': 'C8H12N',
    '(BZA)': 'C7H10N',
    '(5-AVA)': 'C5H12NO2',
    'DA': 'C12H28N',
    'EA': 'CH3CH2NH3',
    'GU': 'C(NH2)3',
    'GA': 'CH6N3',
    '(DMA)': 'C2H8N',
    '(pFPEA)': 'C16H22N2F2',
    '(PTA)': 'C9H14N',
    '(PEI)': 'CH2CH2',
    'Aa': 'C15H21',
    'MA': 'CH3NH3',
    '(PyrEA)': 'C7H12N2',
    'PEA': 'C6H5C2H4NH3',
    'HA': 'CH3(CH2)4CH2NH2',
    '(BEA)': 'C7H10N',
    'PA': 'C3H10N',
    'BA': 'C4H9NH3',
    '(ThMA)': 'C5H7NS',
    '(HTAB)': 'C9H22NBr',
    '(PDA)': 'NH2(CH2)3NH2',
    '(DPA)': 'C10H16N2',
    '(A43)': '(CF3)3CO(CH2)3NH3',
    '(N-EtPy)': 'C7H10N',
    '(PPA)': 'C9H11N',
    '(5-AVAI)': 'C5H12INO2',
    '(f-PEA)': 'C8H11FN',
    '(TBA)': 'C16H36N',
    'EDA': 'C2H8N2',
    '(DAP)': 'C3H10N2',
    '(CPEA)': 'ClC6H4C2H4NH3',
    '(CIEA)': 'C4H10Cl2N',
    '(TFEA)': 'C2H4F3N',
    'NEA': 'C12H11NH3',
    '(Cl-PEA)': 'C8H11ClN',
    '(HdA)': 'C6H16N2',
    'BDA': 'C8H12N2',
    '(Ace)': 'C2H7N2',
    '(oF1PEA)': 'C8H11FIN',
    '(ODA)': 'C8H21N2',
    '(BDA)': 'C8H12N2',
    '(EU-pyP)': 'C46H39EuN8O2',
    '(EDA)': 'C2H8N2',
    '(4AMP)': 'C6H14N2',
    '(3AMP)': 'C6H8N2',
    'HDA': 'C6H16N2',
    '(4ApyH)': 'C5H6N2',
    '(iso-BA)': 'C4H9NH3',
    '(AVA)': 'HOOC(CH2)4NH3I',
    '(PEA)': 'C6H5C2H4NH3',
    '(6-ACA)': 'C6H13NO2',
    '(iPA)': 'C3H8O',
    '(MIC3)': 'C5H8INS',
    '(Ada)': 'C10H17N',
    '(BdA)': 'C8H12N2',
    '(NMA)': 'C10H7CH2NH2',
    'OA': 'CH3(CH2)7NH3)',
    '(4FPEA)': 'C8H11FN',
    'NMABr': 'C11H12BrN',
    '(F-PEA)': 'C8H11FIN',
    '(ImEA)': 'C3H5N2C2H4NH2',
    '(F3EA)': 'C2H4F3N',
    '(DAT)': 'C16H10N2',
    '(TEA)': 'C6H9NS',
    'PMA': 'C8H12N',
    '(MTEA)': 'C3H9NS',
    '(3AMPY)': 'C6H8N2',
    '(BIM)': 'C6H6N4',
}

cation_dict_miss = {
    'IM': '',
    'Bn': '',
    '(CHMA)': '',
    '(THM)': '',
    '(Br-PEA)': '',
    '(EPA)': '',
    '(IEA)': '',
    '(PGA)': '',
    'BU': '',
    '(MIC2)': '',
    '(BzDA)': '',
    'IA': '',
    '(PPEA)': '',
    '(1.3-Pr(NH3)2)': '',
    'PN': '',
    'CA': '',
    '(PBA)': '',
    'DI': '',
    'PR': '',
    '(APMim)': '',
    '(FEA)': '',
    '(PyEA)': '',
    '(n-C3H7NH3)': '',
    '(BI)': '',
    'BE': '',
    '(MIC1)': '',
    '(HEA)': '',
    '(BYA)': '',
    'AN': '',
    '(H-PEA)': '',
    '(ALA)': '',
    'TA': '',
    'TN': '',
    '(HAD)': '',
    'PDA': '',
    '(F5PEA)': '',
    '(OdA)': '',
    '(ThFA)': ''}


class PerovskiteFormulaCleaner():

    def __init__(self, input_formula: str):
        """
        """
        self.input_formula = input_formula
        self.cation_dict = cation_dict
        self.cation_dict_miss = cation_dict_miss

    def replace_formula(self):
        item = self.input_formula
        if any(key in item for key in self.cation_dict_miss.keys()):
            print('''
                  The given perovskite composition contains an undefined abbreviation.
                  The composition could not be parsed.
                  ''')
        else:
            for word, replacement in sorted(self.cation_dict.items(), key=lambda x: len(x[1])):
                item = re.sub(word, replacement, item)
            output_formula = item
            return output_formula

    def clean_formula(self):
        '''
        Takes a formula and formats it into a nomad `chemical_formula_reduced`

        :return: It returns a list containing the `reduced_formula`,
        and the `elements`.
        '''
        replaced_formula = self.replace_formula()
        if replaced_formula is not None:
            try:
                composition = Composition(replaced_formula)
                int_formula = composition.get_integer_formula_and_factor()[0]
                composition_final = Composition(int_formula)
                clean_formulas_no_brackets = composition_final.get_reduced_composition_and_factor()[0]
                composition_final_int = Composition(clean_formulas_no_brackets)
                # hill_formula = composition_final_int.hill_formula
                reduced_formula = composition_final_int.get_reduced_composition_and_factor()[0].to_pretty_string()
                # reduced_formula = Composition(hill_formula).reduced_formula
                elements = composition_final_int.chemical_system.split('-')

            except ValueError:
                print('Perovskite formula with a cation abbreviation could not be parsed')

            return reduced_formula, elements


if __name__ == "__main__":

    # formulas_raw = ['Cs0.075FA0.75MA0.175PbBr0.33I2.67', 'Cs0.2FA0.8PbI3', 'Cs0.025FA0.475MA0.5Pb0.5Sn0.5I3', 'MA2CuBr3.5Cl0.5', 'FA0.95MA0.05PbBr0.15I2.85 | DA2PbI4', '(PDA)MAPb2I7', 'MAPb0.2Sn0.8Br0.4I2.6', 'CsPb1.0Br1.5I1.5', 'FA0.8MA0.2PbBr0.095I0.905', 'MASnF0.4I2.6', 'BU2FA8Pb9I28', 'MASnBr2.64I0.36', '(4AMP)MA2Pb3I10', 'Cs0.15FA0.85PbBr0.45I2.55', 'BA2MA3Pb4I12', 'Cs0.05FA0.07MA0.25PbI3', 'Cs0.04FA0.96PbBr3']
    # strValue = "FA0.95MA0.05PbBr0.15I2.85"
    # for item in formulas_raw:
    formula_cleaner = PerovskiteFormulaCleaner('(ThMA)2MA2Pb3I10')
    replaced_formula = formula_cleaner.replace_formula()
    # replaced_formula = formula_cleaner.replace_formula()
    final_formula = formula_cleaner.clean_formula()
    print(final_formula, replaced_formula)
