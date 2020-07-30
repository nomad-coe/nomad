# Copyright 2018 Fawzi Mohamed, Danio Brambila, Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os.path
import json
import numpy as np

from nomad.normalizing.normalizer import Normalizer

controlIn_basis_set = 'x_fhi_aims_section_controlIn_basis_set'
controlIn_basis_func = 'x_fhi_aims_section_controlIn_basis_func'
controlIn_nucleus = 'x_fhi_aims_controlIn_nucleus'

pure_types_json = dict()

for pure_types_str in ['light', 'really_tight', 'tight']:
    with open(os.path.join(os.path.dirname(__file__), 'data', pure_types_str + '.json')) as f:
        json_data = json.load(f)
        section_method = json_data['sections']['section_run-0']['sections']['section_method-0']
        pure_types_json[pure_types_str] = section_method[controlIn_basis_set]


class FhiAimsBaseNormalizer(Normalizer):

    # Finds out if val is in the array
    def compare_val_list(self, val, list):
        if val in list:
            return 0
        return 1

    # Comparison between two dicts.
    # List structure:
    def compare_dict_dict(self, d1, d2):
        sum2 = np.zeros(len(d2))

        # Loop over the size of dict2
        for k in np.arange(0, len(d2)):
            # Lopp over the elements of each dict2
            for idx, val in d1.items():
                # Excludes the keys that are always different.
                if (idx not in ["gIndex", "references", "uri"]):
                    try:
                        if (val != d2[k][idx]):
                            sum2[k] = sum2[k] + 1
                    except KeyError:  # this exception case arises if the cut off potential is not a number
                        continue
        if (min(sum2) == 0):
            return 0
        else:
            return 1  # sum(sum2)

    def compare_to_defaults(self, dict2_default, dict1):
        # first compare the integration grid
        false_hits_integration_grid = 0
        false_hits_basis = 0

        for key in dict1:
            if key not in ['gIndex', 'uri', controlIn_basis_func]:
                if np.size(dict1[key]) == 1:
                    if(dict1[key] != dict2_default[key]):
                        false_hits_integration_grid += 1
                        false_hits_integration_grid += abs(np.size(dict1[key]) - np.size(dict2_default[key]))
                if np.size(dict1[key]) > 1:
                    for i in dict1[key]:
                        false_hits_integration_grid += self.compare_val_list(i, dict2_default[key])
                        false_hits_integration_grid += abs(np.size(dict1[key]) - np.size(dict2_default[key]))

            elif (key == controlIn_basis_func):
                for i in np.arange(0, len(dict1[key])):
                    false_hits_basis += self.compare_dict_dict(
                        dict1[key][i], dict2_default[key])
                false_hits_basis += abs(len(dict1[key]) - len(dict2_default[key]))

        return [false_hits_integration_grid, false_hits_basis]

    def normalize(self, logger=None) -> None:
        super().normalize(logger)
        if not self.section_run or self.section_run.program_name != 'FHI-aims':
            return

        for method in self.section_run.section_method:
            to_compare = getattr(method, controlIn_basis_set, None)
            if to_compare is None:
                # not fhi aims data
                continue

            matrix_hits_int = dict.fromkeys(pure_types_json, 0)
            matrix_hits_basis = dict.fromkeys(pure_types_json, 0)

            for index, data in enumerate(to_compare):
                atom_index = int(data[controlIn_nucleus])
                for key, val in pure_types_json.items():
                    results = self.compare_to_defaults(val[atom_index], to_compare[index])
                    matrix_hits_int[key] += results[0]
                    matrix_hits_basis[key] += results[1]

                    # matrix_hits[key]=matrix_hits[key]+CompareToDefaults(val[AtomIndex],to_compare[i])

            closest_base_int = min(matrix_hits_int, key=matrix_hits_int.get)
            if (matrix_hits_basis[min(matrix_hits_basis, key=matrix_hits_basis.get)] == 0):
                closest_base_base = ''
            else:
                closest_base_base = '+'

            if (matrix_hits_int[closest_base_int] == 0):
                method.basis_set = closest_base_int + closest_base_base
            elif(matrix_hits_int[closest_base_int] <= 5):
                method.basis_set = '~' + closest_base_int + closest_base_base
            elif(matrix_hits_int[closest_base_int] > 5):
                method.basis_set = 'custom-' + closest_base_int


# import setup_paths
# import json
# from  numpy import *
# import os.path, glob
# from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
# from nomadcore.parser_backend import JsonParseEventsWriterBackend
# from nomadcore.parse_streamed_dicts import *
# import sys

# import logging

# base_path = os.path.abspath(os.path.dirname(__file__))

# files = glob.glob(os.path.join(base_path,"json_default/*.json"))

# # Open json files
# def open_json(path):
#   default=dict()
#   with open(path) as data_file:
#     default = json.load(data_file)
#   return default

# # Finds out if val is in the array
# def CompareList(val,list):
#   if val in list:
#     return 0
#   return 1

# # Comparison between two dicts.
# # List structure:
# def CompareTwoDicts(d1,d2):
#   sum2=zeros(len(d2))

# # Loop over the size of dict2
#   for k in arange(0,len(d2)):
# #   Lopp over the elements of each dict2
#     for idx,val in d1.items():
#       if (idx not in ["gIndex","references", "uri"]) : # Excludes the keys that are always different.
#         try:
#           if (val!=d2[k][idx]):
#             sum2[k]=sum2[k]+1
#         except KeyError: # this exception case arises if the cut off potential is not a number
#           continue
#   if (min(sum2)==0):
#     return 0
#   else:
#     return 1# sum(sum2)

# def CompareToDefaults(dict2_default,dict1):

# #first compare the integration grid
#   false_hits_integration_grid=0
#   false_hits_basis=0
#   for key in dict1:
#     if key not in ["gIndex", "uri", "x_fhi_aims_section_controlIn_basis_func"]:
#       if size(dict1[key])==1:
#         if(dict1[key]!=dict2_default[key]):
#           false_hits_integration_grid+=1
#         false_hits_integration_grid+=abs(size(dict1[key])-size(dict2_default[key]))
#       if size(dict1[key])>1:
#         for i in dict1[key]:
#            false_hits_integration_grid+=CompareList(i,dict2_default[key])
#         false_hits_integration_grid+=abs(size(dict1[key])-size(dict2_default[key]))
#     elif (key == "x_fhi_aims_section_controlIn_basis_func"):
#       for i in arange(0,len(dict1[key])):
#           false_hits_basis+=CompareTwoDicts(dict1[key][i],dict2_default[key])
#       false_hits_basis+=abs(len(dict1[key])-len(dict2_default[key]))
#   return [false_hits_integration_grid,false_hits_basis]

# import sys, getopt
# def main():
#   metapath = '../../../../nomad-meta-info/meta_info/nomad_meta_info/' +\
#         'common.nomadmetainfo.json'
#   metaInfoPath = os.path.normpath(
#     os.path.join(os.path.dirname(os.path.abspath(__file__)), metapath))

#   metaInfoEnv, warns = loadJsonFile(filePath=metaInfoPath,
#                                     dependencyLoader=None,
#                                     extraArgsHandling=InfoKindEl.ADD_EXTRA_ARGS,
#                                     uri=None)
#   backend = JsonParseEventsWriterBackend(metaInfoEnv)
#   calcContext = sys.argv[1]
#   backend.startedParsingSession(
#     calcContext,
#     parserInfo = {'name':'FhiAimsBasisNormalizer', 'version': '1.0'})

#   pure_types_json={}
#   matrix_hits_int={}
#   matrix_hits_basis={}


#   A=ParseStreamedDicts(sys.stdin)

#   while True:
#     to_compare=A.readNextDict()
#     if to_compare is None:
#       break
#     context=to_compare["context"]  #uncomment for prod
#     try:
#       to_compare=to_compare["section_method"]["x_fhi_aims_section_controlIn_basis_set"]
#     except:
#       continue

#     matrix_hits_int = dict.fromkeys(matrix_hits_int, 0)
#     matrix_hits_basis = dict.fromkeys(matrix_hits_basis, 0)

#     for i,d in enumerate(to_compare):
#        AtomIndex=int(d["x_fhi_aims_controlIn_nucleus"])
#        for key,val in pure_types_json.items():
#           res=CompareToDefaults(val[AtomIndex],to_compare[i])
#           matrix_hits_int[key]+=res[0]
#           matrix_hits_basis[key]+=res[1]

# #          matrix_hits[key]=matrix_hits[key]+CompareToDefaults(val[AtomIndex],to_compare[i])


#     closest_base_int=min(matrix_hits_int, key=matrix_hits_int.get)
#     if (matrix_hits_basis[min(matrix_hits_basis, key=matrix_hits_basis.get)] ==0):
#        closest_base_base=''
#     else:
#        closest_base_base='+'
#     if (matrix_hits_int[closest_base_int]==0):
# #      print(closest_base_int +closest_base_base)
#       backend.addValue("basis_set",closest_base_int +closest_base_base)
#     elif(matrix_hits_int[closest_base_int]<=5):
# #      print('~'+closest_base_int+closest_base_base)
#       backend.addValue("basis_set","~"+closest_base_int +closest_base_base)
#     elif(matrix_hits_int[closest_base_int]>5):
#       backend.addValue("basis_set",'custom-'+closest_base_int)
# #      print('custom-'+closest_base_int)

#   backend.finishedParsingSession("ParseSuccess", None)
#   sys.stdout.flush()
#   return

# if __name__ == "__main__":
#     main()
