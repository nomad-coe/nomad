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

import numpy as np

from nomad.datamodel.metainfo.simulation.system import Descriptors, SOAP

from .normalizer import SystemBasedNormalizer


class SoapNormalizer(SystemBasedNormalizer):
    def normalize_system(self, system, is_representative):
        try:
            from quippy import descriptors
        except ImportError:
            descriptors = None

        # Only store SOAP for representative system to start with
        if not is_representative:
            return True

        # structures used to be stored in results, adding them back like this for now...
        if system.atoms is None:
            return False

        if not descriptors:
            self.logger.warning('SOAP normalizer runs, but quippy is not installed.')
            return False

        # TODO compute descriptors for primitive system, need to discuss how to store this stuff.
        # repr_symmetry = system.symmetry[0]
        # print(repr_symmetry)
        # symmetry_analyzer = repr_symmetry.m_cache.get("symmetry_analyzer")
        # prim_atoms = symmetry_analyzer.get_primitive_system()
        # print("prim_atoms is", type(prim_atoms), prim_atoms)

        soap = SOAP()
        atoms = system.atoms.to_ase(raise_exp=True)

        # setup params to be used by quippy
        params, _ = dataset_to_params(atoms)
        N, S, L = params["n_max"], params["n_Z"], params["l_max"]
        soap.n_max = N
        soap.l_max = L
        soap.r_cut = np.float64(params["soap cutoff"])
        soap.atom_sigma = np.float64(params["atom_sigma"])

        # regular soap
        quippy_str = params_to_quippy_str(params)
        desc = descriptors.Descriptor(quippy_str)
        output = desc.calc(atoms)
        ps_array_list = [flat_to_array(p, N, L, S) for p in output["data"]]
        soap.soap = np.array(
            [np.reshape(d, (S, S, -1)) for d in ps_array_list]
        )

        # global soap
        quippy_str += " average=T"
        desc = descriptors.Descriptor(quippy_str)
        output = desc.calc(atoms)
        A = flat_to_array(output["data"][0], N, L, S)
        soap.global_soap = np.reshape(A, (S, S, -1))

        # tensor reduced soap
        quippy_str = params_to_quippy_str(params)
        quippy_str += "sym_mix=F Z_mix=T R_mix=T K=50 coupling=F"
        desc = descriptors.Descriptor(quippy_str)
        output = desc.calc(atoms)
        soap.tr_soap = output["data"]

        # gloabl tensor reduced soap
        quippy_str += " average=T"
        desc = descriptors.Descriptor(quippy_str)
        output = desc.calc(atoms)
        soap.global_tr_soap = output["data"][0]

        # add the soap descriptor to the system
        if not system.descriptors:
            system.descriptors = Descriptors()
        system.descriptors.soap = soap

        return True


def params_to_quippy_str(params):
    qs = ""
    for key, val in params.items():
        qs += str(key) + "=" + str(val) + " "
    return qs


def dataset_to_params(atoms):
    """convert the dataset into params"""
    Zs = set(atoms.numbers)
    Zs = sorted(list(Zs), key=lambda x: x)
    Zstr = "{"
    for Z in Zs:
        Zstr += str(Z) + " "
    Zstr = Zstr[:-1] + "}"

    params = {}
    params["n_Z"] = len(Zs)
    params["n_species"] = len(Zs)
    params["Z"] = Zstr
    params["species_Z"] = Zstr
    params["n_max"] = 8
    params["l_max"] = 4
    params["soap cutoff"] = 5
    params["atom_sigma"] = 0.4
    params["cutoff_transition_width"] = 0.5
    params["central_weight"] = 1
    return params, Zs


def flat_to_array(ps_full, N, L, S, check=False):
    # generate rs_index
    rs_index = []
    for alpha in range(0, S):
        for a in range(0, N):
            rs_index.append([a, alpha])

    # initiialise power spec
    d = np.zeros(shape=(S, S, L + 1, N, N))

    # assign entries
    ipow = 0
    for ia in range(0, N * S):
        a, i_species = rs_index[ia]
        for jb in range(0, ia + 1):
            b, j_species = rs_index[jb]
            for ll in range(0, L + 1):
                val = ps_full[ipow]
                if ia != jb:
                    val /= 2**0.5
                    d[i_species, j_species, ll, a, b] = val
                    d[j_species, i_species, ll, b, a] = val
                else:
                    d[i_species, j_species, ll, a, b] = val
                ipow += 1
    if check:
        d_sort = sorted(d.flatten())
        fp_sort = sorted(ps_full)
        assert np.allclose(d_sort, fp_sort)
    return d
