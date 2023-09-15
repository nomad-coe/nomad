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

import os
import re
from typing import Dict

import numpy as np
from nomad.units import ureg
try:
    from pyzeo.netstorage import AtomNetwork
    from pyzeo.area_volume import volume, surface_area
except ImportError:
    # NOTE: pyzeo is an optional dependency
    pass
from pymatgen.io.ase import AseAtomsAdaptor

from nomad.normalizing.topology import add_system, add_system_info, get_topology_original
from nomad.normalizing.normalizer import SystemBasedNormalizer
from nomad.normalizing.common import ase_atoms_from_nomad_atoms
from nomad.normalizing import mof_deconstructor
from nomad.datamodel.results import System, Material, Relation

species_re = re.compile(r'^([A-Z][a-z]?)(\d*)$')


class PorosityNormalizer(SystemBasedNormalizer):
    '''
    This normalizer creates general parameters that are needed
    to fully describe the geometric proporties of a metal-organic framework (MOFs).
    The normalizer takes MOFs, deconstruct it into its building units.
    These building units correspond to the metal cluster, organic ligand,
    metal secondary building units (sbu) and organic sbu.
    1) The sbus represents the building units as nodes and edges. This information
    can be relevant for automatic building of hypthetical MOFs
    2) The organic linkers are relevant for understanding the chemistry of the MOF

    Moreover, the normalizer also performs a geometrical analysis of the MOF
    and provides information regarding the porosity. These are the accessible surface area,
    pore limiting diamter (pld), local cavity diamter (lcd), accessible volume (AV),
    number of channels and void fraction.

    In the future, the it will also compute the rcsr tological code.
    '''

    def __init__(self, archive):
        super().__init__(archive, only_representatives=True)
        return

    def normalize_system(self, system, is_representative):
        if not is_representative:
            return False

        if system.atoms is None:
            return

        if not system.atoms or system.atoms.positions is None or system.atoms.lattice_vectors is None or system.atoms.labels is None:
            return

        if len(system.atoms.lattice_vectors) != 3:
            return

        if len(system.atoms.positions) != len(system.atoms.labels):
            return

        topology: Dict[str, System] = {}
        original_atoms = ase_atoms_from_nomad_atoms(system.atoms)

        # TODO: We enforce full periodicity for the system because the GULP
        # parser does not currenlty parse it correctly.
        original_atoms.set_pbc(True)

        if len(original_atoms) < 20:
            return
        number_of_metals, _ = mof_deconstructor.metal_coordination_number(
            original_atoms)
        if len(number_of_metals) == 0 and len(original_atoms) > 500:
            return

        if len(original_atoms) > 5000:
            return

        indices_in_parent = mof_deconstructor.remove_unbound_guest(
            original_atoms)

        if len(indices_in_parent) > 0:
            atoms = original_atoms[indices_in_parent]
        else:
            atoms = original_atoms
        guest_indices = [
            i.index for i in original_atoms if i.index not in indices_in_parent]

        if not mof_deconstructor.inter_atomic_distance_check(atoms):
            return

        number_of_carbon = [i.symbol for i in atoms if i.symbol in ['C']]
        aluminium_silicon = [
            i.symbol for i in atoms if i.symbol in ['Al', 'Si']]
        if len(number_of_carbon) == 0 and len(aluminium_silicon) == 0:
            return

        try:
            porosity_data = zeo_calculation(atoms)
        except Exception:
            return

        if len(porosity_data) != 7:
            return

        if porosity_data['Number_of_channels'] == 0:
            return

        if porosity_data['PLD'] < 1.86:
            return

        parent_system = get_topology_original(system.atoms, self.entry_archive)
        no_guest = porosity_system(
            original_atoms, porosity_data, indices_in_parent)
        add_system(parent_system, topology)
        add_system_info(parent_system, topology)
        add_system(no_guest, topology, parent_system)
        add_system_info(no_guest, topology)

        if no_guest.label == 'MOF':
            try:
                mof_properties(original_atoms, indices_in_parent,
                               parent_system, topology)
            except Exception:
                pass

        if len(guest_indices) > 0:
            guest = System(
                indices=[guest_indices],
                label="guest",
                structural_type='group',
                system_relation=Relation(type='group'),
                method='porosity',
            )
            add_system(guest, topology, parent_system)
            add_system_info(guest, topology)
        # Add topology to the archive
        material = self.entry_archive.m_setdefault('results.material')
        for system in topology.values():
            material.m_add_sub_section(Material.topology, system)


def porosity_system(original_atoms, porosity_data, indices_in_parent):
    atoms = original_atoms[indices_in_parent]
    porous_system_label = 'Porous system'
    number_of_metals, metal_coordination = mof_deconstructor.metal_coordination_number(
        atoms)
    number_of_carbon = [i.symbol for i in atoms if i.symbol in ['C']]
    organic_entity_on = [i.symbol for i in atoms if i.symbol in ['N', 'O']]
    aluminium_silicon = [i.symbol for i in atoms if i.symbol in ['Al', 'Si']]
    oxygen = [i.symbol for i in atoms if i.symbol in ['O']]

    if len(number_of_metals) > 0 and len(number_of_carbon) > 0 and len(organic_entity_on) > 0:
        porous_system_label = 'MOF'
    elif len(aluminium_silicon) > 0 and len(oxygen) > 0:
        porous_system_label = 'Zeolite'
    elif len(number_of_metals) == 0 and len(number_of_carbon) > 0 and len(organic_entity_on) > 0:
        porous_system_label = 'COF'
    no_guest = System(
        structural_type='bulk',
        system_relation=Relation(type='group'),
        label=porous_system_label,
        indices=[indices_in_parent],
        # Porous Unique properties
        largest_cavity_diameter=porosity_data['LCD'] * ureg.angstrom,
        pore_limiting_diameter=porosity_data['PLD'] * ureg.angstrom,
        largest_included_sphere_along_free_sphere_path=porosity_data['lfpd'] * ureg.angstrom,
        accessible_surface_area=porosity_data['ASA'] * ureg.angstrom ** 2,
        accessible_volume=porosity_data['AV'] * ureg.angstrom ** 3,
        void_fraction=porosity_data['AV_Volume_fraction'],
        n_channels=porosity_data['Number_of_channels'],
        method='porosity',
    )
    if len(number_of_metals) > 0:
        no_guest.metal_coordination = metal_coordination
    return no_guest


def mof_properties(original_atoms, indices_in_parent, parent_system, topology):
    atoms = original_atoms[indices_in_parent]
    number_of_metals, _ = mof_deconstructor.metal_coordination_number(
        atoms)
    if len(number_of_metals) == 0:
        return

    organic_entity_carbon = [i.symbol for i in atoms if i.symbol in ['C']]
    if len(organic_entity_carbon) == 0:
        return

    organic_entity_on = [i.symbol for i in atoms if i.symbol in ['N', 'O']]
    if len(organic_entity_on) == 0:
        return
    indices_new_atom = [i.index for i in original_atoms]
    map_indices = dict(zip(indices_new_atom, indices_in_parent))

    list_of_connected_components, atom_pairs_at_breaking_point, porphyrin_checker, all_regions = mof_deconstructor.\
        secondary_building_units(atoms)
    metal_sbus, organic_sbus, _ = mof_deconstructor.\
        find_unique_building_units(
            list_of_connected_components,
            atom_pairs_at_breaking_point,
            atoms,
            porphyrin_checker,
            all_regions
        )

    list_of_connected_components, atom_pairs_at_breaking_point, porphyrin_checker, all_regions = mof_deconstructor.\
        ligands_and_metal_clusters(atoms)
    _, organic_linker, _ = mof_deconstructor.\
        find_unique_building_units(
            list_of_connected_components,
            atom_pairs_at_breaking_point,
            atoms,
            porphyrin_checker,
            all_regions
        )
    # SBU Groups
    group_counter = 1
    for i, sbu in enumerate(metal_sbus):
        sbu_group = System(
            structural_type='group',
            system_relation=Relation(type='group'),
            method='porosity',
        )
        sbu_group.label = f'metal_sbu_{i+1}'
        sbu_indices = [map_indices[i]
                       for i in sum(sbu.info['atom_indices_mapping'], [])]
        sbu_group.indices = [sbu_indices]
        group_counter += 1
        add_system(sbu_group, topology, parent_system)
        add_system_info(sbu_group, topology)

        mof_group_indices = [map_indices[i]
                             for i in sbu.info['atom_indices_mapping'][0]]
        mof_group = System(
            indices=[mof_group_indices],
            n_atoms=len(mof_group_indices),
            label=f'metal_sbu_{i+1}',
            structural_type='molecule',
            method='porosity',
        )
        mof_group.sbu_type = sbu.info['sbu_type']
        mof_group.sbu_coordination_number = len(sbu.info['point_of_extension'])
        add_system(mof_group, topology, sbu_group)
        add_system_info(mof_group, topology)
    # -------------------------------------------------------------
    for i, linker in enumerate(organic_sbus):
        sbu_group = System(
            structural_type='group',
            system_relation=Relation(type='group'),
            method='porosity',
        )
        sbu_indices = [map_indices[i]
                       for i in sum(linker.info['atom_indices_mapping'], [])]
        sbu_group.label = f'organic_sbu_{i+1}'
        sbu_group.indices = [sbu_indices]
        sbu_group.index = group_counter
        sbu_group.type = 'linker'
        add_system(sbu_group, topology, parent_system)
        add_system_info(sbu_group, topology)
        # Single SBU
        mof_group_indices = [map_indices[i]
                             for i in linker.info['atom_indices_mapping'][0]]
        mof_group = System(
            indices=[mof_group_indices],
            label=f'organic_sbu_{i+1}',
            structural_type='molecule',
            method='porosity',
        )
        mof_group.sbu_coordination_number = len(
            linker.info['point_of_extension'])
        add_system(mof_group, topology, sbu_group)
        add_system_info(mof_group, topology)
    # ------------------------------------------------------------
    for i, ligand in enumerate(organic_linker):
        sbu_group = System(
            structural_type='group',
            system_relation=Relation(type='group'),
            method='porosity',
        )
        mof_group_indices = [map_indices[i]
                             for i in ligand.info['atom_indices_mapping'][0]]
        sbu_indices = [map_indices[i]
                       for i in sum(ligand.info['atom_indices_mapping'], [])]
        sbu_group.label = f'ligand_{i+1}'
        sbu_group.indices = [sbu_indices]
        sbu_group.type = 'ligand'
        add_system(sbu_group, topology, parent_system)
        add_system_info(sbu_group, topology)
        mof_group = System(
            indices=[mof_group_indices],
            label=f'organic_ligand_{i+1}',
            structural_type='molecule',
            method='porosity',
        )
        sbu_group.label = f'ligand_{i+1}'
        add_system(mof_group, topology, sbu_group)
        add_system_info(mof_group, topology)


def put_contents(filename, output):
    '''
    Create a new writable object
    Args:
        filename : path of file to write
    Return:
        output : list of strings to write
    '''
    with open(filename, 'w', encoding='utf-8') as f:
        f.writelines(output)
    return


def get_contents(filename):
    '''
    A python function to read a text and return a list of lines
    Args:
        filename
    Return:
        list of lines
    '''
    with open(filename, 'r', encoding='utf-8') as f:
        contents = f.readlines()
    return contents


def zeo_calculation(ase_atom, probe_radius=1.8, number_of_steps=5000, high_accuracy=True):
    '''
    Main script to compute geometric structure of porous systems.
    The focus here is on MOF, but the script can run on any porous periodic
    system. The script computes the accesible surface area, accessible volume
    and the pore geometry. There are many more outputs which can be extracted
    from ,vol_str and sa_str. Moreover there are also other computation that can be done.
    Check out the test directory in dependencies/pyzeo/test. Else contact bafgreat@gmail.com
    if you need more output and can't figure it out.
    Main parameter:
        ase_atom: ase atom object
        probe_radius: The radius of the probe. Here 1.86 is used as default
        number_of_steps: Number of GCMC simulation cycles
        high_accuracy: key to determine where to perform high accuracy computation
    return
        python dictionary containing
        1) Accessible volume void fraction
        2) Accessible volume (A^3)
        3) Accessible surface area (A^2)
        4) Number_of_channels: Number of channels present in the porous system, which correspond to the number of
                               pores within the system
        5) LCD_A: The largest cavity diameter is the largest sphere that can be inserted in a porous
                  system without overlapping with any of the atoms in the system.
        6) lfpd_A:The largest included sphere along free sphere path is
                  largest sphere that can be inserted in the pore
        7) PLD_A:The pore limiting diameter is the largest sphere that can freely
                 diffuse through the porous network without overlapping with any
                 of the atoms in the system
    '''
    tmp_cssr = 'tmp.cssr'
    tmp_out = 'tmp.res'
    tmp = ase_to_zeoobject(ase_atom)
    put_contents(tmp_cssr, tmp)
    parameters = {}
    atmnet = AtomNetwork.read_from_CSSR(tmp_cssr)
    vol_str = volume(
        atmnet, probe_radius, probe_radius, number_of_steps, high_accuracy=high_accuracy)
    if high_accuracy is True:
        vol_str = vol_str[0].decode("utf-8").split()
    else:
        vol_str = vol_str.decode("utf-8").split()
    parameters['AV_Volume_fraction'] = np.float64(vol_str[10])
    parameters['AV'] = np.float64(vol_str[8])
    sa_str = surface_area(atmnet, probe_radius, probe_radius,
                          number_of_steps, high_accuracy=False)
    sa_str = sa_str.decode("utf-8").split()
    parameters['ASA'] = np.float64(sa_str[8])
    parameters['Number_of_channels'] = np.int64(sa_str[20])
    atmnet.calculate_free_sphere_parameters(tmp_out)
    outlines = get_contents(tmp_out)
    data = outlines[0].split()
    parameters['LCD'] = np.float64(data[1])
    parameters['lfpd'] = np.float64(data[3])
    parameters['PLD'] = np.float64(data[2])
    if os.path.exists(tmp_cssr):
        os.remove(tmp_cssr)
    if os.path.exists(tmp_out):
        os.remove(tmp_out)
    return parameters


def ase_to_zeoobject(ase_atom):
    '''
    Converts an ase atom type to a zeo++ Cssr object
    In zeo++ the xyz coordinate system is rotated to a zyx format.

    Args:
        ase atom
    Return:
        zeo++ cssr object
    '''
    pymol = AseAtomsAdaptor.get_structure(ase_atom)
    a, b, c = ase_atom.cell.lengths()
    alpha, beta, gama = ase_atom.cell.angles()
    load = [
        f"{c:.4f} {b:.4f} {a:.4f}",
        f"{gama:.2f} {beta:.2f} {alpha:.2f} SPGR =  1 P 1    OPT = 1",
        f"{len(ase_atom)} 0",
        f"{pymol.formula}"
    ]
    for index, atom in enumerate(ase_atom):
        charge = pymol[index].charge if hasattr(pymol[index], "charge") else 0
        element = atom.symbol
        position = ase_atom.get_scaled_positions()[index]
        load.append(
            f"{index+1} {element} { position[2]:.4f} {position[1]:.4f}  {position[0]:.4f} 0 0 0 0 0 0 0 0 {charge:.4f}")
    return "\n".join(load)
