# How to work with the H5MD-NOMAD schema

!!! warning "Attention"

    The H5MD-NOMAD functionalities are still in the beta testing phase.
    Please contact us with any issues or questions.

## Writing an HDF5 file according to H5MD-NOMAD with python

You can write to an HDF5 file via a python interface, using the [h5py](https://docs.h5py.org/en/stable/quick.html) package. This section walks you through the creation of each section of the H5MD-NOMAD schema, using practical examples to help you get started.

### Imports
```python
import numpy as np
import json

import h5py
import parmed as chem
import MDAnalysis as mda
from pint import UnitRegistry
ureg = UnitRegistry()
```

[h5py](https://docs.h5py.org/en/stable/quick.html)
: module for reading and writing HDF5 files.

[UnitRegistry](https://pint.readthedocs.io/en/0.10.1/tutorial.html)
: object from the [pint](https://pint.readthedocs.io/en/stable/) package that provides assistance for working with units. We suggest using this package for easiest compatibility with NOMAD. If you have `nomad-lab` installed, you can alternatively import `ureg` with `from nomad.units import ureg`.

[MDAnalysis](https://www.mdanalysis.org/)
: a library to analyze trajectories from molecular dynamics simulations stored in various formats.

[ParmEd](https://parmed.github.io/ParmEd/html/index.html)
: a tool for aiding in investigations of biomolecular systems using popular molecular simulation packages.

### Example Data

For concreteness, we consider a fictitious set of "vanilla" molecular dynamics simulations, run with the [OpenMM](https://openmm.org/) software. The following definitions set the dimensionality, periodicity, and the units for this simulation.

```python
dimension = 3
periodicity = [True, True, True]

time_unit = 1.0 * ureg.picosecond
length_unit = 1.0 * ureg.angstrom
energy_unit = 1000. * ureg.joule
mass_unit = 1.0 * ureg.amu
charge_unit = 1.0 * ureg.e
temperature_unit = 1.0 * ureg.K
custom_unit = 1.0 * ureg.newton / length_unit**2
acceleration_unit = 1.0 * length_unit / time_unit**2
```

In this example, we will assume that the relevant simulation data is compatible with [MDAnalysis](https://www.mdanalysis.org/), such that a **universe** containing the trajectory and topology information can be created.

!!! Note

    Knowledge of the MDAnalysis package is not necessary for understanding this example. The dimensions of the supplied quantities will be made clear in each case.

Create a universe by supplying a `pdb` structure file and corresponding `dcd` trajectory file ([MDAnalysis supports many different file formats](https://userguide.mdanalysis.org/stable/formats/index.html)):
```python
universe = mda.Universe('initial_structure.pdb', 'trajectory.dcd')

n_frames = len(universe.trajectory)
n_atoms = universe.trajectory[0].n_atoms
```
Some topologies can be loaded directly into MDAnalysis. However, for simulations from OpenMM, one can read the topology using `parmed` and then import it to MDanalysis:
```python
pdb = app.PDBFile('initial_structure.pdb')
forcefield = app.ForceField('force_field.xml')
system = forcefield.createSystem(pdb.topology)
struct = chem.openmm.load_topology(pdb.topology, system)
universe_toponly = mda.Universe(struct)
```

### [H5MD Group](h5md_expl.md#the-h5md-group)

Create an HDF5 file called `test_h5md-nomad.h5` and create the group `h5md` under `root`:
```python
h5_write = h5py.File('test_h5md-nomad.h5', 'w')
h5md = h5_write.create_group('h5md')
```

Add the h5md version (1.0.x in this case) as an attribute of the `h5md` group:
```python
h5md.attrs['version'] = [1, 0]
```

Create the `author` group and add the associated metadata:
```python
author = h5md.create_group('author')
author.attrs['name'] = 'author name'
author.attrs['email'] = 'author-name@example-domain.edu'
```

Create the `program` group and add the associated metadata:
```python
program = h5md.create_group('program')
program.attrs['name'] = 'OpenMM'
program.attrs['version'] = '7.7.0'
```

Create the `creator` group and add the associated metadata:
```python
program = h5md.create_group('creator')
program.attrs['name'] = h5py.__name__
program.attrs['version'] = str(h5py.__version__)
```

### [Particles Group](h5md_expl.md#the-particles-group)

Create the `particles` group and the underlying `all` group to hold the relevant particle data:
```python
particles = h5_write.create_group('particles')
particles_group_all = particles.create_group('all')
```

Get the steps, times, positions, and lattice vectors (i.e., box dimensions) from the MDA universe:
```python
# quantities extracted from MDAnalysis
steps = []
times = []
positions = []
lattice_vectors = []
for i_frame, frame in enumerate(universe.trajectory):
    times.append(frame.time)
    steps.append(frame.frame)
    positions.append(frame.positions)
    lattice_vectors.append(frame.triclinic_dimensions)
```

Set the positions and corresponding metadata:
```python
position_group_all = particles_group_all.create_group('position')
position_group_all['step'] = steps  # shape = (n_frames)
position_group_all['time'] = times  # shape = (n_frames)
position_group_all['time'].attrs['unit'] = str(time_unit.units)
position_group_all['time'].attrs['unit_factor'] = time_unit.magnitude
position_group_all['value'] = positions  # shape = (n_frames, n_atoms, dimension)
position_group_all['value'].attrs['unit'] = str(length_unit.units)
position_group_all['value'].attrs['unit_factor'] = length_unit.magnitude
```

Set the particle-specific metadata:
```python
particles_group_all['species_label'] = universe_toponly.atoms.types  # shape = (n_atoms)
particles_group_all['force_field_label'] = universe_toponly.atoms.names  # shape = (n_atoms)
particles_group_all['mass'] = universe_toponly.atoms.masses  # shape = (n_atoms)
particles_group_all['mass'].attrs['unit'] = str(mass_unit.units)
particles_group_all['mass'].attrs['unit_factor'] = mass_unit.magnitude
particles_group_all['charge'] = universe_toponly.atoms.charges  # shape = (n_atoms)
particles_group_all['charge'].attrs['unit'] = str(charge_unit.units)
particles_group_all['charge'].attrs['unit_factor'] = charge_unit.magnitude
```
<!-- TODO - Should I say something about the species label here? Original code was:  [re.sub('[1-9]', '', atom_type) for atom_type in universe_toponly.atoms.types] -->

Create the `box` group under `particles.all` and write corresponding data:
```python
box_group = particles_group_all.create_group('box')
box_group.attrs['dimension'] = dimension
box_group.attrs['boundary'] = periodicity
edges = box_group.create_group('edges')
edges['step'] = steps
edges['time'] = times
edges['time'].attrs['unit'] = str(time_unit.units)
edges['time'].attrs['unit_factor'] = time_unit.magnitude
edges['value'] = lattice_vectors
edges['value'].attrs['unit'] = str(length_unit.units)
edges['value'].attrs['unit_factor'] = length_unit.magnitude

```

### [Connectivity Group](h5md_expl.md#the-connectivity-group)

Create the `connectivity` group under `root` and add the tuples of bonds, angles, and dihedrals:
```python
connectivity = h5_write.create_group('connectivity')
connectivity['bonds'] = universe_toponly.bonds._bix  # shape = (n_bonds, 2)
connectivity['angles'] = universe_toponly.angles._bix  # shape = (n_angles, 3)
connectivity['dihedrals'] = universe_toponly.dihedrals._bix  # shape = (n_dihedrals, 4)
connectivity['impropers'] = universe_toponly.impropers._bix  # shape = (n_impropers, 4)
```
Here `n_bonds`, `n_angles`, `n_dihedrals`, and `n_impropers` represent the corresponding number of instances of each interaction within the force field.

You can read more about the creation of the hierarchical `particles_group` in [Creating a topology](#creating-a-topology-particles_group).


### [Observables Group](h5md_expl.md#the-observables-group)

For this section, we will consider sets of fabricated observable data for clarity. First, create the `observables` group under root:
```python
observables = h5_write.create_group('observables')
```

There are 3 types of support observables:
```python
types = ['configurational', 'ensemble_average', 'correlation_function']
```

#### [Configurational Observables](h5md_expl.md#configurational_observable_anchor)

Fabricated data:
```python
temperatures = 300. * np.ones(n_frames)
potential_energies = 1.0 * np.ones(n_frames)
kinetic_energies = 2.0 * np.ones(n_frames)
```

Create a `temperature` group and populate the associated metadata:

```python
temperature = observables.create_group('temperature')
temperature.attrs['type'] = types[0]
temperature['step'] = steps
temperature['time'] = times
temperature['time'].attrs['unit'] = str(time_unit.units)
temperature['time'].attrs['unit_factor'] = time_unit.magnitude
temperature['value'] = temperatures
temperature['value'].attrs['unit'] = str(temperature_unit.units)
temperature['value'].attrs['unit_factor'] = temperature_unit.magnitude
```

Create an `energy` group to hold various types of energies. Add :
```python
energies = observables.create_group('energy')

potential_energy = energies.create_group('potential')
potential_energy.attrs['type'] = types[0]
potential_energy['step'] = steps
potential_energy['time'] = times
potential_energy['time'].attrs['unit'] = str(time_unit.units)
potential_energy['time'].attrs['unit_factor'] = time_unit.magnitude
potential_energy['value'] = potential_energies
potential_energy['value'].attrs['unit'] = str(energy_unit.units)
potential_energy['value'].attrs['unit_factor'] = energy_unit.magnitude

kinetic_energy = energies.create_group('kinetic')
kinetic_energy.attrs['type'] = types[0]
kinetic_energy['step'] = steps
kinetic_energy['time'] = times
kinetic_energy['time'].attrs['unit'] = str(time_unit.units)
kinetic_energy['time'].attrs['unit_factor'] = time_unit.magnitude
kinetic_energy['value'] = kinetic_energies
kinetic_energy['value'].attrs['unit'] = str(energy_unit.units)
kinetic_energy['value'].attrs['unit_factor'] = energy_unit.magnitude
```
#### [Ensemble Average Observables](h5md_expl.md#ensemble_average_observable_anchor)

Fabricated data - the following represents radial distribution function (rdf) data calculated between molecule types `X` and `Y`, stored in `rdf_MOLX-MOLY.xvg`:
```
      0.24 0.000152428
     0.245 0.00457094
      0.25  0.0573499
     0.255   0.284764
      0.26   0.842825
     0.265    1.64705
      0.27    2.37243
     0.275    2.77916
      0.28    2.80622
     0.285    2.60082
      0.29    2.27182
      ...
```

Store the rdf data in a dictionary along with some relevant metadata:

```python
rdf_XX = np.loadtxt('rdf_MOLX-MOLX.xvg')
rdf_XY = np.loadtxt('rdf_MOLX-MOLY.xvg')
rdf_YY = np.loadtxt('rdf_MOLY-MOLY.xvg')
rdfs = {
    'MOLX-MOLX': {
        'n_bins': len(rdf_XX[:, 0]),
        'bins': rdf_XX[:, 0],
        'value': rdf_XX[:, 1],
        'type': 'molecular',
        'frame_start': 0,
        'frame_end': n_frames-1
    },
    'MOLX-MOLY': {
        'n_bins': len(rdf_XY[:, 0]),
        'bins': rdf_XY[:, 0],
        'value': rdf_XY[:, 1],
        'type': 'molecular',
        'frame_start': 0,
        'frame_end': n_frames-1
    },
    'MOLY-MOLY': {
        'n_bins': len(rdf_YY[:, 0]),
        'bins': rdf_YY[:, 0],
        'value': rdf_YY[:, 1],
        'type': 'molecular',
        'frame_start': 0,
        'frame_end': n_frames-1
    }
}
```

Now create the `radial_distribution_functions` group under `observables` and store each imported rdf:
```python
radial_distribution_functions = observables.create_group('radial_distribution_functions')
for key in rdfs.keys():
    rdf = radial_distribution_functions.create_group(key)
    rdf.attrs['type'] = types[1]
    rdf['type'] = rdfs[key]['type']
    rdf['n_bins'] = rdfs[key]['n_bins']
    rdf['bins'] = rdfs[key]['bins']
    rdf['bins'].attrs['unit'] = str(length_unit.units)
    rdf['bins'].attrs['unit_factor'] = length_unit.magnitude
    rdf['value'] = rdfs[key]['value']
    rdf['frame_start'] = rdfs[key]['frame_start']
    rdf['frame_end'] = rdfs[key]['frame_end']
```

We can also store scalar ensemble average observables. Let's consider some fabricated diffusion constant data:
```python
Ds = {
    'MOLX': {
        'value': 1.0,
        'error_type': 'Pearson_correlation_coefficient',
        'errors': 0.98
    },
    'MOLY': {
        'value': 2.0,
        'error_type': 'Pearson_correlation_coefficient',
        'errors': 0.95
    }
}
```

Create the `diffusion constants` group under `observables` and store the correspond (meta)data:

```python
diffusion_constants = observables.create_group('diffusion_constants')
for key in Ds.keys():
    diffusion_constant = diffusion_constants.create_group(key)
    diffusion_constant.attrs['type'] = types[1]
    diffusion_constant['value'] = Ds[key]['value']
    diffusion_constant['value'].attrs['unit'] = str(diff_unit.units)
    diffusion_constant['value'].attrs['unit_factor'] = diff_unit.magnitude
    diffusion_constant['error_type'] = Ds[key]['error_type']
    diffusion_constant['errors'] = Ds[key]['errors']
```

#### [Time Correlation Observables](h5md_expl.md#time_correlation_observable_anchor)

Fabricated data - the following represents mean squared displacement (msd) data calculated for molecule type `X`, stored in `msd_MOLX.xvg`:
```
         0           0
         2   0.0688769
         4    0.135904
         6    0.203573
         8    0.271162
        10    0.339284
        12    0.410115
        14    0.477376
        16    0.545184
        18     0.61283
        ...
```

Store the msd data in a dictionary along with some relevant metadata:

```python
msd_X = np.loadtxt('msd_MOLX.xvg')
msd_Y = np.loadtxt('msd_MOLY.xvg')
msds = {
    'MOLX': {
        'n_times': len(msd_X[:, 0]),
        'times': msd_X[:, 0],
        'value': msd_X[:, 1],
        'type': 'molecular',
        'direction': 'xyz',
        'error_type': 'standard_deviation',
        'errors': np.zeros(len(msd_X[:, 0])),
    },
    'MOLY': {
        'n_times': len(msd_Y[:, 0]),
        'times': msd_Y[:, 0],
        'value': msd_Y[:, 1],
        'type': 'molecular',
        'direction': 'xyz',
        'error_type': 'standard_deviation',
        'errors': np.zeros(len(msd_Y[:, 0])),
    }
}
```

Now create the `mean_squared_displacements` group under `observables` and store each imported rdf:

```python
mean_squared_displacements = observables.create_group('mean_squared_displacements')
msd_unit = length_unit * length_unit
diff_unit = msd_unit / time_unit
for key in msds.keys():
    msd = mean_squared_displacements.create_group(key)
    msd.attrs['type'] = types[2]
    msd['type'] = msds[key]['type']
    msd['direction'] = msds[key]['direction']
    msd['error_type'] = msds[key]['error_type']
    msd['n_times'] = msds[key]['n_times']
    msd['times'] = msds[key]['times']
    msd['times'].attrs['unit'] = str(time_unit.units)
    msd['times'].attrs['unit_factor'] = time_unit.magnitude
    msd['value'] = msds[key]['value']
    msd['value'].attrs['unit'] = str(msd_unit.units)
    msd['value'].attrs['unit_factor'] = msd_unit.magnitude
    msd['errors'] = msds[key]['errors']
```

### [Parameter Group](h5md_expl.md#the-parameters-group)

Using the json templates for [force calculations](./h5md_expl.md#force_calculation_template_anchor) and [molecular dynamics workflows](h5md_expl.md#md_workflow_template_anchor), the (meta)data can be written to the H5MD-NOMAD file using the following code:

First, import the data extracted from the JSON templates:
```python
with open('force_calculations_metainfo.json') as json_file:
    force_calculation_parameters = json.load(json_file)

with open('workflow_metainfo.json') as json_file:
    workflow_parameters = json.load(json_file)
```

Then, create the appropriate container groups:
```python
parameters = h5_write.create_group('parameters')
force_calculations = parameters.create_group('force_calculations')
workflow = parameters.create_group('workflow')
```


Now, recursively write the (meta)data:
```python
def get_parameters_recursive(parameter_group, parameter_dict):
    # Store the parameters from parameter dict into an hdf5 file
    for key, val in parameter_dict.items():
        if type(val) == dict:
            param = val.get('value')
            if param is not None:
                parameter_group[key] = param
                unit = val.get('unit')
                if unit is not None:
                    parameter_group[key].attrs['unit'] = unit
            else:  # This is a subsection
                subsection = parameter_group.create_group(key)
                subsection = get_parameters_recursive(subsection, val)
        else:
            if val is not None:
                parameter_group[key] = val

    return parameter_group


force_calculations = get_parameters_recursive(force_calculations, force_calculation_parameters)
workflow = get_parameters_recursive(workflow, workflow_parameters)
```

It's as simple as that! Now, we can [upload our H5MD-NOMAD file directly to NOMAD](./uploading.md#drag-and-drop-uploads) and all the written (meta)data will be stored according to the standard NOMAD schema.



## Accessing an H5MD-NOMAD file


The following functions are useful for accessing data from your H5MD-NOMAD file:
```python
def apply_unit(quantity, unit, unit_factor):
    from pint import UnitRegistry
    ureg = UnitRegistry()

    if quantity is None:
        return
    if unit:
        unit = ureg(unit)
        unit *= unit_factor
        quantity *= unit

    return quantity

def decode_hdf5_bytes(dataset):
    if dataset is None:
        return
    elif type(dataset).__name__ == 'ndarray':
        if dataset == []:
            return dataset
        dataset = np.array([val.decode("utf-8") for val in dataset]) if type(dataset[0]) == bytes else dataset
    else:
        dataset = dataset.decode("utf-8") if type(dataset) == bytes else dataset
    return dataset

def hdf5_attr_getter(source, path, attribute, default=None):
    '''
    Extracts attribute from object based on path, and returns default if not defined.
    '''
    section_segments = path.split('.')
    for section in section_segments:
        try:
            value = source.get(section)
            source = value[-1] if isinstance(value, list) else value
        except Exception:
            return
    value = source.attrs.get(attribute)
    source = value[-1] if isinstance(value, list) else value
    source = decode_hdf5_bytes(source) if source is not None else default
    return source

def hdf5_getter(source, path, default=None):
    '''
    Extracts attribute from object based on path, and returns default if not defined.
    '''
    section_segments = path.split('.')
    for section in section_segments:
        try:
            value = source.get(section)
            unit = hdf5_attr_getter(source, section, 'unit')
            unit_factor = hdf5_attr_getter(source, section, 'unit_factor', default=1.0)
            source = value[-1] if isinstance(value, list) else value
        except Exception:
            return

    if source is None:
        source = default
    elif type(source) == h5py.Dataset:
        source = source[()]
        source = apply_unit(source, unit, unit_factor)

    return decode_hdf5_bytes(source)
```

Open your H5MD-NOMAD file with `h5py`:
```python
import h5py

h5_read = h5py.File('test_h5md-nomad.h5', 'r')
```

Access a particular data set:
```python
potential_energies = h5_read['observables']['energy']['potential']['value']
print(potential_energies[()])
```
result:
```
array([1., 1., 1., 1., 1.])
```

Get the unit information for this quantity:
```python
unit = potential_energies.attrs['unit']
unit_factor = potential_energies.attrs['unit_factor']

print(unit)
print(unit_factor)
```

results:
```
joule
1000.0
```

Alternatively, the above functions will return the dataset as python arrays, i.e., already applying `[()]` to the HDF5 element, and also apply the appropriate units where applicable:
```python
potential_energies = hdf5_getter(h5_read, 'observables.energy.potential.value')
print(potential_energies)
```

result:
```
Magnitude
[1000.0 1000.0 1000.0 1000.0 1000.0]
Units	joule
```


<!-- TODO Separate topology creation and put into the general docs -->
## Creating a topology (`particles_group`)

This page demonstrates how to create a "standard" topology in H5MD-NOMAD. The demonstrated organization of molecules and monomers is identical to what other NOMAD parsers do to create a topology from native simulation files (e.g., outputs from [GROMACS](https://www.gromacs.org/) or [LAMMPS](https://www.lammps.org)). However, the user is free to deviate from this standard to create arbitrary organizations of particles, as described in [Connectivity](h5md_expl.md#the-connectivity-group).
<!-- TODO add hrefs to "writing a parser plugin" section when you say "NOMAD parsers" -->

### Standard topology structure for bonded force fields

```
topology
├── molecule_group_1
│    ├── molecule_1
│    │      ├── monomer_group_1
│    │      │       ├── monomer_1
│    │      │       │       └── metadata for monomer_1
│    │      │       ├── monomer_2
│    │      │       │       └── metadata for monomer_2
│    │      │       ├── ...
│    │      ├── monomer_group_2
│    │      └── ...
│    ├── molecule_2
│    │      ├── ...
│    └── ...
└── molecule_group_2
│    ├── molecule_1
│    │      └── metadata for molecule_1
│    ├── molecule_2
│    │      └── metadata for molecule_2
│    │      └── ...
│    └── ...
└── ...
```
Here, the first level of organization is the "molecule group". Molecule groups contain molecules of the same type. In other words, `molecule_group_1` and `molecule_group_2` represent distinct molecule types. At the next level of the hierarchy, each molecule within this group is stored (i.e., `molecule_1`, `molecule_2`, etc.). In the above example, `molecule_group_1` represents a polymer (or protein). Thus, below the molecule level, there is a "monomer group level". Similar to the molecule group, the monomer group organizes all monomers (of the parent molecule) that are of the same type. Thus, for `molecule_1` of `molecule_group_1`, `monomer_group_1` and `monomer_group_2` represent distinct types of monomers existing within the polymer. Then, below `monomer_group_1`, each monomer within this group is stored. Finally, beneath these individual monomers, only the metadata for that monomer is stored (i.e., no further organization levels). Note however, that metadata can be (and is) stored at each level of the hierarchy, but is left out of the illustration for clarity. Notice also that `molecule_group_2` is not a polymer. Thus, each molecule within this group stores only the corresponding metadata, and no further levels of organization.
<!-- TODO - add image of topology bar in overview page? Here or maybe when the topology is introduced on the connectivity page?
add full explanation and reference about what Topology/Groups mean in NOMAD -->

### Creating the standard hierarchy from an MDAnalysis universe
We start from the perspective of the [Writing an HDF5 file according to H5MD-NOMAD with python](#writing-an-hdf5-file-according-to-h5md-nomad-with-python) section, with identical imports and assuming that an MDAnalysis `universe` is already instantiated from the raw simulation files. As in the previous example, the `universe` containing the topology information is called `universe_topology`.

The following functions will be useful for creating the topology:

```python
def get_composition(children_names):
    '''
    Given a list of children, return a compositional formula as a function of
    these children. The format is <child_1>(n_child_1)<child_2>(n_child_2)...
    '''
    children_count_tup = np.unique(children_names, return_counts=True)
    formula = ''.join([f'{name}({count})' for name, count in zip(*children_count_tup)])
    return formula


def get_molecules_from_bond_list(n_particles: int, bond_list: List[int], particle_types: List[str] = None, particles_typeid=None):
    '''
    Returns a dictionary with molecule info from the list of bonds
    '''

    import networkx

    system_graph = networkx.empty_graph(n_particles)
    system_graph.add_edges_from([(i[0], i[1]) for i in bond_list])
    molecules = [system_graph.subgraph(c).copy() for c in networkx.connected_components(system_graph)]
    mol_dict = []
    for i_mol, mol in enumerate(molecules):
        mol_dict.append({})
        mol_dict[i_mol]['indices'] = np.array(mol.nodes())
        mol_dict[i_mol]['bonds'] = np.array(mol.edges())
        mol_dict[i_mol]['type'] = 'molecule'
        mol_dict[i_mol]['is_molecule'] = True
        if particles_typeid is None and len(particle_types) == n_particles:
            mol_dict[i_mol]['names'] = [particle_types[int(x)] for x in sorted(np.array(mol_dict[i_mol]['indices']))]
        if particle_types is not None and particles_typeid is not None:
            mol_dict[i_mol]['names'] = [particle_types[particles_typeid[int(x)]] for x in sorted(np.array(mol_dict[i_mol]['indices']))]
        mol_dict[i_mol]['formula'] = get_composition(mol_dict[i_mol]['names'])

    return mol_dict


def is_same_molecule(mol_1: dict, mol_2: dict):
    '''
    Checks whether the 2 input molecule dictionaries represent the same
    molecule type, i.e., same particle types and corresponding bond connections.
    '''

    if sorted(mol_1['names']) == sorted(mol_2['names']):
        mol_1_shift = np.min(mol_1['indices'])
        mol_2_shift = np.min(mol_2['indices'])
        mol_1_bonds_shift = mol_1['bonds'] - mol_1_shift
        mol_2_bonds_shift = mol_2['bonds'] - mol_2_shift

        bond_list_1 = [sorted((mol_1['names'][i], mol_1['names'][j])) for i, j in mol_1_bonds_shift]
        bond_list_2 = [sorted((mol_2['names'][i], mol_2['names'][j])) for i, j in mol_2_bonds_shift]

        bond_list_names_1, bond_list_counts_1 = np.unique(bond_list_1, axis=0, return_counts=True)
        bond_list_names_2, bond_list_counts_2 = np.unique(bond_list_2, axis=0, return_counts=True)

        bond_list_dict_1 = {bond[0] + '-' + bond[1]: bond_list_counts_1[i_bond] for i_bond, bond in enumerate(bond_list_names_1)}
        bond_list_dict_2 = {bond[0] + '-' + bond[1]: bond_list_counts_2[i_bond] for i_bond, bond in enumerate(bond_list_names_2)}
        if bond_list_dict_1 == bond_list_dict_2:
            return True

        return False

    return False
```

Then, we can create the topology structure from the MDAnalysis universe:

```python
bond_list = universe_toponly.bonds._bix
molecules = get_molecules_from_bond_list(n_atoms, bond_list, particle_types=universe_toponly.atoms.types, particles_typeid=None)

# create the topology
mol_groups = []
mol_groups.append({})
mol_groups[0]['molecules'] = []
mol_groups[0]['molecules'].append(molecules[0])
mol_groups[0]['type'] = 'molecule_group'
mol_groups[0]['is_molecule'] = False
for mol in molecules[1:]:
    flag_mol_group_exists = False
    for i_mol_group in range(len(mol_groups)):
        if is_same_molecule(mol, mol_groups[i_mol_group]['molecules'][0]):
            mol_groups[i_mol_group]['molecules'].append(mol)
            flag_mol_group_exists = True
            break
    if not flag_mol_group_exists:
        mol_groups.append({})
        mol_groups[-1]['molecules'] = []
        mol_groups[-1]['molecules'].append(mol)
        mol_groups[-1]['type'] = 'molecule_group'
        mol_groups[-1]['is_molecule'] = False


for i_mol_group, mol_group in enumerate(mol_groups):
    mol_groups[i_mol_group]['formula'] = molecule_labels[i_mol_group] + '(' + str(len(mol_group['molecules'])) + ')'
    mol_groups[i_mol_group]['label'] = 'group_' + str(molecule_labels[i_mol_group])
    mol_group_indices = []
    for i_molecule, molecule in enumerate(mol_group['molecules']):
        molecule['label'] = molecule_labels[i_mol_group]
        mol_indices = molecule['indices']
        mol_group_indices.append(mol_indices)
        mol_resids = np.unique(universe_toponly.atoms.resindices[mol_indices])
        if mol_resids.shape[0] == 1:
            continue

        res_dict = []
        for i_resid, resid in enumerate(mol_resids):
            res_dict.append({})
            res_dict[i_resid]['indices'] = np.where( universe_toponly.atoms.resindices[mol_indices] == resid)[0]
            res_dict[i_resid]['label'] = universe_toponly.atoms.resnames[res_dict[i_resid]['indices'][0]]
            res_dict[i_resid]['formula'] = get_composition(universe_toponly.atoms.names[res_dict[i_resid]['indices']])
            res_dict[i_resid]['is_molecule'] = False
            res_dict[i_resid]['type'] = 'monomer'

        res_groups = []
        res_groups.append({})
        res_groups[0]['residues'] = []
        res_groups[0]['residues'].append(res_dict[0])
        res_groups[0]['label'] = 'group_' + res_dict[0]['label']
        res_groups[0]['type'] = 'monomer_group'
        res_groups[0]['is_molecule'] = False
        for res in res_dict[1:]:
            flag_res_group_exists = False
            for i_res_group in range(len(res_groups)):
                if res['label'] == res_groups[i_res_group]['label']:
                    res_groups[i_res_group]['residues'].append(res)
                    flag_res_group_exists = True
                    break
            if not flag_res_group_exists:
                res_groups.append({})
                res_groups[-1]['residues'] = []
                res_groups[-1]['residues'].append(res)
                res_groups[-1]['label'] = 'group_' + res['label']
                res_groups[-1]['formula'] = get_composition(universe_toponly.atoms.names[res['indices']])
                res_groups[-1]['type'] = 'monomer_group'
                res_groups[-1]['is_molecule'] = False

        molecule['formula'] = ''
        for res_group in res_groups:
            res_group['formula'] = res_group['residues'][0]['label'] + '(' + str(len(res_group['residues'])) + ')'
            molecule['formula'] += res_group['formula']
            res_group_indices = []
            for res in res_group['residues']:
                res_group_indices.append(res['indices'])
            res_group['indices'] = np.concatenate(res_group_indices)

        mol_group['indices'] = np.concatenate(mol_group_indices)

        molecule['residue_groups'] = res_groups

```

### Writing the topology to an H5MD-NOMAD file

Here we assume an H5MD-NOMAD file has already been created, as demonstrated on the [Writing an HDF5 file according to H5MD-NOMAD with python](#writing-an-hdf5-file-according-to-h5md-nomad-with-python) section, and that the `connectivity` group was created under the root level.

Now, create the `particles_group` group under `connectivity` within our HDF5-NOMAD file:
```python
topology_keys = ['type', 'formula', 'particles_group', 'label', 'is_molecule', 'indices']
custom_keys = ['molecules', 'residue_groups', 'residues']

topology = connectivity.create_group('particles_group')

for i_mol_group, mol_group in enumerate(mol_groups):
    hdf5_mol_group = topology.create_group('group_' + molecule_labels[i_mol_group])
    for mol_group_key in mol_group.keys():
        if mol_group_key not in topology_keys + custom_keys:
            continue
        if mol_group_key != 'molecules':
            hdf5_mol_group[mol_group_key] = mol_group[mol_group_key]
        else:
            hdf5_molecules = hdf5_mol_group.create_group('particles_group')
            for i_molecule, molecule in enumerate(mol_group[mol_group_key]):
                hdf5_mol = hdf5_molecules.create_group('molecule_' + str(i_molecule))
                for mol_key in molecule.keys():
                    if mol_key not in topology_keys + custom_keys:
                        continue
                    if mol_key != 'residue_groups':
                        hdf5_mol[mol_key] = molecule[mol_key]
                    else:
                        hdf5_residue_groups = hdf5_mol.create_group('particles_group')
                        for i_res_group, res_group in enumerate(molecule[mol_key]):
                            hdf5_res_group = hdf5_residue_groups.create_group('residue_group_' + str(i_res_group))
                            for res_group_key in res_group.keys():
                                if res_group_key not in topology_keys + custom_keys:
                                    continue
                                if res_group_key != 'residues':
                                    hdf5_res_group[res_group_key] = res_group[res_group_key]
                                else:
                                    hdf5_residues = hdf5_res_group.create_group('particles_group')
                                    for i_res, res in enumerate(res_group[res_group_key]):
                                        hdf5_res = hdf5_residues.create_group('residue_' + str(i_res))
                                        for res_key in res.keys():
                                            if res_key not in topology_keys:
                                                continue
                                            if res[res_key] is not None:
                                                hdf5_res[res_key] = res[res_key]
```


