# H5MD-NOMAD: A flexible data-storage schema for uploading molecular simulations to NOMAD

## Overview

Most computational data in NOMAD is harvested with code-specific parsers that recognize the output files from a particular software and retrieve the appropriate (meta)data accordingly.
However, this approach is not possible for many modern molecular simulation engines that use fully-flexible scriptable input and non-fixed output files.
["HDF5 for molecular data" (H5MD)](http://h5md.nongnu.org/) is a data schema for storage of molecular simulation data, based on the HDF5 file format.
This page describes an extension of the H5MD schema, denoted H5MD-NOMAD, which adds specificity to several of the H5MD guidelines while also retaining reasonable flexibility. This enables simulation data stored according to the H5MD-NOMAD schema to be stored in the NOMAD.

**Due to the nature of extending upon the original H5MD schema, portions of this doc page was duplicated, extended, or summarized from the [H5MD webpage](http://h5md.nongnu.org/).**

## Introduction to the H5MD storage format

H5MD was originally proposed by P. de Buyl, P. H. Colberg and F. Höfling in [H5MD: A structured, efficient, and portable file format for molecular data](http://dx.doi.org/10.1016/j.cpc.2014.01.018), Comp. Phys. Comm. 185, 1546–1553 (2014) [[arXiv:1308.6382](http://arxiv.org/abs/1308.6382)]. The schema is maintained, along with associated tools, in a GitHub repository: [H5MD GitHub](https://github.com/h5md).

This section provides the basic nomenclature of the H5MD schema relevant for understanding H5MD-NOMAD, and was duplicated or summarized from the [H5MD webpage](http://h5md.nongnu.org/).

### File format

H5MD structures are stored in the
[HDF5 file format](http://www.hdfgroup.org/HDF5/doc/H5.format.html) version 0
or later. It is recommended to use the HDF5 file format version 2, which
includes the implicit tracking of the creation and modification times of the
file and of each of its objects.

### Notation and naming

HDF5 files are organized into groups and datasets, summarized as *objects*,
which form a tree structure with the datasets as leaves. Attributes can be
attached to each object. The H5MD specification adopts this naming and uses the
following notation to depict the tree or its subtrees:

`\-- item`
:   An object within a group, that is either a dataset or a group. If it is a
    group itself, the objects within the group are indented by five spaces with
    respect to the group name.

`+-- attribute`
:   An attribute, that relates either to a group or a dataset.

`\-- data: <type>[dim1][dim2]`
:   A dataset with array dimensions `dim1` by `dim2` and of type `<type>`. The
    type is taken from `Enumeration`, `Integer`, `Float` or `String` and follows
    the HDF5 Datatype classes. If the type is not mandated by H5MD, `<type>` is
    indicated. A scalar dataspace is indicated by `[]`.

`(identifier)`
:   An optional item.

`<identifier>`
:   An optional item with unspecified name.

H5MD defines a structure called *H5MD element* (or *element* whenever there
is no confusion). An element is either a time-dependent group or a single
dataset (see time-dependent data below), depending on the situation.

<!-- ### General organization

H5MD defines an organization of the HDF5 file or a part thereof into groups,
datasets, and attributes. The root level of the H5MD structure may coincide
with the root of the HDF5 file or be an arbitrary group inside the HDF5 tree. A
number of groups are defined at the H5MD root level. Several levels of
subgroups may exist inside the H5MD structure, allowing the storage and
description of subsystems.

The H5MD structure is allowed to possess non-specified groups, datasets, or
attributes that contain additional information such as application-specific
parameters or data structures, leaving scope for future extensions. Only the
`h5md` group is mandatory at the H5MD root level. All other root groups are
optional, allowing the user to store only relevant data. Inside each group,
every group or dataset is again optional, unless specified differently.

H5MD supports equally the storage of time-dependent and time-independent data,
i.e., data that change in the course of the simulation or that do not. The
choice between those storage types is not made explicit for the elements in the
specification, it has to be made according to the situation. For instance, the
species and mass of the particles are often fixed in time, but in chemically
reactive systems this might not be appropriate. -->

### Time-dependent data

Time-dependent data consists of a series of samples (or frames) referring to
multiple time steps. Such data are found inside a single dataset and are
accessed via dataset slicing. In order to link the samples to the time axis of
the simulation, H5MD defines a *time-dependent H5MD element* as a group that
contains, in addition to the actual data, information on the corresponding
integer time step and on the physical time. The structure of such a group is:

    <element>
     \-- step
     \-- (time)
     \-- value: <type>[variable][...]

`value`
:   A dataset that holds the data of the time series. It uses a simple
    dataspace whose rank is given by 1 plus the tensor rank of the data stored.
    Its shape is the shape of a single data item prepended by a `[variable]`
    dimension that allows the accumulation of samples during the course of
    time. For instance, the data shape of scalars has the form `[variable]`,
    `D`-dimensional vectors use `[variable][D]`, etc. The first dimension of
    `value` must match the unique dimension of `step` and `time`.

If several H5MD elements are sampled at equal times, `step` and `time` of one
element may be hard links to the `step` and `time` datasets of a different
element. If two elements are sampled at different times (for instance, if one
needs the positions more frequently than the velocities), `step` and `time` are
unique to each of them.

The storage of step and time information follows one of the two modes below,
depending on the dataset layout of `step`.

#### Explicit step and time storage

    <element>
     \-- step: Integer[variable]
     \-- (time: type[variable])
     \-- value: <type>[variable][...]

`step`
:   A dataset with dimensions `[variable]` that contains the time steps at
    which the corresponding data were sampled. It is of `Integer` type to allow
    exact temporal matching of data from one H5MD element to another. The
    values of the dataset are in monotonically increasing order.

`time`
:   An optional dataset that is the same as the `step` dataset, except it is
    `Float` or `Integer`-valued and contains the simulation time in physical units. The
    values of the dataset are in monotonically increasing order.

#### Fixed step and time storage **(currently not supported in H5MD-NOMAD)**

    <element>
     \-- step: Integer[]
	     +-- (offset: type[])
     \-- (time: type[])
	     +-- (offset: type[])
     \-- value: <type>[variable][...]

`step`
:   A scalar dataset of `Integer` type that contains the increment of the
    time step between two successive rows of data in `value`.

    `offset`
	: A scalar attribute of type `Integer` corresponding to the first sampled
    value of `step`.

`time`
:   An optional scalar dataset that is the same as the `step` dataset, except that
    it is `Float` or `Integer`-valued and contains the increment in simulation
    time, in physical units.

`offset`
	: A scalar attribute of the same type as `time` corresponding to the first
    sampled value of `time`.

For this storage mode, the explicit value $s(i)$ of the step corresponding to
the $i$-th row of the dataset `value` is $s(i) = i\times\mathrm{step} +
\mathrm{offset}$ where $\mathrm{offset}$ is set to zero if absent.
The corresponding formula for the time $t(i)$ is identical: $t(i) =
i\times\mathrm{time} + \mathrm{offset}$.
The index $i$ is zero-based.

### Time-independent data

H5MD defines a *time-independent H5MD element* as a dataset. As for the
`value` dataset in the case of time-dependent data, data type and array shape
are implied by the stored data, where the `[variable]` dimension is omitted.

### Storage order of arrays

All arrays are stored in C-order as enforced by the HDF5 file format.
A C or C++ program may thus declare `r[N][D]` for the array
of particle coordinates while the Fortran program will declare a `r(D,N)` array
(appropriate index ordering for a system of `N` particles in `D` spatial
dimensions), and the HDF5 file will be the same.

### Storage of particles and tuples lists

#### Storage of a list of particles

A list of particles is an H5MD element:

    <list_name>: Integer[N]
     +-- particles_group: Object reference

where `list_name` is a dataset of `Integer` type and dimensions `[N]`, N being
the number of particle indices stored in the list. `particles_group` is an
attribute containing an HDF5 Object Reference as defined by the HDF5 file format. `particles_group`
must refer to one of the groups in `/particles`.

If a *fill value* is defined for `list_name`, the particles indices in
`list_name` set to this value are ignored.

If the corresponding `particles_group` does not possess the `id` element, the
values in `list_name` correspond to the indexing of the elements in
`particles_group`. Else, the values in `list_name` must be put in correspondence
with the equal values in the `id` element.

#### Storage of tuples

A list of tuples is an H5MD element:

    <tuples_list_name>: Integer[N,T]
     +-- particles_group: Object reference

where `N` is the length of the list and `T` is the size of the tuples.  Both `N`
and `T` may indicate variable dimensions. `particles_group` is an attribute
containing an HDF5 Object Reference, obeying the same rules as for the lists of
particles.

The interpretation of the values stored within the tuples is done as for a list
of particles.

If a *fill value* is defined, tuples with at least one entry set to this
value are ignored.

#### Time-dependence **(time-dependent particle lists currently not supported in H5MD-NOMAD)**

As the lists of particles and tuples above are H5MD elements, they can be stored
either as time-dependent groups or time-independent datasets.

As an example, a time-dependent list of pairs is stored as:

    <pair_list_name>
       +-- particles_group: Object reference
       \-- value: Integer[variable,N,2]
       \-- step: Integer[variable]

The dimension denoted by `N` may be variable.



## The root level

The root level of H5MD-NOMAD structure is organized as follows (identical to the original H5MD specification):

    H5MD-NOMAD root
     \-- h5md
     \-- (particles)
     \-- (observables)
     \-- (connectivity)
     \-- (parameters)

`h5md`
:   A group that contains metadata and information on the H5MD structure
    itself. It is the only mandatory group at the root level of H5MD.

`particles`
:   An optional group that contains information on each particle in the system,
    e.g., a snapshot of the positions or the full trajectory in phase space.

`observables`
:   An optional group that contains other quantities of interest, e.g.,
    physical observables that are derived from the system state at given points
    in time.

`connectivity`
:   An optional group that contains information about the connectivity between particles.

`parameters`
:   An optional group that contains application-specific (meta)data such as
    control parameters or simulation scripts.



## The H5MD Group

A set of global metadata describing the H5MD structure is stored in the `h5md`
group as attributes. The contents of the group are:

    h5md
     +-- version: Integer[2]
     \-- author
     |    +-- name: String[]
     |    +-- (email: String[])
     \-- creator
     |    +-- name: String[]
     |    +-- version: String[]
     \-- program
          +-- name: String[]
          +-- version: String[]

`version`
:   An attribute, of `Integer` datatype and of simple dataspace of rank 1 and
    size 2, that contains the major version number and the minor version number
    of the H5MD specification the H5MD structure conforms to.

The version *x.y.z* of the H5MD specification follows
[semantic versioning](https://semver.org/spec/v2.0.0.html): A change of the major
version number *x* indicates backward-incompatible changes to the file
structure. A change of the minor version number *y* indicates
backwards-compatible changes to the file structure. A change of the patch
version number *z* indicates changes that have no effect on the file
structure and serves to allow for clarifications or minor text editing of
the specification.

As the *z* component has no impact on the content of an H5MD file, the
`version` attribute contains only *x* and *y*.

`author`
:   A group that contains metadata on the person responsible for the simulation
    (or the experiment) as follows:

* `name`
:   An attribute, of fixed-length string datatype and of scalar
    dataspace, that holds the author's real name.

* `email`
:   An optional attribute, of fixed-length string datatype and
    of scalar dataspace, that holds the author's email address of
    the form `email@domain.tld`.

`creator`
:   A group that contains metadata on the program that created the H5MD
    structure as follows:

* `name`
:   An attribute, of fixed-length string datatype and of scalar
    dataspace, that stores the name of the program.

* `version`
:   An attribute, of fixed-length string datatype and of scalar
    dataspace, that yields the version of the program.

`program`
:   A group that contains metadata on the code/package that created the simulation data contained within this H5MD structure:

* `name`
:   An attribute, of fixed-length string datatype and of scalar
    dataspace, that stores the name of the program.

* `version`
:   An attribute, of fixed-length string datatype and of scalar
    dataspace, that yields the version of the program.


#### Modules **(currently unused in H5MD-NOMAD)**

The original H5MD specification allowed the definition of modules under the h5md group.
Such modules are currently ignored when uploading to NOMAD, although they of course will
remain present in the raw uploaded hdf5 file.

<!-- The H5MD specification can be complemented by modules specific to a
domain of research.  A module may define additional data elements within the
H5MD structure, add conditions that the data must satisfy, or define rules for
their semantic interpretation. Multiple modules may be present, as long as
their prescriptions are not contradictory. Each module is identified by a name
and a version number.

The modules that apply to a specific H5MD structure are stored as subgroups
within the group `h5md/modules`. Each module holds its version number as an
attribute, further module-specific information may be stored:

    h5md
     \-- (modules)
          \-- <module1>
          |    +-- version: Integer[2]
          \-- <module2>
          |    +-- version: Integer[2]
          \-- ...

`version`
:   An attribute, of `Integer` datatype and of simple dataspace of rank 1 and
    size 2, that contains the major version number and the minor version number
    of the module.

    The version *x.y.z* of an H5MD module follows [semantic versioning][semver]
    [@semantic_versioning] and again only the components *x* and *y* are
    stored, see `h5md/version` in "[H5MD metadata]."

[semver]: http://semver.org/spec/v2.0.0.html -->




## The particles group

Particle attributes, i.e., information about each particle in the system, are stored within the `particles` group.
According to the original H5MD schema, the `particles` group is a container for subgroups that
represent different subsets of the system under consideration.
For simplicity of parsing, H5MD-NOMAD currently requires one such group, labeled `all`, to contain all the particles and corresponding attributes to be stored in the NOMAD archive.
**Additional particle groups will be ignored**.

For each dataset, the ordering of indices (whenever relevant) is as follows: frame index, particle index, dimension index.
Thus, the contents of the `particles` group for a trajectory with `N_frames` frames and `N_part` particles in a `D`-dimensional space can be represented:

    particles
     \-- all
     |    \-- box
     |    \-- (<time-dependent_vector_attribute>)
     |    |    \-- step: Integer[N_frames]
     |    |    \-- time: Float[N_frames]
     |    |    \-- value: <type>[N_frames][N_part][D]
     |    \-- (<time-dependent_scalar_attribute>)
     |    |    \-- step: Integer[N_frames]
     |    |    \-- time: Float[N_frames]
     |    |    \-- value: <type>[N_frames][N_part]
     |    \-- (<time-independent_vector_attribute>): <type>[N_part][D]
     |    \-- (<time-independent_scalar_attribute>): <type>[N_part]
     |    \-- ...
     \-- <group2>
          \-- ...

### Standardized H5MD elements for particles group

`position`
:   **(required for parsing other particle attributes)** An element that describes the particle positions as coordinate vectors of `Float` or `Integer` type.

<!-- If the component $k$ of `box/boundary` (see [below](#simulation-box)) is set
to `none`, the data indicate for each particle the component $k$ of its
absolute position in space. If the component $k$ of `box/boundary` is set to
`periodic`, the data indicate for each particle the component $k$ of the
absolute position in space of an *arbitrary* periodic image of that particle. -->

`velocity`
:   An element that contains the velocities for each particle as a vector of
    `Float` or `Integer` type.

`force`
:   An element that contains the total forces (i.e., the accelerations
    multiplied by the particle mass) for each particle as a vector of `Float`
    or `Integer` type.

`mass`
:   An element that holds the mass for each particle as a scalar of `Float`
    type.

- `image`
:   <a id="image_anchor"></a>**(currently unused in H5MD-NOMAD)**

# TODO can we make these admonitions indented somehow or more obviously connected with the members of this list?
??? details

    An element that represents periodic images of the box as coordinate vectors
    of `Float` or `Integer` type and allows one to compute for each particle its
    absolute position in space. If `image` is present, `position` must be
    present as well. For time-dependent data, the `step` and `time` datasets of
    `image` must equal those of `position`, which must be accomplished by
    hard-linking the respective datasets.

<!-- If the component $k$ of `box/boundary` (see [below](#simulation-box)) is set
to `none`, the values of the corresponding component $k$ of `image` serve as
placeholders. If the component $k$ of `box/boundary` is set to `periodic`,
for a cuboid box, the component $k$ of the absolute position of particle $i$
is computed as $R_{ik} = r_{ik} + L_k a_{ik}$, where $\vec r_i$ is taken
from `position`, $\vec a_i$ is taken from `image`, and $\vec L$ from
`box/edges`. -->

`species`
:   **(currently unused in H5MD-NOMAD)**

??? details

    An element that describes the species for each particle, i.e., its
    atomic or chemical identity, as a scalar of `Enumeration` or `Integer`
    data type. Particles of the same species are assumed to be identical with
    respect to their properties and unbonded interactions.

`id`
: **(currently unused in H5MD-NOMAD)**
??? details

    An element that holds a scalar identifier for each particle of `Integer`
    type, which is unique within the given particle subgroup. The `id` serves
    to identify particles over the course of the simulation in the case when
    the order of the particles changes, or when new particles are inserted and
    removed. If `id` is absent, the identity of the particles is given by their
    index in the `value` datasets of the elements within the same subgroup.

<!-- A *fill value* (see
[§ 6.6](http://www.hdfgroup.org/HDF5/doc/UG/11_Datatypes.html#Fvalues) in
[@HDF5_users_guide]) may be defined for `id/value` upon dataset creation.
When the identifier of a particle is equal to this user-defined value,
the particle is considered non-existing, the entry serves as a
placeholder. This permits the storage of subsystems whose number of
particles varies in time. For the case of varying particle number, the
dimension denoted by `[N]` above may be variable. -->

`charge`
:   An element that contains the charge associated to each particle as a
    scalar, of `Integer` or `Float` type.

<!-- `charge` has the optional attribute `type` of fixed-length string datatype
and of scalar dataspace, possible values are `effective` and `formal`. In
the case `effective`, the charge is part of an effective description of the
interactions with the precise meaning depending on the underlying empirical
force fields or coarse-grained models.

In the case `formal`, the charge is the so-called "formal charge" assigned
to an atom (see <http://en.wikipedia.org/wiki/Formal_charge>) and must be
of `Integer` type. This case corresponds to the entries in PDB files (see
definition in the PDBx/mmCIF dictionary
<http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v40.dic/Items/_atom_site.pdbx_formal_charge.html>).

If none of `effective` or `formal` describes the data properly, the
attribute `type` may be omitted. -->

### Standardized H5MD-NOMAD elements for particles group

`species_label`
:   <a id="species_label_anchor"></a> An element that holds a label (fixed-length string datatype) for each particle. This label denotes the fundamental species type of the particle (e.g., the chemical element label for atoms), regardless of its given interactions within the model. **Both** time-independent and time-dependent `species_label` elements are supported.

`model_label`
:   An element that holds a label (fixed-length string datatype) for each particle. This label denotes the type of particle with respect to the given interactions within the model (e.g., force field) **Currently only time-independent species labels are supported.**

### Non-standard elements in particles group

**All non-standard elements within the particles group are currently ignored by the NOMAD H5MD parser.** In principle, one can store additional custom attributes as configuration-specific observables (see [The observables group](#the-observables-group)).


### The simulation box subgroup

Information about the simulation box is stored in a subgroup named `box`, within the relevant particles group (`all` in our case).
**Both** time-independent and time-dependent box information are supported (i.e. via the `edges` element).
Because the `box` group is specific to a particle group of particles, time-dependent boxes **must** contain `step` and `time` datasets that exactly match those of the corresponding `position` group.
In principal, this should be accomplished by hard-linking the respective datasets.
In practice, H5MD-NOMAD currently assumes that this is the case (i.e., **the box group `step` and `time` information is unused**), and simply checks that `edges.value` has the same leading dimension as `position`.

The structure of the `box` group is as follows:

    particles
     \-- all
          \-- box
               +-- dimension: Integer[]
               +-- boundary: String[D]
               \-- (edges)

`dimension`
:   An attribute that stores the spatial dimension `D` of the simulation box
    and is of `Integer` datatype and scalar dataspace.

`boundary`
:   <a id="boundary_anchor"></a> An attribute, of boolean datatype (**changed from string to boolean in H5MD-NOMAD**) and of simple dataspace of rank 1 and size `D`, that specifies the boundary condition of the box along each dimension, i.e., `True` implies periodic boundaries are applied in the corresponding dimension. If all values in `boundary` are `False`, `edges` may be omitted.

<!-- Information on the geometry of the box edges is stored as an H5MD element,
allowing for the box to be fixed in time or not.  Supported box shapes are the
cuboid and triclinic unit cell, for other shapes a transformation to the
triclinic shape may be considered [@Bekker:1997]. -->

`edges`
:   A `D`-dimensional vector or a `D` × `D` matrix, depending on the geometry
of the box, of `Float` or `Integer` type. Only cuboid and triclinic boxes are allowed.
If `edges` is a vector, it specifies the space diagonal of a cuboid-shaped box. If `edges` is a
matrix, the box is of triclinic shape with the edge vectors given by the
rows of the matrix. For a time-dependent box, a cuboid geometry is encoded by a dataset `value`
(within the H5MD element) of rank 2 (1 dimension for the time and 1 for the
vector) and a triclinic geometry by a dataset `value` of rank 3 (1
dimension for the time and 2 for the matrix). For a time-independent box, a cuboid geometry is encoded by a dataset `edges` of rank 1 and a triclinic geometry by a dataset of rank 2.

For instance, a cuboid box that changes in time would appear as:

    particles
     \-- all
          \-- box
               +-- dimension: Integer[]
               +-- boundary: String[D]
               \-- edges
                    \-- step: Integer[variable]
                    \-- time: Float[variable]
                    \-- value: <type>[variable][D]

where `dimension` is equal to `D`. A triclinic box that is fixed in time would
appear as:

    particles
     \-- all
          \-- box
               +-- dimension: Integer[]
               +-- boundary: String[D]
               \-- edges: <type>[D][D]

where `dimension` is equal to `D`.

<!-- TODO - double check that both shapes are supported in the parser! -->




## The observables group

The initial H5MD proposed a simple and flexible schema for the general storage of observable info, defined roughly as "macroscopic observables" or "averages of a property over many particles", as H5MD elements:

    observables
     \-- <observable1>
     |    \-- step: Integer[N_frames]
     |    \-- time: Float[N_frames]
     |    \-- value: <type>[N_frames]
     \-- <observable2>
     |    \-- step: Integer[N_frames]
     |    \-- time: Float[N_frames]
     |    \-- value: <type>[N_frames][D]
     \-- <group1>
     |    \-- <observable3>
     |         \-- step: Integer[N_frames]
     |         \-- time: Float[N_frames]
     |         \-- value: <type>[N_frames][D][D]
     \-- <observable4>: <type>[]
     \-- ...

<a id="obs_para2"></a>

As depicted above, observables representing only a subset of the particles may be stored in appropriate subgroups similar to the `particles` tree. H5MD-NOMAD **does** support the organization of observables into subgroups (as discussed in more detail below). However, **grouping by particle groups is not fully supported** in the sense that there is currently no metadata storing the corresponding indices of the relevant particles subgroup. Additionally, since [only the `all` particles group is parsed](#the-particles-group), information about the named subgroup will not be stored anywhere in the archive. *Thus, we recommend for now that only observables relevant to the `all` particles subgroup are stored within this section.*
<!-- TODO - not sure about this, it might be fine if you can add additional metadata that is stored -->

### H5MD-NOMAD observables

H5MD-NOMAD extends H5MD observable storage by 1. specifying standard observable types with associated metadata and 2. providing standardized specifications for some common observables. In contrast to the schema above, a more restrictive structure is required:

    observables
     \-- <observable_type_1>
     |    \-- <observable_1_label_1>
     |    |    +-- type: String[]
     |    |    \-- ...
     \-- <observable_type_2>
     |    \-- <observable_2_label_1>
     |    |    +-- type: String[]
     |    |    \-- ...
     |    \-- <observable_2_label_2>
     |    |    +-- type: String[]
     |    |    \-- ...
     |    \-- ...
     \-- ...

Here, each `observable_type` corresponds to a particular group of observables, e.g., to be plotted together in a single plot. The given name for this group could be generic, e.g., `radial distribution function`, or more specific, e.g., `molecular radial distribution function for solvents`. The latter may be useful in case multiple groupings of a single type of observable are needed.
Each `observable_label` then corresponds to a specific name for an individual instance of this observable type. For example, for a radial distribution function between particles of type `A` and `B`, `observable_label` might be set to `A-B`.

Finally, H5MD-NOMAD has added the observable `type` as an attribute of each observable:
The following observable types are supported:

<a id="configurational_observable_anchor"></a>

`configurational`
:   An observable that is computed for each individual configuration, with the following general structure:

    observables
     \-- <configurational_subgroup>
     |    \-- <label_1>
     |    |    +-- type: "configurational"
     |    |    \-- step: Integer[N_frames]
     |    |    \-- time: Float[N_frames]
     |    |    \-- value: <type>[N_frames][M]
     |    \-- ...
     \-- ...
 where `M` is the dimension of the observable. This section may also be used to store per-particle quantities/attributes that are not currently supported as [standardized H5MD-NOMAD elements for particles group](#standardized-h5md-elements-for-particles-group), in which case `value` will have dimensions `[N_frames][N_part][M]`.

<a id="ensemble_average_observable_anchor"></a>

`ensemble_average`
:   An observable that is computed by averaging over multiple configurations, with the following generic structure:

    observables
     \-- <ensemble_average_subgroup>
     |    \-- <label_1>
     |    |    +-- type: "ensemble_average"
     |    |    \-- (n_variables): Integer
     |    |    \-- (variables_name): String[n_variables][]
     |    |    \-- (n_bins): Integer[]
     |    |    \-- bins: Float[n_bins][]
     |    |    \-- value: <type>[n_bins][]
     |    |    \-- (frame_start): Integer
     |    |    \-- (frame_end): Integer
     |    |    \-- (n_smooth): Integer
     |    |    \-- (type): String[]
     |    |    \-- (error_type): String[]
     |    |    \-- (errors): Float[n_bins]
     |    |    \-- (error_labels): String[]
     |    |    \-- (frame_end): Integer
     |    |    \-- (<custom_dataset>): <type>[]
     |    \-- ...
     \-- ...

* `n_variables`
:   dimensionality of the observable. Can also be inferred from leading dimension of `bins`.

* `variables_name`
:   name/description of the independent variables along which the observable is defined.

* `n_bins`
:   number of bins along each dimension of the observable. Either single Integer for 1-D observables, or a list of Integers for multi-dimensional observable. Can also be inferred from dimensions of `bins`.

* `bins`
:   value of the bins used for calculating the observable along each dimension of the observable.

* `value`
:   value of the calculated ensemble average at each bin.

* `frame_start`
:   trajectory frame index at which the averaging begins. **This index must correspond to the list of steps and times in `particles.all.position`.**

* `frame_end`
:   trajectory frame index at which the averaging ends. **This index must correspond to the list of steps and times in `particles.all.position`.**

* `n_smooth`
:   number of bins over which the running average was computed for `value`.

* `type`
:   Allowed values of `molecular` or `atomic`. Categorizes if the observable is calculated at the molecular or atomic level.
<!-- TODO - not sure if this is useful -->

* `error_type`
:   describes the type of error reported for this observable. Examples: `Pearson correlation coefficient`, `mean squared error`.

* `errors`
:   value of the error at each bin. Can be multidimensional with corresponding label stored in `error_labels`.

* `error_labels`
:   describes the error along individual dimensions for multi-D errors.

* `<custom_dataset>`
:   additional metadata may be given as necessary.
<!-- TODO - Is this really parsed?! -->

<a id="time_correlation_observable_anchor"></a>

`time_correlation`
:   An observable that is computed by calculating correlations between configurations in time, with the following general structure:

    observables
     \-- <time_correlation_subgroup>
     |    \-- <label_1>
     |    |    +-- type: "time_correlation"
     |    |    \-- (direction): String[]
     |    |    \-- (n_times): Integer[]
     |    |    \-- times: Float[n_times][]
     |    |    \-- value: <type>[n_bins][]
     |    |    \-- (type): String[]
     |    |    \-- (error_type): String[]
     |    |    \-- (errors): Float[n_bins]
     |    |    \-- (error_labels): String[]
     |    |    \-- (<custom_dataset>): <type>[]
     |    \-- ...
     \-- ...

* `label`
:   describes the particles involved in determining the property. For example, for a radial distribution function between particles of type `A` and `B`, `label` might be set to `A-B`

* `direction`
:   allowed values of `x`, `y`, `z`, `xy`, `yz`, `xz`, `xyz`. The direction/s used for calculating the correlation function.

* `n_times`
:   number of times windows for the calculation of the correlation function. Can also be inferred from dimensions of `times`.

* `times`
:   time values used for calculating the correlation function (i.e., &Delta;t values).

* `value`
:   value of the calculated correlation function at each time.

* `type`
:   Allowed values of `molecular` or `atomic`. Categorizes if the observable is calculated at the molecular or atomic level.
<!-- TODO - not sure if this is useful -->

* `error_type`
:   describes the type of error reported for this observable. Examples: `Pearson correlation coefficient`, `mean squared error`.

* `errors`
:   value of the error at each bin. Can be multidimensional with corresponding label stored in `error_labels`.

* `error_labels`
:   describes the error along individual dimensions for multi-D errors.

* `<custom_dataset>`
:   additional metadata may be given as necessary.
<!-- TODO - Is this really parsed?! -->

A list of standardized observables can be found in [Reference - H5MD-NOMAD > Standardized observables in H5MD-NOMAD](h5md_ref.md#standardized-observables-in-h5md-nomad).




## The connectivity group

The initial H5MD proposed a simple and flexible schema for the storage of "connectivity" information, e.g., to be used in conjunction with a molecular mechanics force field.
The connectivity information is stored as tuples in the group
`/connectivity`. The tuples are pairs, triples, etc. as needed and may be either
time-independent or time-dependent.
As with other elements, connectivity elements can be defined for particular particle groups. However, H5MD-NOMAD focuses on the storage of connectivity elements for the entire system (i.e., the `all` particles group).

### Standardized H5MD-NOMAD connectivity

The general structure of the `connectivity` group is as follows:

    connectivity
     \-- (bonds): Integer[N_part][2]
     \-- (angles): Integer[N_part][3]
     \-- (dihedrals): Integer[N_part][4]
     \-- (impropers): Integer[N_part][4]
     \-- (<custom_interaction>): Integer[N_part][m]
     \-- (particles_group)
          \-- ...

`N_part` corresponds to the number of particles stored in the `particles/all` group.

* `bonds`
: a list of 2-tuples specifying the indices of particles containing a "bond interaction".

* `angles`
: a list of 3-tuples specifying the indices of particles containing an "angle interaction".

* `dihedrals`
: a list of 4-tuples specifying the indices of particles containing a "dihedral interaction".

* `impropers`
: a list of 4-tuples specifying the indices of particles containing an "improper dihedral interaction".

* `<custom_interaction>`
: a list of m-tuples specifying the indices of particles containing an arbitrary interaction. `m` denotes the number of particles involved in the interaction.

* `particles_group`
: See below.

<a id="connectivity_support_anchor"></a>
**Currently only time-independent connectivity elements are supported.**

### The particles_group subgroup

Despite not fully utilizing the organization of arbitrary groups of particles within the `particles` group, H5MD-NOMAD allows for the user to provide an arbitrary hierarchy of particle groupings, also referred to as a "topology", within the `connectivity` subgroup called `particles_group`. This information will be used by NOMAD to facilitate visualizations of the system, through the "topology bar" in the overview page. The general structure of the topology group is as follows:

    connectivity
     \-- particles_group
          \-- <group_1>
          |    \-- (type): String[]
          |    \-- (formula): String[]
          |    \-- indices: Integer[]
          |    \-- (is_molecule): Bool
     |    |    \-- (<custom_dataset>): <type>[]
          |    \-- (particles_group):
          |        \-- ...
          \-- <group_2>
              \-- ...

The initial `particles_group` subgroup, directly under `connectivity`, is a container for the entire topology. `particles_group` contains a series of subgroups with arbitrary names, which denote the first level of organization within the topology. The name of each subgroup will become the group label within the NOMAD metadata. Each of these subgroups then contain a series of datasets:

* `type`
: describes the type of particle group. There exists a list of standardized types: `molecule_group`, `molecule`, `monomer_group`, `monomer`. However, arbitrary types can be given. We suggest that you 1. use the standardized types when appropriate (note that protein residues should be generically typed as `monomer`) and 2. use the general format `<type>_group` for groups of a distinct type (see further description of suggested hierarchy below).

* `formula`
: a "chemical-like" formula that describes the particle group with respect to its underlying components. The format for the formula is `<child_1>(n_child_1)<child_2>(n_child_2)...`, where `<child_x>` is the name/label of the underlying component, and `n_child_x` is the number of such components found within this particle group. Example: A particles group containing 100 water molecules named `water` has the formula `water(100)`, whereas each underlying water molecule has the standard chemical formula `H2O`.

* `indices`
: a list of integer indices corresponding to all particles belonging to this group. Indices should correspond to the list of particles stored in the `particles/all` group.

* `is_molecule`
: indicator of individual molecules (typically with respect to the bond connections defined by a force field).

* `custom_dataset`
: arbitrary additional metadata for this particle group may be given.


Each subgroup may also contain a (nested) `particles_group` subgroup, in order to subdivide the group of particles into an organizational hierarchy. As with the overall `particles_group` container, the groups contained within `particles_group` must not *partition* the particles within this group (i.e., overlapping or non-complete groupings are allowed). However, particle groups *must* contain particles already contained within the parent `particles_group` (i.e., subgroups must be a subset of the grouping at the previous level of the hierarchy).

Note that typically the `particles_group` hierarchy ends at the level of individual particles (i.e., individual particles are not stored, since this information is already contained within the `particles` group).



## The parameters group

The initial H5MD proposed a simple and flexible schema for the storage of general "parameter" information within the `parameters` group, with the following structure:

    parameters
     +-- <user_attribute1>
     \-- <user_data1>
     \-- <user_group1>
     |    \-- <user_data2>
     |    \-- ...
     \-- ...

In contrast, the H5MD-NOMAD schema calls for very specific structures to be used when storing parameter information. While the previous groups have attempted to stay away from enforcing NOMAD-specific data structures on the user, instead opting for more intuitive and generally-convenient structures, the `parameters` group utilizes already-existing metadata and structures within NOMAD to efficiently import simulation parameters in a way that is searchable and comparable to simulations performed by other users.

In this way, the H5MD-NOMAD `parameters` group has the following structure:

    parameters
     \-- <parameter_subgroup_1>
     |    \-- ...
     \-- <parameter_subgroup_2>
     |    \-- ...
     \-- ...

The subgroups `force_calculations` and `workflow` are supported. The following describes the detailed data structures for these subgroups, using the NOMAD MetaInfo definitions for each underlying `Quantity`. Please note that:
<!-- TODO add href to glossary -->

1. Quantities with `type=MEnum()` are restricted to the provided allowed values.

2. The unit given in the MetaInfo definition does not have to be used within the H5MD-NOMAD file, however, the dimensionality of the unit should match.


### Force calculations

The `force_calculations` group contains the parameters for force calculations according to the force field during a molecular dynamics run.

<a id="force_calculation_template_anchor"></a>

The following json template illustrates the structure of the `force_calculations` group, with example values for clarity:

```json
{
    "vdw_cutoff": {"value": 1.2, "unit": "nm"},
    "coulomb_type": "particle_mesh_ewald",
    "coulomb_cutoff": {"value": 1.2, "unit": "nm"},
    "neighbor_searching": {
        "neighbor_update_frequency": 1,
        "neighbor_update_cutoff": {"value": 1.2, "unit": "nm"}
        }
    }
```

In the following, we provide the NOMAD definitions for each of these quantities:

* `vdw_cutoff`
:

        Quantity(
                type=np.float64,
                shape=[],
                unit='m',
                description='''
                Cutoff for calculating VDW forces.
                ''')

* `coulomb_type`
:

        Quantity(
            type=MEnum('cutoff', 'ewald', 'multilevel_summation', 'particle_mesh_ewald',
                    'particle_particle_particle_mesh', 'reaction_field'),
            shape=[],
            description='''
            Method used for calculating long-ranged Coulomb forces.

            Allowed values are:

            | Barostat Name          | Description                               |

            | ---------------------- | ----------------------------------------- |

            | `""`                   | No thermostat               |

            | `"Cutoff"`          | Simple cutoff scheme. |

            | `"Ewald"` | Standard Ewald summation as described in any solid-state physics text. |

            | `"Multi-Level Summation"` |  D. Hardy, J.E. Stone, and K. Schulten,
            [Parallel. Comput. **35**, 164](https://doi.org/10.1016/j.parco.2008.12.005)|

            | `"Particle-Mesh-Ewald"`        | T. Darden, D. York, and L. Pedersen,
            [J. Chem. Phys. **98**, 10089 (1993)](https://doi.org/10.1063/1.464397) |

            | `"Particle-Particle Particle-Mesh"` | See e.g. Hockney and Eastwood, Computer Simulation Using Particles,
            Adam Hilger, NY (1989). |

            | `"Reaction-Field"` | J.A. Barker and R.O. Watts,
            [Mol. Phys. **26**, 789 (1973)](https://doi.org/10.1080/00268977300102101)|
            ''')


* `coulomb_cutoff`
:

        Quantity(
            type=np.float64,
            shape=[],
            unit='m',
            description='''
            Cutoff for calculating short-ranged Coulomb forces.
            ''')

* `neighbor_searching`
:
Section containing the parameters for neighbor searching/lists during a molecular dynamics run.

* `neighbor_update_frequency`
:

        Quantity(
            type=int,
            shape=[],
            description='''
            Number of timesteps between updating the neighbor list.
            ''')

* `neighbor_update_cutoff`
:

        Quantity(
            type=np.float64,
            shape=[],
            unit='m',
            description='''
            The distance cutoff for determining the neighbor list.
            ''')


### The molecular dynamics workflow

The `workflow` group contains the parameters for any type of workflow. Here we describe the specific case of the well-defined `molecular_dynamics` workflow. Custom workflows are described in detail in [Workflows in NOMAD](../overview.md).

<a id="md_workflow_template_anchor"></a>

The following json template illustrates the structure of the `molecular_dynamics` subsection of the `workflow` group, with example values for clarity:

```json
{
    "molecular_dynamics": {
        "thermodynamic_ensemble": "NPT",
        "integrator_type": "langevin_leap_frog",
        "integration_timestep": {"value": 2e-15, "unit": "ps"},
        "n_steps": 20000000,
        "coordinate_save_frequency": 10000,
        "velocity_save_frequency": null,
        "force_save_frequency": null,
        "thermodynamics_save_frequency": null,
        "thermostat_parameters": {
            "thermostat_type": "langevin_leap_frog",
            "reference_temperature": {"value": 300.0, "unit": "kelvin"},
            "coupling_constant": {"value": 1.0, "unit": "ps"}},
        "barostat_parameters": {
            "barostat_type": "berendsen",
            "coupling_type": "isotropic",
            "reference_pressure": {"value": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], "unit": "bar"},
            "coupling_constant": {"value": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]},
            "compressibility": {"value": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]}
            }
    }
}
```

In the following, we provide the NOMAD definitions for each of these quantities:

* `thermodynamic_ensemble`
:

        Quantity(
            type=MEnum('NVE', 'NVT', 'NPT', 'NPH'),
            shape=[],
            description='''
            The type of thermodynamic ensemble that was simulated.

            Allowed values are:

            | Thermodynamic Ensemble          | Description                               |

            | ---------------------- | ----------------------------------------- |

            | `"NVE"`           | Constant number of particles, volume, and energy |

            | `"NVT"`           | Constant number of particles, volume, and temperature |

            | `"NPT"`           | Constant number of particles, pressure, and temperature |

            | `"NPH"`           | Constant number of particles, pressure, and enthalpy |
            ''')

* `integrator_type`
:
        Quantity(
            type=MEnum(
                'brownian', 'conjugant_gradient', 'langevin_goga',
                'langevin_schneider', 'leap_frog', 'rRESPA_multitimescale', 'velocity_verlet'
            ),
            shape=[],
            description='''
            Name of the integrator.

            Allowed values are:

            | Integrator Name          | Description                               |

            | ---------------------- | ----------------------------------------- |

            | `"langevin_goga"`           | N. Goga, A. J. Rzepiela, A. H. de Vries,
            S. J. Marrink, and H. J. C. Berendsen, [J. Chem. Theory Comput. **8**, 3637 (2012)]
            (https://doi.org/10.1021/ct3000876) |

            | `"langevin_schneider"`           | T. Schneider and E. Stoll,
            [Phys. Rev. B **17**, 1302](https://doi.org/10.1103/PhysRevB.17.1302) |

            | `"leap_frog"`          | R.W. Hockney, S.P. Goel, and J. Eastwood,
            [J. Comp. Phys. **14**, 148 (1974)](https://doi.org/10.1016/0021-9991(74)90010-2) |

            | `"velocity_verlet"` | W.C. Swope, H.C. Andersen, P.H. Berens, and K.R. Wilson,
            [J. Chem. Phys. **76**, 637 (1982)](https://doi.org/10.1063/1.442716) |

            | `"rRESPA_multitimescale"` | M. Tuckerman, B. J. Berne, and G. J. Martyna
            [J. Chem. Phys. **97**, 1990 (1992)](https://doi.org/10.1063/1.463137) |
            ''')

* `integration_timestep`
:
        Quantity(
            type=np.float64,
            shape=[],
            unit='s',
            description='''
            The timestep at which the numerical integration is performed.
            ''')

* `n_steps`
:
        Quantity(
            type=int,
            shape=[],
            description='''
            Number of timesteps performed.
            ''')

* `coordinate_save_frequency`
:
        Quantity(
            type=int,
            shape=[],
            description='''
            The number of timesteps between saving the coordinates.
            ''')

* `velocity_save_frequency`
:
        Quantity(
            type=int,
            shape=[],
            description='''
            The number of timesteps between saving the velocities.
            ''')

* `force_save_frequency`
:
        Quantity(
            type=int,
            shape=[],
            description='''
            The number of timesteps between saving the forces.
            ''')

* `thermodynamics_save_frequency`
:
        Quantity(
            type=int,
            shape=[],
            description='''
            The number of timesteps between saving the thermodynamic quantities.
            ''')

* `thermostat_parameters`
:  Section containing the parameters pertaining to the thermostat for a molecular dynamics run.

* `thermostat_type`
:

        Quantity(
            type=MEnum('andersen', 'berendsen', 'brownian', 'langevin_goga', 'langevin_schneider', 'nose_hoover', 'velocity_rescaling',
                    'velocity_rescaling_langevin'),
            shape=[],
            description='''
            The name of the thermostat used for temperature control. If skipped or an empty string is used, it
            means no thermostat was applied.

            Allowed values are:

            | Thermostat Name        | Description                               |

            | ---------------------- | ----------------------------------------- |

            | `""`                   | No thermostat               |

            | `"andersen"`           | H.C. Andersen, [J. Chem. Phys.
            **72**, 2384 (1980)](https://doi.org/10.1063/1.439486) |

            | `"berendsen"`          | H. J. C. Berendsen, J. P. M. Postma,
            W. F. van Gunsteren, A. DiNola, and J. R. Haak, [J. Chem. Phys.
            **81**, 3684 (1984)](https://doi.org/10.1063/1.448118) |

            | `"brownian"`           | Brownian Dynamics |

            | `"langevin_goga"`           | N. Goga, A. J. Rzepiela, A. H. de Vries,
            S. J. Marrink, and H. J. C. Berendsen, [J. Chem. Theory Comput. **8**, 3637 (2012)]
            (https://doi.org/10.1021/ct3000876) |

            | `"langevin_schneider"`           | T. Schneider and E. Stoll,
            [Phys. Rev. B **17**, 1302](https://doi.org/10.1103/PhysRevB.17.1302) |

            | `"nose_hoover"`        | S. Nosé, [Mol. Phys. **52**, 255 (1984)]
            (https://doi.org/10.1080/00268978400101201); W.G. Hoover, [Phys. Rev. A
            **31**, 1695 (1985) |

            | `"velocity_rescaling"` | G. Bussi, D. Donadio, and M. Parrinello,
            [J. Chem. Phys. **126**, 014101 (2007)](https://doi.org/10.1063/1.2408420) |

            | `"velocity_rescaling_langevin"` | G. Bussi and M. Parrinello,
            [Phys. Rev. E **75**, 056707 (2007)](https://doi.org/10.1103/PhysRevE.75.056707) |
            ''')

* `reference_temperature`
:

        Quantity(
            type=np.float64,
            shape=[],
            unit='kelvin',
            description='''
            The target temperature for the simulation.
            ''')

* `coupling_constant`
:

        Quantity(
            type=np.float64,
            shape=[],
            unit='s',
            description='''
            The time constant for temperature coupling. Need to describe what this means for the various
            thermostat options...
            ''')

* `effective_mass`
:

        Quantity(
            type=np.float64,
            shape=[],
            unit='kilogram',
            description='''
            The effective or fictitious mass of the temperature resevoir.
            ''')

* `barostat_parameters`
: Section containing the parameters pertaining to the barostat for a molecular dynamics run.

* `barostat_type`
:

        Quantity(
            type=MEnum('berendsen', 'martyna_tuckerman_tobias_klein', 'nose_hoover', 'parrinello_rahman', 'stochastic_cell_rescaling'),
            shape=[],
            description='''
            The name of the barostat used for temperature control. If skipped or an empty string is used, it
            means no barostat was applied.

            Allowed values are:

            | Barostat Name          | Description                               |

            | ---------------------- | ----------------------------------------- |

            | `""`                   | No thermostat               |

            | `"berendsen"`          | H. J. C. Berendsen, J. P. M. Postma,
            W. F. van Gunsteren, A. DiNola, and J. R. Haak, [J. Chem. Phys.
            **81**, 3684 (1984)](https://doi.org/10.1063/1.448118) |

            | `"martyna_tuckerman_tobias_klein"` | G.J. Martyna, M.E. Tuckerman, D.J. Tobias, and M.L. Klein,
            [Mol. Phys. **87**, 1117 (1996)](https://doi.org/10.1080/00268979600100761);
            M.E. Tuckerman, J. Alejandre, R. López-Rendón, A.L. Jochim, and G.J. Martyna,
            [J. Phys. A. **59**, 5629 (2006)](https://doi.org/10.1088/0305-4470/39/19/S18)|

            | `"nose_hoover"`        | S. Nosé, [Mol. Phys. **52**, 255 (1984)]
            (https://doi.org/10.1080/00268978400101201); W.G. Hoover, [Phys. Rev. A
            **31**, 1695 (1985) |

            | `"parrinello_rahman"`        | M. Parrinello and A. Rahman,
            [J. Appl. Phys. **52**, 7182 (1981)](https://doi.org/10.1063/1.328693);
            S. Nosé and M.L. Klein, [Mol. Phys. **50**, 1055 (1983) |

            | `"stochastic_cell_rescaling"` | M. Bernetti and G. Bussi,
            [J. Chem. Phys. **153**, 114107 (2020)](https://doi.org/10.1063/1.2408420) |
            ''')

* `coupling_type`
:

        Quantity(
            type=MEnum('isotropic', 'semi_isotropic', 'anisotropic'),
            shape=[],
            description='''
            Describes the symmetry of pressure coupling. Specifics can be inferred from the `coupling constant`

            | Type          | Description                               |

            | ---------------------- | ----------------------------------------- |

            | `isotropic`          | Identical coupling in all directions. |

            | `semi_isotropic` | Identical coupling in 2 directions. |

            | `anisotropic`        | General case. |
            ''')

* `reference_pressure`
:

        Quantity(
            type=np.float64,
            shape=[3, 3],
            unit='pascal',
            description='''
            The target pressure for the simulation, stored in a 3x3 matrix, indicating the values for individual directions
            along the diagonal, and coupling between directions on the off-diagonal.
            ''')

* `coupling_constant`
:

        Quantity(
            type=np.float64,
            shape=[3, 3],
            unit='s',
            description='''
            The time constants for pressure coupling, stored in a 3x3 matrix, indicating the values for individual directions
            along the diagonal, and coupling between directions on the off-diagonal. 0 values along the off-diagonal
            indicate no-coupling between these directions.
            ''')

* `compressibility`
:

        Quantity(
            type=np.float64,
            shape=[3, 3],
            unit='1 / pascal',
            description='''
            An estimate of the system's compressibility, used for box rescaling, stored in a 3x3 matrix indicating the values for individual directions
            along the diagonal, and coupling between directions on the off-diagonal. If None, it may indicate that these values
            are incorporated into the coupling_constant, or simply that the software used uses a fixed value that is not available in
            the input/output files.
            ''')



## Units

In the original H5MD schema, units were given as string attributes of datasets, e.g., ``60 m s-2``.
H5MD-NOMAD amends the treatment of units in 2 ways:

1. If needed, the leading prefactor is stored as a separate attribute of `float` datatype called `unit_factor`.

2. The string that describes the unit should be compatible with the `UnitRegistry` class of the `pint` python module.

Generic representation of unit storage in H5MD-NOMAD:

    <group>
        \-- <dataset>
            +-- (unit: String[])
            +-- (unit_factor: Float)