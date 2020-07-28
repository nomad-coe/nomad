from abc import ABC, abstractmethod
from collections import OrderedDict
import numpy as np
from typing import Tuple, List
from nomad.units import ureg

from nomad.metainfo import Section
from nomad.utils import RestrictedDict


def get_basis_set(context, entry_archive, logger) -> RestrictedDict:
    """Decide which type of basis set settings are applicable to the entry and
    return a corresponding settings as a RestrictedDict.

    Args:
        context: The calculation context.
        entry_archive: EntryArchive from which values are extracted.
        logger: Shared logger.

    Returns:
        RestrictedDict or None: Returns the extracted settings as a
        RestrictedDict. If no suitable basis set settings could be identified,
        returns None.
    """
    settings: BasisSet = None
    program_name = entry_archive.section_run[0].program_name
    if program_name == "exciting":
        settings = BasisSetExciting(context, entry_archive, logger)
    elif program_name == "FHI-aims":
        settings = BasisSetFHIAims(context, entry_archive, logger)
    else:
        return None

    return settings.to_dict()


class BasisSet(ABC):
    """Abstract base class for basis set settings. The idea is to create
    subclasses that inherit this class and hierarchically add new mandatory and
    optional settings with the setup()-function.
    """
    def __init__(self, context, entry_archive, logger):
        """
        """
        self._ctx = context
        self._entry_archive = entry_archive
        self._logger = logger
        mandatory, optional = self.setup()
        self.settings = RestrictedDict(mandatory, optional, forbidden_values=[None])

    @abstractmethod
    def to_dict(self) -> RestrictedDict:
        """Used to extract basis set settings from the backend and returning
        them as a RestrictedDict.
        """
        pass

    @abstractmethod
    def setup(self) -> Tuple:
        """Used to define a list of mandatory and optional settings for a
        subclass.

        Returns:
            Should return a tuple of two lists: the first one defining
            mandatory keys and the second one defining optional keys.
        """
        mandatory: List = []
        optional: List = []
        return mandatory, optional


class BasisSetFHIAims(BasisSet):
    """Basis set settings for 'FHI-Aims' (code-dependent).
    """
    def setup(self) -> Tuple:
        # Get previously defined values from superclass
        mandatory, optional = super().setup()

        # Add new values
        mandatory += ["fhiaims_basis"]

        return mandatory, optional

    def to_dict(self):
        # Get basis set settings for each species
        aims_bs = self._ctx.representative_method.x_fhi_aims_section_controlIn_basis_set
        if aims_bs is not None:
            bs_by_species = {}
            for this_aims_bs in aims_bs:
                this_bs_dict = self._values_to_dict(this_aims_bs, level=2)
                this_species = this_aims_bs['x_fhi_aims_controlIn_species_name'][0]
                bs_by_species[this_species] = this_bs_dict

            # Sort alphabetically by species label
            if bs_by_species:
                basis = OrderedDict()
                for k in sorted(bs_by_species.keys()):
                    basis[k] = bs_by_species[k]
                self.settings["fhiaims_basis"] = basis

        return self.settings

    @classmethod
    def _values_to_dict(cls, data, level=0):
        result = None
        if data is None:
            return None
        elif isinstance(data, (Section, dict)):
            result = OrderedDict()
            for k in sorted(cls._filtered_section_keys(data)):
                v = data.get(k, None)
                result[k] = cls._values_to_dict(v, level=level + 1)
        elif isinstance(data, (list)):
            result = []
            for k in range(len(data)):
                v = data[k]
                result.append(cls._values_to_dict(v, level=level + 1))
        elif isinstance(data, (np.ndarray)):
            result = data.tolist()
        else:
            result = data
        return result

    @classmethod
    def _filtered_section_keys(cls, section):
        for k in section.keys():
            # skip JSON-specific keys
            if k == '_gIndex':
                continue
            if k == '_name':
                continue
            else:
                # json values and subsections
                yield k


class BasisSetExciting(BasisSet):
    """Basis set settings for 'Exciting' (code-dependent).
    """
    def setup(self) -> Tuple:
        # Get previously defined values from superclass
        mandatory, optional = super().setup()

        # Add new values
        mandatory += [
            "muffin_tin_settings",
            "rgkmax",
            "gkmax",
            "lo",
            "lmaxapw",
        ]

        return mandatory, optional

    def to_dict(self):
        """Special case of basis set settings for Exciting code. See list at:
        https://gitlab.mpcdf.mpg.de/nomad-lab/encyclopedia-general/wikis/FHI-visit-preparation
        """
        # Add the muffin-tin settings for each species ordered alphabetically by atom label
        try:
            groups = self._ctx.representative_system.x_exciting_section_atoms_group
            groups = sorted(groups, key=lambda group: group.x_exciting_geometry_atom_labels)
            muffin_tin_settings = OrderedDict()
            for group in groups:
                label = group.x_exciting_geometry_atom_labels
                try:
                    muffin_tin_settings["{}_muffin_tin_radius".format(label)] = "%.6f" % (group.x_exciting_muffin_tin_radius.to(ureg.angstrom).magnitude)
                except Exception:
                    muffin_tin_settings["{}_muffin_tin_radius".format(label)] = None
                try:
                    muffin_tin_settings["{}_muffin_tin_points".format(label)] = "%d" % group.x_exciting_muffin_tin_points
                except Exception:
                    muffin_tin_settings["{}_muffin_tin_points".format(label)] = None
            self.settings["muffin_tin_settings"] = muffin_tin_settings
        except Exception:
            pass

        # Other important method settings
        system = self._ctx.representative_system
        try:
            self.settings['rgkmax'] = "%.6f" % (system.x_exciting_rgkmax.magnitude)
        except Exception:
            pass
        try:
            self.settings['gkmax'] = "%.6f" % (1e-10 * system.x_exciting_gkmax.magnitude)
        except Exception:
            pass
        try:
            self.settings['lo'] = "%d" % (system.x_exciting_lo)
        except Exception:
            pass
        try:
            self.settings['lmaxapw'] = "%d" % (system.x_exciting_lmaxapw)
        except Exception:
            pass

        return self.settings
