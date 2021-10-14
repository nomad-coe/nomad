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

from typing import Any

from nomad.normalizing.normalizer import Normalizer
from nomad.datamodel import EntryMetadata
from nomad.datamodel.encyclopedia import (
    EncyclopediaMetadata,
    Material,
    Method,
    Properties,
    Calculation,
)
from nomad.normalizing.encyclopedia.context import Context
from nomad.normalizing.encyclopedia.material import MaterialBulkNormalizer, Material2DNormalizer, Material1DNormalizer
from nomad.normalizing.encyclopedia.method import MethodDFTNormalizer, MethodGWNormalizer
from nomad.normalizing.encyclopedia.properties import PropertiesNormalizer
from nomad import config


class EncyclopediaNormalizer(Normalizer):
    """
    This normalizer emulates the functionality of the old Encyclopedia backend.
    The data used by the encyclopedia have been assigned under new metainfo
    within a new section called "encyclopedia". In the future these separate
    metainfos could be absorbed into the existing metainfo hiearchy.
    """
    def calc_type(self, calc: Calculation) -> str:
        """Decides what type of calculation this is: single_point, md,
        geometry_optimization, etc.
        """
        calc_enums = Calculation.calculation_type.type
        calc_type = calc_enums.unavailable

        # Primarily try to determine the calculation type from workflow
        # information
        try:
            workflow = self.entry_archive.workflow[-1]
            workflow_map = {
                "molecular_dynamics": calc_enums.molecular_dynamics,
                "geometry_optimization": calc_enums.geometry_optimization,
                "phonon": calc_enums.phonon_calculation,
            }
            workflow_enum = workflow_map.get(workflow.type)
            if workflow_enum is not None:
                calc.calculation_type = workflow_enum
                return workflow_enum
        except Exception:
            pass

        # Fall back to old frame sequence data
        try:
            sccs = self.section_run.calculation
        except Exception:
            sccs = []

        n_scc = len(sccs)

        # No sequences, only a few calculations
        n_workflow = len(self.entry_archive.workflow)
        if n_scc <= 3 and n_workflow == 0:
            calc_type = calc_enums.single_point

        # One sequence. Currently calculations with multiple sequences are
        # unsupported.
        elif n_workflow == 1:
            workflow_type = self.entry_archive.workflow[0].type
            if workflow_type == "molecular_dynamics":
                calc_type = calc_enums.molecular_dynamics
            if workflow_type == "geometry_optimization":
                calc_type = calc_enums.geometry_optimization
            if workflow_type == "phonon":
                calc_type = calc_enums.phonon_calculation
            if workflow_type == 'elastic':
                calc_type = calc_enums.elastic_constants
            if workflow_type == 'single_point':
                calc_type = calc_enums.single_point
        calc.calculation_type = calc_type
        return calc_type

    def material_type(self, material: Material) -> tuple:
        # Try to fetch representative system
        system = None
        material_type = config.services.unavailable_value
        material_enums = Material.material_type.type
        try:
            system_idx = self.section_run.m_cache["representative_system_idx"]
        except (AttributeError, KeyError):
            pass
        else:
            # Try to find system type information from archive for the selected system.
            try:
                system = self.section_run.system[system_idx]
                stype = system.type
            except KeyError:
                pass
            else:
                if stype == material_enums.one_d or stype == material_enums.two_d:
                    material_type = stype
                # For bulk systems we also ensure that the symmetry information is available
                if stype == material_enums.bulk:
                    try:
                        system.symmetry[0]
                    except (KeyError, IndexError):
                        self.logger.info("Symmetry information is not available for a bulk system. No Encylopedia entry created.")
                    else:
                        material_type = stype

        material.material_type = material_type
        return system, material_type

    def method_type(self, method: Method) -> tuple:
        repr_method = None
        method_id = config.services.unavailable_value
        methods = self.section_run.method
        n_methods = len(methods)

        if n_methods == 1:
            repr_method = methods[0]
            method_id = repr_method.electronic.method if repr_method.electronic else None
            if method_id is None:
                method_id = config.services.unavailable_value
            elif method_id in {"G0W0", "scGW"}:
                method_id = "GW"
        elif n_methods > 1:
            for sec_method in methods:
                # GW
                if sec_method.gw is not None:
                    repr_method = sec_method
                    method_id = "GW"
                    break

                # Methods linked to each other through references. Get all
                # linked methods, try to get electronic_structure_method from
                # each.
                linked_methods = [sec_method]
                if sec_method.core_method_ref is not None:
                    linked_methods.append(sec_method.core_method_ref)
                for i_method in linked_methods:
                    electronic_structure_method = i_method.electronic.method if i_method.electronic else None
                    if electronic_structure_method is not None:
                        repr_method = sec_method
                        method_id = electronic_structure_method
        method.method_type = method_id
        return repr_method, method_id

    def fill(self, context: Context):
        # Fill structure related metainfo
        struct: Any = None
        if context.material_type == Material.material_type.type.bulk:
            struct = MaterialBulkNormalizer(self.entry_archive, self.logger)
        elif context.material_type == Material.material_type.type.two_d:
            struct = Material2DNormalizer(self.entry_archive, self.logger)
        elif context.material_type == Material.material_type.type.one_d:
            struct = Material1DNormalizer(self.entry_archive, self.logger)
        if struct is not None:
            struct.normalize(context)

        # Fill method related metainfo
        method = None
        if context.method_type == Method.method_type.type.DFT or context.method_type == Method.method_type.type.DFTU:
            method = MethodDFTNormalizer(self.entry_archive, self.logger)
        elif context.method_type == Method.method_type.type.GW:
            method = MethodGWNormalizer(self.entry_archive, self.logger)
        if method is not None:
            method.normalize(context)

        # Fill properties related metainfo
        properties = PropertiesNormalizer(self.entry_archive, self.logger)
        properties.normalize(context)

    def normalize(self, logger=None) -> None:
        """The caller will automatically log if the normalizer succeeds or ends
        up with an exception.
        """
        if self.entry_archive.metadata is None:
            self.entry_archive.m_create(EntryMetadata)
        sec_enc = self.entry_archive.metadata.m_create(EncyclopediaMetadata)
        status_enums = EncyclopediaMetadata.status.type
        calc_enums = Calculation.calculation_type.type

        # Do nothing if section_run is not present
        if self.section_run is None:
            status = status_enums.invalid_metainfo
            sec_enc.status = status
            self.logger.info(
                "required metainfo is missing or is invalid.",
                enc_status=status,
                invalid_metainfo="section_run",
            )
            return

        try:
            super().normalize(logger)
            # Initialise metainfo structure
            material = sec_enc.m_create(Material)
            method = sec_enc.m_create(Method)
            sec_enc.m_create(Properties)
            calc = sec_enc.m_create(Calculation)

            # Determine run type, stop if unknown
            calc_type = self.calc_type(calc)
            if calc_type == config.services.unavailable_value:
                status = status_enums.unsupported_calculation_type
                sec_enc.status = status
                self.logger.info(
                    "unsupported calculation type for encyclopedia",
                    enc_status=status,
                )
                return

            # Get the system type.
            material_enums = Material.material_type.type
            representative_system, material_type = self.material_type(material)
            if material_type != material_enums.bulk and material_type != material_enums.two_d and material_type != material_enums.one_d:
                status = status_enums.unsupported_material_type
                sec_enc.status = status
                self.logger.info(
                    "unsupported material type for encyclopedia",
                    enc_status=status,
                )
                return

            # Get the method type. For now, we allow unknown method type for
            # phonon calculations, as the method information is resolved at a
            # later stage.
            representative_method, method_type = self.method_type(method)
            if method_type == config.services.unavailable_value:
                sec_enc.status = status_enums.unsupported_method_type
                self.logger.info(
                    "unsupported method type for encyclopedia",
                    enc_status=status_enums.unsupported_method_type,
                )
                if calc_type != calc_enums.phonon_calculation:
                    return

            # Get representative scc
            try:
                representative_scc_idx = self.section_run.m_cache["representative_scc_idx"]
                representative_scc = self.section_run.calculation[representative_scc_idx]
            except Exception:
                representative_scc = None
                representative_scc_idx = None

            # Create one context that holds all details
            context = Context(
                material_type=material_type,
                method_type=method_type,
                calc_type=calc_type,
                representative_system=representative_system,
                representative_method=representative_method,
                representative_scc=representative_scc,
                representative_scc_idx=representative_scc_idx,
            )

            # Put the encyclopedia section into archive
            self.fill(context)

            # Check that the necessary information is in place
            functional_type = method.functional_type
            if functional_type is None:
                sec_enc.status = status_enums.unsupported_method_type
                self.logger.info(
                    "unsupported functional type for encyclopedia",
                    enc_status=status_enums.unsupported_method_type,
                )
                return

        except Exception:
            status = status_enums.failure
            sec_enc.status = status
            self.logger.error(
                "failed to process encyclopedia data due to an unhandlable exception",
                enc_status=status,
            )
            raise  # Reraise for the caller to log the exception as well
        else:
            status = status_enums.success
            sec_enc.status = status
            self.logger.info(
                "successfully created metainfo for encyclopedia.",
                enc_status=status,
            )
