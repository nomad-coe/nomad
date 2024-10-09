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
from nomad.datamodel.data import ArchiveSection
from nomad.metainfo import Quantity, Package

m_package = Package()


class ActionSection(ArchiveSection):
    """
    Use this section to define actions that can be triggered from the ELN
    interface.

    Usage:
    - In any of your schemas, use a a section that inherits from this one.
    - Provide a `description` and an eln annotation for the `action_trigger`
      quantity. To show a button that sets the trigger and performs a save
      operation, use the `ActionEditQuantity` component in the eln annotation .
    - Define a `perform_action` function that contains the code to trigger.

    Example Usage:
    ```python
    class MyAction(ActionSection):
        action_trigger = Quantity(
            description='Description for this action. Will be shown as a tooltip.',
            a_eln=dict(component='ActionEditQuantity', label='Button label'),
        )

        def perform_action(self, archive, logger):
            # First check that was this action requested. Note that you do not
            # have to modify this quantity: it will be set and reset
            # automatically.
            if self.action_trigger:
                # Perform your action here
                pass
    ```
    """

    action_trigger = Quantity(type=bool)

    def normalize(self, archive, logger):
        try:
            if self.action_trigger:
                self.perform_action(archive, logger)
        finally:
            self.action_trigger = False

    def perform_action(self, archive, logger):
        raise NotImplementedError(
            'No implementation provided for the perform_action function of an ActionSection.'
        )


m_package.__init_metainfo__()
