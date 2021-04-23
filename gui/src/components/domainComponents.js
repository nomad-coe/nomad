/*
 * Copyright The NOMAD Authors.
 *
 * This file is part of NOMAD. See https://nomad-lab.eu for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
import DFTEntryDetails from './dft/DFTEntryDetails'
import DFTEntryOverview from './dft/DFTEntryOverview'
import DFTEntryRawView from './dft/DFTEntryRawView'
import EMSEntryDetails from './ems/EMSEntryDetails'
import EMSEntryRawView from './ems/EMSEntryRawView'
import EMSEntryOverview from './ems/EMSEntryOverview'
import QCMSEntryDetails from './qcms/QCMSEntryDetails'
import QCMSEntryOverview from './qcms/QCMSEntryOverview'
import QCMSEntryRawView from './qcms/QCMSEntryRawView'

// Maps the components for each domain. These components are stored
export const domainComponents = ({
  dft: {
    /**
     * A component to render the domain specific quantities in the metadata card of
     * the entry view. Needs to work with props: data (the entry data from the API),
     * loading (a bool with api loading status).
     */
    EntryDetails: DFTEntryDetails,
    /**
     * Determines the layout of the overview page.
     */
    EntryOverview: DFTEntryOverview,
    /**
     * A component to render additional domain specific cards in the
     * the entry view. Needs to work with props: data (the entry data from the API),
     * loading (a bool with api loading status).
     */
    EntryRawView: DFTEntryRawView,
    /**
     * A component to render additional domain specific cards in the
     * the entry view. Needs to work with props: data (the entry data from the API),
     * loading (a bool with api loading status).
     */
    searchTabs: ['entries', 'materials', 'datasets', 'groups', 'uploads']
  },
  ems: {
    /**
     * A component to render the domain specific quantities in the metadata card of
     * the entry view. Needs to work with props: data (the entry data from the API),
     * loading (a bool with api loading status).
     */
    EntryDetails: EMSEntryDetails,
    /**
     * Determines the layout of the overview page.
     */
    EntryOverview: EMSEntryOverview,
    /**
     * A component to render additional domain specific cards in the
     * the entry view. Needs to work with props: data (the entry data from the API),
     * loading (a bool with api loading status).
     */
    EntryRawView: EMSEntryRawView,
    /**
     * Names of the possible search tabs for this domain
     */
    searchTabs: ['entries', 'datasets', 'uploads']
  },
  qcms: {
    /**
     * A component to render the domain specific quantities in the metadata card of
     * the entry view. Needs to work with props: data (the entry data from the API),
     * loading (a bool with api loading status).
     */
    EntryDetails: QCMSEntryDetails,
    /**
     * Determines the layout of the overview page.
     */
    EntryOverview: QCMSEntryOverview,
    /**
     * A component to render additional domain specific cards in the
     * the entry view. Needs to work with props: data (the entry data from the API),
     * loading (a bool with api loading status).
     */
    EntryRawView: QCMSEntryRawView,
    /**
     * Names of the possible search tabs for this domain
     */
    searchTabs: ['entries', 'datasets', 'uploads']
  }
})
