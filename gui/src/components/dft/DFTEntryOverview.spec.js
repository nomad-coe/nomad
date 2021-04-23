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

import React from 'react'
import 'regenerator-runtime/runtime'
import DFTEntryOverview from './DFTEntryOverview'
import { renderWithAPIRouter } from '../../testutils'
import { screen, waitForElementToBeRemoved } from '@testing-library/react'
import { within } from '@testing-library/dom'
import { repoDftBulk } from '../../../tests/repo'
import '@testing-library/jest-dom/extend-expect'

describe('<DFTEntryOverview />', () => {
  it('correctly renders basic data', async () => {
    const repo = repoDftBulk
    renderWithAPIRouter(
      <DFTEntryOverview
        data={repo}
      ></DFTEntryOverview>
    )

    // The data from repository should show up immediately, as it has been
    // loaded already.
    const code_name = screen.getByTitle('The name of the used program.') // We need to wait a bit for the component to render
    expect(within(code_name).getByText('program name')).toBeInTheDocument()
    expect(within(code_name).getByText(repo.results.method.simulation.program_name)).toBeInTheDocument()
    const code_version = screen.getByTitle('The version of the used program.')
    expect(within(code_version).getByText('program version')).toBeInTheDocument()
    expect(within(code_version).getByText(repo.results.method.simulation.program_version)).toBeInTheDocument()
    const xc_family = screen.getByTitle('The libXC based xc functional classification used in the simulation.')
    expect(within(xc_family).getByText('xc functional family')).toBeInTheDocument()
    expect(within(xc_family).getByText(repo.results.method.simulation.dft.xc_functional_type)).toBeInTheDocument()
    const xc_names = screen.getByTitle('The list of libXC functional names that where used in this entry.')
    expect(within(xc_names).getByText('xc functional names')).toBeInTheDocument()
    expect(within(xc_names).getByText(repo.results.method.simulation.dft.xc_functional_names.join(', '))).toBeInTheDocument()
    const basis_set_type = screen.getByTitle('The used basis set functions.')
    expect(within(basis_set_type).getByText('basis set type')).toBeInTheDocument()
    expect(within(basis_set_type).getByText(repo.results.method.simulation.dft.basis_set_type)).toBeInTheDocument()
    const comment = screen.getByTitle('A user provided comment for this entry')
    expect(within(comment).getByText('comment')).toBeInTheDocument()
    expect(within(comment).getByText(repo.comment)).toBeInTheDocument()
    const references = screen.getByTitle('User provided references (URLs) for this entry')
    expect(within(references).getByText('references')).toBeInTheDocument()
    expect(within(references).getByText(repo.references[0])).toBeInTheDocument()
    const authors = screen.getByTitle('All authors (uploader and co-authors)')
    expect(within(authors).getByText('authors')).toBeInTheDocument()
    expect(within(authors).getByText(repo.authors[0].name)).toBeInTheDocument()
    const datasets = screen.getByTitle('A list of user curated datasets this entry belongs to.')
    expect(within(datasets).getByText('datasets')).toBeInTheDocument()
    expect(within(datasets).getByText(repo.datasets[0].name)).toBeInTheDocument()
    const mainfile = screen.getByTitle('The path to the mainfile from the root directory of the uploaded files')
    expect(within(mainfile).getByText('mainfile')).toBeInTheDocument()
    expect(within(mainfile).getByText(repo.mainfile)).toBeInTheDocument()
    const entry_id = screen.getByTitle('The unique primary id for this entry.')
    expect(within(entry_id).getByText('entry id')).toBeInTheDocument()
    expect(within(entry_id).getByText(repo.calc_id)).toBeInTheDocument()
    const material_id = screen.getByTitle('A fixed length, unique material identifier in the form of a hash digest.')
    expect(within(material_id).getByText('material id')).toBeInTheDocument()
    expect(within(material_id).getByText(repo.encyclopedia.material.material_id)).toBeInTheDocument()
    const upload_id = screen.getByTitle('The persistent and globally unique identifier for the upload of the entry')
    expect(within(upload_id).getByText('upload id')).toBeInTheDocument()
    expect(within(upload_id).getByText(repo.upload_id)).toBeInTheDocument()
    const upload_time = screen.getByTitle('The date and time this entry was uploaded to nomad')
    expect(within(upload_time).getByText('upload time')).toBeInTheDocument()
    expect(within(upload_time).getByText(new Date(repo.upload_time).toLocaleString())).toBeInTheDocument()
    const last_processing = screen.getByTitle('The datetime of the last processing')
    expect(within(last_processing).getByText('last processing')).toBeInTheDocument()
    expect(within(last_processing).getByText(new Date(repo.last_processing).toLocaleString())).toBeInTheDocument()
    const processing_version = screen.getByTitle('Version used in the last processing')
    expect(within(processing_version).getByText('processing version')).toBeInTheDocument()
    expect(within(processing_version).getByText(`${repo.nomad_version}/${repo.nomad_commit}`)).toBeInTheDocument()
    const formula = screen.getByTitle('The chemical formula for a structure in Hill form with element symbols followed by integer chemical proportion numbers. The proportion number MUST be omitted if it is 1.')
    expect(within(formula).getByText('formula')).toBeInTheDocument()
    expect(within(formula).getByText(repo.results.material.chemical_formula_hill)).toBeInTheDocument()
    const structural_type = screen.getByTitle('Classification based on structural features.')
    expect(within(structural_type).getByText('structural type')).toBeInTheDocument()
    expect(within(structural_type).getByText(repo.results.material.type_structural)).toBeInTheDocument()
    const material_name = screen.getByTitle('Meaningful names for this a material if any can be assigned.')
    expect(within(material_name).getByText('material name')).toBeInTheDocument()
    expect(within(material_name).getByText(repo.results.material.material_name)).toBeInTheDocument()
    const crystal_system = screen.getByTitle('Name of the crystal system. Can be one of the following: triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal or cubic.')
    expect(within(crystal_system).getByText('crystal system')).toBeInTheDocument()
    expect(within(crystal_system).getByText(repo.results.material.symmetry.crystal_system)).toBeInTheDocument()
    const space_group = screen.getByTitle('Space group symbol and number')
    expect(within(space_group).getByText('space group')).toBeInTheDocument()
    expect(within(space_group).getByText(`${repo.results.material.symmetry.space_group_symbol} (${repo.results.material.symmetry.space_group_number})`)).toBeInTheDocument()
    const method_name = screen.getByTitle('Common name for the used method.')
    expect(within(method_name).getByText('method name')).toBeInTheDocument()
    expect(within(method_name).getByText('DFT')).toBeInTheDocument()

    // The test DOM does not support canvas or WebGL, and trying to add mocks
    // for them does not seem to work ATM. Thus we expect a message saying that
    // the 3D viewers are disabled.
    expect(screen.getByText('Could not display the visualization as your browser does not support WebGL content.')).toBeInTheDocument()

    // The data from archive should show placeholders until it is loaded
    const dos_electronic = screen.getByTestId('dos-electronic-overview')
    const dos_electronic_placeholder = screen.getByTestId('dos-electronic-overview-placeholder')
    expect(dos_electronic).toBeInTheDocument()
    expect(dos_electronic_placeholder).toBeInTheDocument()

    // When the archived data is loaded, check that the placeholders are removed
    // and the data is displayed correctly.
    waitForElementToBeRemoved(dos_electronic_placeholder).then(() => {
      expect(dos_electronic).not.toBeInTheDocument()
    })
  })
})
