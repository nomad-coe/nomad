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
import { render, screen } from './conftest.spec'
import {CodeList} from './About'

test('list of codes renders correctly', async () => {
  render(<CodeList />)
  // Check that the list of parsers is printed correctly within the correct categories.
  // This list should be updated if new codes are added or the code names are updated.
  // This test value is hardcoded so that any unintentional changes to the code names and
  // categories can be avoided.
  const list = /Atomistic codes: ABACUS, ABINIT, AMS, ASAP, Amber, BOPfox, BigDFT, CASTEP, CHARMM, CP2K, CPMD, CRYSTAL, DFTB\+, DL_POLY, DMol3, Elk, FHI-aims, FLEUR, FPLO, GAMESS, GPAW, GROMACS, GROMOS, GULP, Gaussian, LAMMPS, MOPAC, Molcas, NAMD, NWChem, OCEAN, ONETEP, ORCA, Octopus, OpenMX, Psi4, Qball, Qbox, QuantumATK, QuantumESPRESSO, SIESTA, TURBOMOLE, Tinker, VASP, WIEN2k, Wannier90, YAMBO, exciting, libAtoms, soliddmft, w2dynamics, xTB, Workflow managers: AFLOW, ASR, Atomate, ElaStic, FHI-vibes, LOBSTER, MOFStructures, QuantumESPRESSOXSpectra, QuantumEspressPhonon, QuantumEspressoEPW, phonopy, Database managers: EELSDB, NeXus, OpenKIM$/
  expect(screen.getByTestId('code-list')).toHaveTextContent(list)
})
