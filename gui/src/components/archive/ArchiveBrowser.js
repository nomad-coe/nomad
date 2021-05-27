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
import React, { useMemo, useState, useCallback } from 'react'
import PropTypes from 'prop-types'
import { atom, useRecoilState, useRecoilValue } from 'recoil'
import { Box, FormGroup, FormControlLabel, Checkbox, TextField, Typography, makeStyles, Tooltip, FormControl, RadioGroup, Radio } from '@material-ui/core'
import { useRouteMatch, useHistory } from 'react-router-dom'
import Autocomplete from '@material-ui/lab/Autocomplete'
import Browser, { Item, Content, Compartment, List, Adaptor } from './Browser'
import { resolveRef, rootSections } from './metainfo'
import { Title, metainfoAdaptorFactory, DefinitionLabel } from './MetainfoBrowser'
import { Matrix, Number } from './visualizations'
import Structure from '../visualization/Structure'
import BrillouinZone from '../visualization/BrillouinZone'
import BandStructure from '../visualization/BandStructure'
import DOS from '../visualization/DOS'
import Markdown from '../Markdown'
import { getHighestOccupiedEnergy } from '../../utils'
import { convertSI, useUnits } from '../../units'
import { electronicRange } from '../../config'

export const configState = atom({
  key: 'config',
  default: {
    'showMeta': false,
    'showCodeSpecific': false,
    'showAllDefined': false,
    'energyUnit': 'joule'
  }
})

// Contains details about the currently visualized system. Used to detect if a
// reload is needed for the StructureViewer.
const visualizedSystem = {}

export default function ArchiveBrowser({data}) {
  const searchOptions = useMemo(() => archiveSearchOptions(data), [data])
  return (
    <Browser
      adaptor={archiveAdaptorFactory(data)}
      form={<ArchiveConfigForm searchOptions={searchOptions} />}
    />
  )
}
ArchiveBrowser.propTypes = ({
  data: PropTypes.object.isRequired
})

function ArchiveConfigForm({searchOptions}) {
  const [config, setConfig] = useRecoilState(configState)

  const handleConfigChange = event => {
    const changes = {[event.target.name]: event.target.checked}
    if (changes.showCodeSpecific) {
      changes.showAllDefined = !changes.showCodeSpecific
    } else if (changes.showAllDefined) {
      changes.showCodeSpecific = !changes.showAllDefined
    }
    setConfig({...config, ...changes})
  }

  const history = useHistory()
  const { url } = useRouteMatch()

  return (
    <Box marginTop={-3} padding={0}>
      <FormGroup row style={{alignItems: 'center'}}>
        <Box style={{width: 350, height: 60}}>
          <Autocomplete
            options={searchOptions}
            getOptionLabel={(option) => option.name}
            style={{ width: 350, marginTop: -20 }}
            onChange={(_, value) => {
              if (value) {
                history.push(url + value.path)
              }
            }}
            renderInput={(params) => <TextField {...params} label="search" margin="normal" />}
          />
        </Box>
        <Box flexGrow={1} />
        <Tooltip title="Enable to also show all code specific data">
          <FormControlLabel
            control={
              <Checkbox
                checked={config.showCodeSpecific}
                onChange={handleConfigChange}
                name="showCodeSpecific"
              />
            }
            label="code specific"
          />
        </Tooltip>
        <Tooltip title="Enable to also show metadata that is in principle available, but not within this entry">
          <FormControlLabel
            control={
              <Checkbox
                checked={config.showAllDefined}
                onChange={handleConfigChange}
                name="showAllDefined"
              />
            }
            label="all defined"
          />
        </Tooltip>
        <Tooltip title="Show the Metainfo definition on the bottom of each lane">
          <FormControlLabel
            control={
              <Checkbox
                checked={config.showMeta}
                onChange={handleConfigChange}
                name="showMeta" />
            }
            label="definitions"
          />
        </Tooltip>
      </FormGroup>
    </Box>
  )
}
ArchiveConfigForm.propTypes = ({
  searchOptions: PropTypes.arrayOf(PropTypes.object).isRequired
})

function archiveAdaptorFactory(data, sectionDef) {
  return new SectionAdaptor(data, sectionDef || rootSections.find(def => def.name === 'EntryArchive'), undefined, {archive: data})
}

function archiveSearchOptions(data) {
  const options = []
  const optionDefs = {}
  const optionKeys = {}
  function traverse(data, def, parentPath) {
    for (let key in data) {
      const childDef = def._properties[key]
      if (!childDef) {
        continue
      }

      const child = data[key]
      if (!child) {
        continue
      }

      let path = `${parentPath}/${key}`
      if (childDef.m_def === 'SubSection') {
        const sectionDef = resolveRef(childDef.sub_section)
        if (Array.isArray(child) && child.length > 0 && child[0]) {
          if (child.length > 1) {
            child.forEach((value, index) => traverse(value, sectionDef, `${path}:${index}`))
          } else {
            traverse(child[0], sectionDef, path)
          }
        } else {
          traverse(child, sectionDef, path)
        }
      }

      if (optionDefs[childDef._qualifiedName]) {
        continue
      }
      optionDefs[childDef._qualifiedName] = childDef

      const option = {
        name: key,
        data: data,
        def: childDef,
        path: path
      }
      options.push(option)

      if (optionKeys[key]) {
        const addPath = option => {
          const parents = option.path.split('/')
          const parent = parents[parents.length - 2].replace(/:[0-9]+$/, '')
          option.name += ` (${parent})`
        }
        if (!optionKeys[key].name.includes('(')) {
          addPath(optionKeys[key])
        }
        addPath(option)
      } else {
        optionKeys[key] = option
      }
    }
  }
  traverse(data, rootSections.find(def => def.name === 'EntryArchive'), '')
  return options
}

class ArchiveAdaptor extends Adaptor {
  constructor(obj, def, parent, context) {
    super(obj)
    this.def = def
    this.parent = parent
    this.context = context
  }

  adaptorFactory(obj, def, parent, context) {
    if (def.m_def === 'Section') {
      return new SectionAdaptor(obj, def, parent, context || this.context)
    } else if (def.m_def === 'Quantity') {
      return new QuantityAdaptor(obj, def, parent, context || this.context)
    }
  }

  itemAdaptor(key) {
    if (key === '_metainfo') {
      return metainfoAdaptorFactory(this.def)
    } else {
      throw new Error('Unknown item key')
    }
  }
}

class SectionAdaptor extends ArchiveAdaptor {
  itemAdaptor(key) {
    const [name, index] = key.split(':')
    const property = this.def._properties[name]
    const value = this.e[name]
    if (!property) {
      return super.itemAdaptor(key)
    } else if (property.m_def === 'SubSection') {
      const sectionDef = resolveRef(property.sub_section)
      if (property.repeats) {
        return this.adaptorFactory(value[parseInt(index || 0)], sectionDef, this.e)
      } else {
        return this.adaptorFactory(value, sectionDef, this.e)
      }
    } else if (property.m_def === 'Quantity') {
      // References: sections and quantities
      if (property.type.type_kind === 'reference') {
        if (property.shape.length === 0) {
          // some sections cannot be resolved, because they are not part of the archive
          // user_id->user is one example
          const resolved = resolveRef(value, this.context.archive) || {}
          const resolvedDef = resolveRef(property.type.type_data)
          if (resolvedDef.name === 'User' && !resolved.user_id) {
            resolved.user_id = value
          }
          return this.adaptorFactory(resolved, resolvedDef, this.e)
        } else {
        }
      }
      // Regular quantities
      return this.adaptorFactory(value, property, this.e)
    } else {
      throw new Error('Unknown metainfo meta definition')
    }
  }
  render() {
    return <Section section={this.e} def={this.def} parent={this.parent} />
  }
}

class QuantityAdaptor extends ArchiveAdaptor {
  render() {
    return <Quantity value={this.e} def={this.def} />
  }
}

function QuantityItemPreview({value, def}) {
  const units = useUnits()
  if (def.type.type_kind === 'reference') {
    return <Box component="span" fontStyle="italic">
      <Typography component="span">reference ...</Typography>
    </Box>
  }
  if (def.shape.length > 0) {
    const dimensions = []
    let current = value
    for (let i = 0; i < def.shape.length; i++) {
      dimensions.push(current.length)
      current = current[0]
    }
    let typeLabel
    if (def.type.type_kind === 'python') {
      typeLabel = 'list'
    } else {
      if (dimensions.length === 1) {
        typeLabel = 'vector'
      } else if (dimensions.length === 2) {
        typeLabel = 'matrix'
      } else {
        typeLabel = 'tensor'
      }
    }
    return <Box component="span" whiteSpace="nowrap" fontStyle="italic">
      <Typography component="span">
        {dimensions.map((dimension, index) => (
          <span key={index}>
            {index > 0 && <span>&nbsp;&times;&nbsp;</span>}{String(dimension)}
          </span>
        ))}&nbsp;{typeLabel}
      </Typography>
    </Box>
  } else {
    let finalValue = value
    let finalUnit = def.unit
    if (def.unit) {
      [finalValue, finalUnit] = convertSI(value, def.unit, units)
    }
    return <Box component="span" whiteSpace="nowarp">
      <Number component="span" variant="body1" value={finalValue} exp={8} />
      {finalUnit && <Typography component="span">&nbsp;{finalUnit}</Typography>}
    </Box>
  }
}
QuantityItemPreview.propTypes = ({
  value: PropTypes.any,
  def: PropTypes.object.isRequired
})

function QuantityValue({value, def}) {
  // Figure out the units
  const units = useUnits()
  let finalValue = value
  let finalUnit = def.unit
  if (def.unit) {
    [finalValue, finalUnit] = convertSI(value, def.unit, units)
  }

  return <Box
    marginTop={2} marginBottom={2} textAlign="center" fontWeight="bold"
  >
    {def.shape.length > 0 ? <Matrix values={finalValue} shape={def.shape} invert={def.shape.length === 1} /> : <Number value={finalValue} exp={16} variant="body2" />}
    {def.shape.length > 0 &&
      <Typography noWrap variant="caption">
        ({def.shape.map((dimension, index) => <span key={index}>
          {index > 0 && <span>&nbsp;&times;&nbsp;</span>}{String(dimension)}
        </span>)}&nbsp;)
      </Typography>
    }
    {def.unit && <Typography noWrap>{finalUnit}</Typography>}
  </Box>
}
QuantityValue.propTypes = ({
  value: PropTypes.any,
  def: PropTypes.object.isRequired
})

/**
 * An optional overview for a section displayed directly underneath the section
 * title.
 */
function Overview({section, def, parent}) {
  // States
  const [mode, setMode] = useState('bs')
  const units = useUnits()

  // Styles
  const useStyles = makeStyles(
    {
      bands: {
        width: '30rem',
        height: '30rem',
        margin: 'auto'
      },
      structure: {
        width: '28rem',
        margin: 'auto'
      },
      dos: {
        width: '20rem',
        height: '40rem',
        margin: 'auto'
      },
      radio: {
        display: 'flex',
        justifyContent: 'center'
      }
    }
  )
  const style = useStyles()

  const toggleMode = useCallback((event) => {
    setMode(event.target.value)
  }, [setMode])

  // Structure visualization for section_system
  if (def.name === 'System') {
    let url = window.location.href
    let name = 'section_system'
    let rootIndex = url.indexOf(name) + name.length
    let sectionPath = url.substring(0, rootIndex)
    let tmp = url.substring(rootIndex)
    let tmpIndex = tmp.indexOf('/')
    let index = tmpIndex === -1 ? tmp : tmp.slice(0, tmpIndex)

    // The section is incomplete, we leave the overview empty
    if (!section.atom_species) {
      return null
    }
    const nAtoms = section.atom_species.length
    let system = {
      'species': section.atom_species,
      'cell': section.lattice_vectors ? convertSI(section.lattice_vectors, 'meter', {length: 'angstrom'}, false) : undefined,
      'positions': convertSI(section.atom_positions, 'meter', {length: 'angstrom'}, false),
      'pbc': section.configuration_periodic_dimensions
    }
    visualizedSystem.sectionPath = sectionPath
    visualizedSystem.index = index
    visualizedSystem.nAtoms = nAtoms

    return <Structure
      aspectRatio={4 / 3}
      className={style.structure}
      data={system}
    ></Structure>
  // Structure visualization for idealized_structure
  } else if (def.name === 'IdealizedStructure') {
    // The section is incomplete, we leave the overview empty
    if (!section.atom_labels) {
      return null
    }
    const system = {
      species: section.atom_labels,
      cell: section.lattice_vectors ? convertSI(section.lattice_vectors, 'meter', {length: 'angstrom'}, false) : undefined,
      positions: section.atom_positions,
      fractional: true,
      pbc: section.periodicity
    }
    return <Structure
      data={system}
      className={style.structure}
      aspectRatio={1}>
    </Structure>
  // Band structure plot for section_k_band
  } else if (def.name === 'KBand') {
    return section.band_structure_kind !== 'vibrational'
      ? <>
        {mode === 'bs'
          ? <Box>
            <BandStructure
              className={style.bands}
              data={{
                energy_highest_occupied: getHighestOccupiedEnergy(section, parent),
                segments: section.section_k_band_segment,
                reciprocal_cell: section.reciprocal_cell
              }}
              layout={{yaxis: {autorange: false, range: convertSI(electronicRange, 'electron_volt', units, false)}}}
              aspectRatio={1}
              units={units}
            ></BandStructure>
          </Box>
          : <BrillouinZone
            className={style.bands}
            data={section}
            aspectRatio={1}
          ></BrillouinZone>
        }
        <FormControl component="fieldset" className={style.radio}>
          <RadioGroup row aria-label="position" name="position" defaultValue="bs" onChange={toggleMode} className={style.radio}>
            <FormControlLabel
              value="bs"
              control={<Radio color="primary" />}
              label="Band structure"
              labelPlacement="end"
            />
            <FormControlLabel
              value="bz"
              control={<Radio color="primary" />}
              label="Brillouin zone"
              labelPlacement="end"
            />
          </RadioGroup>
        </FormControl>
      </>
      : <Box>
        <BandStructure
          className={style.bands}
          data={{
            segments: section.section_k_band_segment,
            reciprocal_cell: section.reciprocal_cell
          }}
          aspectRatio={1}
          units={units}
          type='vibrational'
        ></BandStructure>
      </Box>
  // DOS plot for section_dos
  } else if (def.name === 'Dos') {
    const isVibrational = section.dos_kind === 'vibrational'
    const layout = isVibrational
      ? undefined
      : {yaxis: {autorange: false, range: convertSI(electronicRange, 'electron_volt', units, false)}}
    return <DOS
      className={style.dos}
      layout={layout}
      data={{
        energies: section.dos_energies_normalized,
        densities: section.dos_values_normalized,
        energy_highest_occupied: 0
      }}
      aspectRatio={1 / 2}
      units={units}
      type={isVibrational ? 'vibrational' : null}
    ></DOS>
  }
  return null
}
Overview.propTypes = ({
  def: PropTypes.object,
  section: PropTypes.object,
  parent: PropTypes.object
})

function Section({section, def, parent}) {
  const config = useRecoilValue(configState)

  if (!section) {
    console.error('section is not available')
    return ''
  }

  const filter = config.showCodeSpecific ? def => true : def => !def.name.startsWith('x_')
  return <Content>
    <Title def={def} data={section} kindLabel="section" />
    <Overview def={def} section={section} parent={parent}></Overview>
    <Compartment title="sub sections">
      {def.sub_sections
        .filter(subSectionDef => section[subSectionDef.name] || config.showAllDefined)
        .filter(filter)
        .map(subSectionDef => {
          const key = subSectionDef.name
          const disabled = section[key] === undefined
          if (!disabled && subSectionDef.repeats && section[key].length > 1) {
            return <List
              key={subSectionDef.name}
              itemKey={subSectionDef.name}
              title={subSectionDef.name} disabled={disabled}
            />
          } else {
            return <Item key={key} itemKey={key} disabled={disabled}>
              <Typography component="span">
                <Box fontWeight="bold" component="span">
                  {subSectionDef.name}
                </Box>
              </Typography>
            </Item>
          }
        })
      }
    </Compartment>
    <Compartment title="quantities">
      {def.quantities
        .filter(quantityDef => section[quantityDef.name] !== undefined || config.showAllDefined)
        .filter(filter)
        .map(quantityDef => {
          const key = quantityDef.name
          const disabled = section[key] === undefined
          return <Item key={key} itemKey={key} disabled={disabled}>
            <Box component="span" whiteSpace="nowrap">
              <Typography component="span">
                <Box fontWeight="bold" component="span">
                  {quantityDef.name}
                </Box>
              </Typography>{!disabled && <span>&nbsp;=&nbsp;<QuantityItemPreview value={section[quantityDef.name]} def={quantityDef} /></span>}
            </Box>
          </Item>
        })
      }
    </Compartment>
    <Meta def={def} />
  </Content>
}
Section.propTypes = ({
  section: PropTypes.object.isRequired,
  def: PropTypes.object.isRequired,
  parent: PropTypes.any
})

function Quantity({value, def}) {
  return <Content>
    <Title def={def} data={value} kindLabel="value" />
    <QuantityValue value={value} def={def} />
    <Meta def={def} />
  </Content>
}
Quantity.propTypes = ({
  value: PropTypes.any,
  def: PropTypes.object.isRequired
})

const useMetaStyles = makeStyles(theme => ({
  description: {
    marginTop: theme.spacing(1)
  },
  graph: {
    marginTop: theme.spacing(3)
  },
  metainfo: {
    marginBottom: theme.spacing(2)
  },
  metainfoItem: {
    fontWeight: 'bold'
  }
}))
export function Meta({def}) {
  const classes = useMetaStyles()
  const config = useRecoilValue(configState)
  if (!config.showMeta) {
    return ''
  }
  return <Compartment title="meta" color="primary">
    <div className={classes.metainfo}>
      <Item itemKey="_metainfo">
        <DefinitionLabel classes={{root: classes.metainfoItem}} def={def} isDefinition component="span" />
      </Item>
    </div>
    <Markdown classes={{root: classes.description}}>{def.description}</Markdown>
  </Compartment>
}
Meta.propTypes = ({
  def: PropTypes.object
})
