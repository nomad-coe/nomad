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
import React, {useCallback} from 'react'
import PropTypes from 'prop-types'
import {Card, makeStyles, Box} from '@material-ui/core'
import {SourceJsonDialogButton} from '../buttons/SourceDialogButton'
import { AutocompleteEditQuantity, BoolEditQuantity, EnumEditQuantity, NumberEditQuantity, StringEditQuantity } from './EditQuantity'

const coatingMethods = [
  'Vapor deposition', 'Chemical vapor deposition', 'Metalorganic vapour phase epitaxy', 'Electrostatic spray assisted vapour deposition (ESAVD)', 'Sherardizing',
  'Some forms of Epitaxy', 'Molecular beam epitaxy', 'Physical vapor deposition', 'Cathodic arc deposition', 'Electron beam physical vapor deposition (EBPVD)',
  'Ion plating', 'Ion beam assisted deposition (IBAD)', 'Magnetron sputtering', 'Pulsed laser deposition', 'Sputter deposition',
  'Vacuum deposition', 'Pulsed electron deposition (PED)', 'Chemical and electrochemical techniques', 'Conversion coating', 'Autophoretic', 'Anodising',
  'Chromate conversion coating', 'Plasma electrolytic oxidation', 'Phosphate', 'Ion beam mixing', 'Pickled and oiled, a type of plate steel coating',
  'Plating', 'Electroless plating', 'Electroplating', 'Spraying', 'Spray painting', 'High velocity oxygen fuel (HVOF)', 'Plasma spraying',
  'Thermal spraying', 'Kinetic metallization (KM)', 'Plasma transferred wire arc thermal spraying', 'The common forms of Powder coating', 'Roll-to-roll coating processes',
  'Common roll-to-roll coating processes include:', 'Air knife coating', 'Anilox coater', 'Flexo coater', 'Gap Coating', 'Knife-over-roll coating', 'Gravure coating',
  'Immersion dip coating', 'Kiss coating', 'Metering rod (Meyer bar) coating', 'Roller coating', 'Forward roller coating', 'Reverse roll coating',
  'Silk Screen coater', 'Rotary screen', 'Lithography', 'Flexography', 'Physical coating processes', 'Langmuir-Blodgett', 'Spin coating', 'Dip coating']

const defs = {
  ShortStringQuantityDef: {
    name: 'name',
    type: {
      type_kind: 'python',
      type_data: 'str'
    },
    m_annotations: {
      'eln': [
        {
          label: 'Name',
          component: 'StringEditQuantity'
        }
      ]
    }
  },
  LongStringQuantityDef: {
    name: 'description',
    description: 'It is a description that explains the whole section',
    type: {
      type_kind: 'python',
      type_data: 'str'
    },
    m_annotations: {
      'eln': [
        {
          label: 'Description',
          component: 'StringEditQuantity',
          props: {
            multiline: true,
            minRows: 5
          }
        }
      ]
    }
  },
  floatQuantityDef: {
    name: 'height',
    unit: 'meter',
    description: 'The value of height in meter',
    default: 1.5e22,
    type: {
      type_kind: 'numpy',
      type_data: 'float64'
    },
    m_annotations: {
      'eln': [
        {
          label: 'Height',
          component: 'NumberEditQuantity'
        }
      ]
    }
  },
  limitedFloatQuantityDef: {
    name: 'mass',
    unit: 'kilogram',
    description: 'The mass in Kg',
    default: 0,
    type: {
      type_kind: 'numpy',
      type_data: 'float64'
    },
    m_annotations: {
      'eln': [
        {
          label: 'Mass',
          component: 'NumberEditQuantity'
        }
      ]
    }
  },
  integerQuantityDef: {
    name: 'spin',
    description: 'The spin',
    unit: 'meter',
    type: {
      type_kind: 'numpy',
      type_data: 'int64'
    },
    m_annotations: {
      'eln': [
        {
          label: 'Spin',
          component: 'NumberEditQuantity',
          props: {
            minValue: -2,
            maxValue: 2
          }
        }
      ]
    }
  },
  limitedIntegerQuantityDef: {
    name: 'count',
    description: 'The count',
    type: {
      type_kind: 'numpy',
      type_data: 'uint64'
    },
    m_annotations: {
      'eln': [
        {
          label: 'Count',
          component: 'NumberEditQuantity',
          props: {
            maxValue: 10
          }
        }
      ]
    }
  },
  EnumQuantityDef: {
    name: 'distance',
    description: 'The distance',
    default: '200',
    type: {
      type_kind: 'Enum',
      type_data: [
        '100',
        '200',
        '300'
      ]
    },
    m_annotations: {
      'eln': [
        {
          label: 'Distance',
          component: 'EnumEditQuantity'
        }
      ]
    }
  },
  EnumQuantityDef2: {
    name: 'polarized',
    description: 'Polarized',
    default: 'No',
    type: {
      type_kind: 'Enum',
      type_data: [
        'Yes',
        'No'
      ]
    },
    m_annotations: {
      'eln': [
        {
          label: 'Polarized',
          component: 'EnumEditQuantity'
        }
      ]
    }
  },
  AutocompleteQuantityDef: {
    name: 'coatingMethod',
    description: 'The coating method',
    default: 'Vapor deposition',
    type: {
      type_kind: 'Enum',
      type_data: coatingMethods
    },
    m_annotations: {
      'eln': [
        {
          label: 'Coating Method',
          component: 'AutocompleteEditQuantity'
        }
      ]
    }
  },
  BoolQuantityDef: {
    name: 'spinPolarized',
    description: 'The spin polarization',
    type: {
      type_kind: 'python',
      type_data: 'bool'
    },
    m_annotations: {
      'eln': [
        {
          label: 'Spin Polarized',
          component: 'BoolEditQuantity'
        }
      ]
    }
  }
}

let section = {
  name: 'unnamed',
  description: 'This is a section',
  count: 10,
  distance: 100,
  coatingMethod: 'Pulsed electron deposition (PED)',
  mass: 1
}

const useStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center'
  },
  card: {
    width: '600px',
    padding: '10px'
  }
}))

const EditQuantity = React.memo((props) => {
  const {quantityDef, section, onChange} = props

  const eln = quantityDef?.m_annotations?.eln
  const component = (eln.length > 0 ? eln[0]?.component : undefined)
  const label = (eln.length > 0 ? eln[0]?.label : quantityDef.name)
  const annotationsProps = (eln.length > 0 ? eln[0]?.props : undefined)

  let otherProps = {'label': label, ...annotationsProps}
  if (component === 'StringEditQuantity') {
    return <StringEditQuantity quantityDef={quantityDef} section={section} onChange={onChange} {...otherProps}/>
  } else if (component === 'NumberEditQuantity') {
    return <NumberEditQuantity quantityDef={quantityDef} section={section} onChange={onChange} {...otherProps}/>
  } else if (component === 'EnumEditQuantity') {
    return <EnumEditQuantity quantityDef={quantityDef} section={section} onChange={onChange} {...otherProps}/>
  } else if (component === 'AutocompleteEditQuantity') {
    return <AutocompleteEditQuantity quantityDef={quantityDef} section={section} onChange={onChange} {...otherProps}/>
  } else if (component === 'BoolEditQuantity') {
    return <BoolEditQuantity quantityDef={quantityDef} section={section} onChange={onChange} {...otherProps}/>
  }
})
EditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}

export function EditQuantityExamples() {
  const styles = useStyles()

  const handleChange = useCallback((value, section, quantityDef) => {
    section[quantityDef.name] = value
  }, [])

  return <div className={styles.root}>
    <Card className={styles.card}>
      <SourceJsonDialogButton title={`Section`} data={section}/>
      <Box margin={1}>
        <EditQuantity quantityDef={defs.ShortStringQuantityDef} section={section} onChange={handleChange}/>
      </Box>
      <Box margin={1}>
        <EditQuantity quantityDef={defs.LongStringQuantityDef} section={section} onChange={handleChange}/>
      </Box>
      <Box margin={1}>
        <EditQuantity quantityDef={defs.floatQuantityDef} section={section} onChange={handleChange}/>
      </Box>
      <Box margin={1}>
        <EditQuantity quantityDef={defs.limitedFloatQuantityDef} section={section} onChange={handleChange}/>
      </Box>
      <Box margin={1}>
        <EditQuantity quantityDef={defs.integerQuantityDef} section={section} onChange={handleChange}/>
      </Box>
      <Box margin={1}>
        <EditQuantity quantityDef={defs.limitedIntegerQuantityDef} section={section} onChange={handleChange}/>
      </Box>
      <Box margin={1}>
        <EditQuantity quantityDef={defs.EnumQuantityDef} section={section} onChange={handleChange}/>
      </Box>
      <Box margin={1}>
        <EditQuantity quantityDef={defs.EnumQuantityDef2} section={section} onChange={handleChange}/>
      </Box>
      <Box margin={1}>
        <EditQuantity quantityDef={defs.AutocompleteQuantityDef} section={section} onChange={handleChange}/>
      </Box>
      <Box margin={1}>
        <EditQuantity quantityDef={defs.BoolQuantityDef} section={section} onChange={handleChange}/>
      </Box>
    </Card>
  </div>
}
