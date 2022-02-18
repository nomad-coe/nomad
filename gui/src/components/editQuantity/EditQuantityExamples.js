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
import {Card, makeStyles, Box} from '@material-ui/core'
import {EditQuantity} from './EditQuantity'
import {SourceJsonDialogButton} from '../buttons/SourceDialogButton'

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
    type: {
      type_kind: 'numpy',
      type_data: 'float64'
    },
    m_annotations: {
      'eln': [
        {
          label: 'Mass',
          component: 'NumberEditQuantity',
          props: {
            defaultValue: 0
          }
        }
      ]
    }
  },
  integerQuantityDef: {
    name: 'spin',
    description: 'The spin',
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
            maxValue: 2,
            defaultValue: 0
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
          component: 'EnumEditQuantity',
          props: {
            defaultValue: 'No'
          }
        }
      ]
    }
  }
}

let section = {
  mass: 100,
  name: 'unnamed',
  description: 'This is a section',
  count: 10,
  distance: 100
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
    </Card>
  </div>
}
