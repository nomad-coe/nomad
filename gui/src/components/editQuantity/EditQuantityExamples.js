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
import {FloatEditQuantity, IntegerEditQuantity, StringEditQuantity} from './EditQuantity'
import {SourceJsonDialogButton} from '../buttons/SourceDialogButton'

const defs = {
  ShortStringQuantityDef: {
    name: 'Name',
    key: 'name'
  },
  LongStringQuantityDef: {
    name: 'Description',
    description: 'It is a description that explains the whole section',
    key: 'description'
  },
  floatQuantityDef: {
    name: 'Height',
    unit: 'meter',
    description: 'The value of height in meter',
    key: 'height'
  },
  limitedFloatQuantityDef: {
    name: 'Mass',
    unit: 'kilogram',
    description: 'The mass in Kg',
    key: 'mass'
  },
  integerQuantityDef: {
    name: 'Spin',
    description: 'The spin',
    key: 'spin'
  },
  limitedIntegerQuantityDef: {
    name: 'Count',
    description: 'The count',
    key: 'count'
  }
}

let section = {
  height: 100,
  mass: 100,
  name: 'unnamed',
  description: 'This is a section',
  spin: 1,
  count: 10
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
    section[quantityDef.key] = value
  }, [])

  return <div className={styles.root}>
    <Card className={styles.card}>
      <SourceJsonDialogButton title={`Section`} data={section}/>
      <Box margin={1}>
        <StringEditQuantity quantityDef={defs.ShortStringQuantityDef} section={section} onChange={handleChange}/>
      </Box>
      <Box margin={1}>
        <StringEditQuantity quantityDef={defs.LongStringQuantityDef} section={section} onChange={handleChange} multiline/>
      </Box>
      <Box margin={1}>
        <FloatEditQuantity quantityDef={defs.floatQuantityDef} section={section} onChange={handleChange}/>
      </Box>
      <Box margin={1}>
        <FloatEditQuantity quantityDef={defs.limitedFloatQuantityDef} section={section} onChange={handleChange} minValue={0} defaultValue={0}/>
      </Box>
      <Box margin={1}>
        <IntegerEditQuantity quantityDef={defs.integerQuantityDef} section={section} onChange={handleChange} minValue={-2} maxValue={2} defaultValue={0}/>
      </Box>
      <Box margin={1}>
        <IntegerEditQuantity quantityDef={defs.limitedIntegerQuantityDef} section={section} onChange={handleChange} minValue={0} defaultValue={0}/>
      </Box>
    </Card>
  </div>
}
