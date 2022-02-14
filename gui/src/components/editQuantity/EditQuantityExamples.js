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
import {Card, makeStyles, Box} from '@material-ui/core'
import {StringEditQantity} from './EditQuantity'

const defs = {
  quantityDef1: {
    name: 'Name',
    key: 'name'
  },
  quantityDef2: {
    name: 'Description',
    description: 'It is a description that explains the whole section',
    key: 'description'
  },
  quantityDef3: {
    name: 'Height',
    unit: 'meter',
    description: 'The value of height in meter',
    key: 'height'
  }
}

let section = {
  height: 100,
  name: 'unnamed',
  description: 'This is a section'
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

  return <div className={styles.root}>
    <Card className={styles.card}>
      <Box margin={1}>
        <StringEditQantity quantityDef={defs.quantityDef1} section={section}/>
      </Box>
      <Box margin={1}>
        <StringEditQantity quantityDef={defs.quantityDef2} section={section} multiline/>
      </Box>
      <Box margin={1}>
        <StringEditQantity quantityDef={defs.quantityDef3} section={section}/>
      </Box>
    </Card>
  </div>
}
