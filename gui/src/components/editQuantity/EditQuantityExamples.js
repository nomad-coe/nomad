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
import React, {useRef, useState} from 'react'
import {Card, Box, Typography, Grid, CardContent} from '@material-ui/core'
import {NumberEditQuantity} from './NumberEditQuantity'
import {StringEditQuantity} from './StringEditQuantity'
import {EnumEditQuantity} from './EnumEditQuantity'
import {AutocompleteEditQuantity} from './AutocompleteEditQuantity'
import {RadioEnumEditQuantity} from './RadioEnumEditQuantity'
import {BoolEditQuantity} from './BoolEditQuantity'
import {SliderEditQuantity} from './SliderEditQuantity'
import { DateEditQuantity, DateTimeEditQuantity, TimeEditQuantity } from './DateTimeEditQuantity'
import RichTextEditQuantity from './RichTextEditQuantity'
import ListEditQuantity from './ListEditQuantity'

const enumValues = [
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

export function EditQuantityExamples() {
  const [, setUpdate] = useState(0)
  const propsRef = useRef({})
  const sectionRef = useRef({})

  const createDefaultProps = (name, props) => {
    if (!propsRef.current[name]) {
      propsRef.current[name] = {
        quantityDef: {
          name: name,
          description: `
            This is **MARKDOWN** help text.
          `,
          ...(props || {})
        },
        onChange: value => {
          sectionRef.current[name] = value
          setUpdate(value => value + 1)
        }
      }
    }
    return {
      ...(propsRef.current[name]),
      value: sectionRef.current[name]
    }
  }

  const float = {
    type_kind: 'python',
    type_data: 'float'
  }

  const int = {
    type_kind: 'python',
    type_data: 'int'
  }

  return <Box margin={3}>
    <Grid container direction="row" spacing={2}>
      <Grid item>
        <Card>
          <CardContent>
            <Box width={800}>
              <Grid container direction="column" spacing={1}>
                <Grid item>
                  <StringEditQuantity {...createDefaultProps('string')} />
                </Grid>
                <Grid item>
                  <NumberEditQuantity {...createDefaultProps('float', {type: float})} />
                </Grid>
                <Grid item>
                  <NumberEditQuantity {...createDefaultProps('float_unit', {type: float, unit: 'meter'})} defaultDisplayUnit={'bohr'}/>
                </Grid>
                <Grid item>
                  <NumberEditQuantity
                    {...createDefaultProps('float_with_bounds', {type: float, unit: 'second'})}
                    minValue={0} maxValue={10}
                  />
                </Grid>
                <Grid item>
                  <NumberEditQuantity {...createDefaultProps('int', {type: int})} />
                </Grid>
                <Grid item>
                  <NumberEditQuantity
                    {...createDefaultProps('int_with_bounds', {type: int})}
                    minValue={0} maxValue={10}
                  />
                </Grid>
                <Grid item>
                  <BoolEditQuantity {...createDefaultProps('bool')} />
                </Grid>
                <Grid item>
                  <RadioEnumEditQuantity {...createDefaultProps('radio_enum', {
                    type: {
                      type_data: ['one', 'two', 'three']
                    }
                  })} />
                </Grid>
                <Grid item>
                  <EnumEditQuantity {...createDefaultProps('select_enum', {
                    type: {
                      type_data: ['one', 'two', 'three']
                    }
                  })} />
                </Grid>
                <Grid item>
                  <EnumEditQuantity {...createDefaultProps('string with suggestions')} suggestions={['one', 'two', 'three']}/>
                </Grid>
                <Grid item>
                  <AutocompleteEditQuantity {...createDefaultProps('autocomplete_enum', {
                    type: {
                      type_data: enumValues
                    }
                  })} />
                </Grid>
                <Grid item>
                  <SliderEditQuantity {...createDefaultProps('slider')} minValue={0} maxValue={10} />
                </Grid>
                <Grid item>
                  <DateTimeEditQuantity {...createDefaultProps('date_time')} />
                </Grid>
                <Grid item>
                  <DateEditQuantity {...createDefaultProps('date')} />
                </Grid>
                <Grid item>
                  <TimeEditQuantity {...createDefaultProps('time')} />
                </Grid>
                <Grid item>
                  <RichTextEditQuantity {...createDefaultProps('richt_text')} />
                </Grid>
                <Grid item>
                  <ListEditQuantity
                    component={StringEditQuantity}
                    {...createDefaultProps('list', {shape: ['*']})}
                  />
                </Grid>
                <Grid item>
                  <ListEditQuantity
                    component={StringEditQuantity}
                    {...createDefaultProps('list_fixed', {shape: [3]})}
                  />
                </Grid>
              </Grid>
            </Box>
          </CardContent>
        </Card>
      </Grid>
      <Grid item>
        <Box margin={2} width={800}>
          <Typography component="pre">
            {JSON.stringify(sectionRef.current, null, 4)}
          </Typography>
        </Box>
      </Grid>
    </Grid>
  </Box>
}
