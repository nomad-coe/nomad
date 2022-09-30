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
import PropTypes from 'prop-types'
import {Card, Box, Typography, Grid, CardContent} from '@material-ui/core'
import {NumberEditQuantity} from './NumberEditQuantity'
import {StringEditQuantity, URLEditQuantity} from './StringEditQuantity'
import {EnumEditQuantity} from './EnumEditQuantity'
import {AutocompleteEditQuantity} from './AutocompleteEditQuantity'
import {RadioEnumEditQuantity} from './RadioEnumEditQuantity'
import {BoolEditQuantity} from './BoolEditQuantity'
import {SliderEditQuantity} from './SliderEditQuantity'
import { DateEditQuantity, DateTimeEditQuantity, TimeEditQuantity } from './DateTimeEditQuantity'
import RichTextEditQuantity from './RichTextEditQuantity'
import ListEditQuantity from './ListEditQuantity'
import { Code } from '../buttons/SourceDialogButton'
import { stripIndent } from '../../utils'
import AuthorEditQuantity from './AuthorEditQuantity'

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

function Example({code, children}) {
  return (
    <Box display="flex" flexDirection="row" alignItems="flext-start" marginBottom={2}>
      <Box width={500} marginRight={2}>
        {children}
      </Box>
      <Code
        code={stripIndent(code)}
      />
    </Box>
  )
}
Example.propTypes = {
  code: PropTypes.string,
  children: PropTypes.any
}

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

  const user = {
    type_kind: 'User',
    type_data: 'User'
  }

  return <Box margin={3}>
    <Grid container direction="row" spacing={2}>
      <Grid item>
        <Card>
          <CardContent>
            <Box width={1100}>
              <Grid container direction="column" spacing={1}>
                <Grid item>
                  <Example
                    code={`
                    string:
                      type: string
                      m_annotations:
                        eln:
                          component: StringEditQuantity`}
                  >
                    <StringEditQuantity {...createDefaultProps('string')} />
                  </Example>
                </Grid>
                <Grid item>
                  <Example
                    code={`
                    url_link:
                      type: string
                      m_annotations:
                        eln:
                          component: URLEditQuantity`}
                  >
                    <URLEditQuantity {...createDefaultProps('url')}/>
                  </Example>
                </Grid>
                <Grid item>
                  <Example
                    code={`
                    float:
                      type: np.float64
                      m_annotations:
                        eln:
                          component: NumberEditQuantity`}
                  >
                    <NumberEditQuantity {...createDefaultProps('float', {type: float})} />
                  </Example>
                </Grid>
                <Grid item>
                  <Example
                    code={`
                    float_unit:
                      type: np.float64
                      unit: second
                      m_annotations:
                        eln:
                          component: NumberEditQuantity
                          defaultDisplayUnit: "fs"`}
                  >
                    <NumberEditQuantity {...createDefaultProps('float_unit', {type: float, unit: 'second', m_annotations: {eln: [{defaultDisplayUnit: 'fs'}]}})}/>
                  </Example>
                </Grid>
                <Grid item>
                  <Example
                    code={`
                    float_derived_unit:
                      type: np.float64
                      unit: joule
                      m_annotations:
                        eln:
                          component: NumberEditQuantity
                          defaultDisplayUnit: "eV"`}
                  >
                    <NumberEditQuantity {...createDefaultProps('float_derived_unit', {type: float, unit: 'joule', m_annotations: {eln: [{defaultDisplayUnit: 'eV'}]}})}/>
                  </Example>
                </Grid>
                <Grid item>
                  <Example
                    code={`
                    float_complex_unit:
                      type: np.float64
                      unit: ampere / second^2 * meter
                      m_annotations:
                        eln:
                          component: NumberEditQuantity
                          defaultDisplayUnit: "milliampere / ms^2 * cm"`}
                  >
                    <NumberEditQuantity {...createDefaultProps('float_complex_unit', {type: float, unit: 'ampere / second^2 * meter', m_annotations: {eln: [{defaultDisplayUnit: 'milliampere / ms^2 * cm'}]}})}/>
                  </Example>
                </Grid>
                <Grid item>
                  <Example
                    code={`
                    float_with_bounds:
                      type: np.float64
                      unit: meter
                      m_annotations:
                        eln:
                          component: NumberEditQuantity
                          minValue: 0
                          maxValue: 10`}
                  >
                    <NumberEditQuantity
                      {...createDefaultProps('float_with_bounds', {type: float, unit: 'meter'})}
                      minValue={0} maxValue={10}
                    />
                  </Example>
                </Grid>
                <Grid item>
                  <Example
                    code={`
                    int:
                      type: int
                      m_annotations:
                        eln:
                          component: NumberEditQuantity`}
                  >
                    <NumberEditQuantity {...createDefaultProps('int', {type: int})} />
                  </Example>
                </Grid>
                <Grid item>
                  <Example
                    code={`
                    int_with_bounds:
                      type: int
                      m_annotations:
                        eln:
                          component: NumberEditQuantity
                          minValue: 0
                          maxValue: 10`}
                  >
                    <NumberEditQuantity
                      {...createDefaultProps('int_with_bounds', {type: int})}
                      minValue={0} maxValue={10}
                    />
                  </Example>
                </Grid>
                <Grid item>
                  <Example
                    code={`
                    bool:
                      type: bool
                      m_annotations:
                        eln:
                          component: BoolEditQuantity`}
                  >
                    <BoolEditQuantity {...createDefaultProps('bool')} />
                  </Example>
                </Grid>
                <Grid item>
                  <Example
                    code={`
                    radio_enum:
                      type:
                        type_kind: enum
                        type_data:
                          - one
                          - two
                          - three
                      m_annotations:
                        eln:
                          component: RadioEnumEditQuantity`}
                  >
                    <RadioEnumEditQuantity {...createDefaultProps('radio_enum', {
                      type: {
                        type_data: ['one', 'two', 'three']
                      }
                    })} />
                  </Example>
                </Grid>
                <Grid item>
                  <Example
                    code={`
                    select_enum:
                      type:
                        type_kind: enum
                        type_data:
                          - one
                          - two
                          - three
                      m_annotations:
                        eln:
                          component: EnumEditQuantity`}
                  >
                    <EnumEditQuantity {...createDefaultProps('select_enum', {
                      type: {
                        type_data: ['one', 'two', 'three']
                      }
                    })} />
                  </Example>
                </Grid>
                <Grid item>
                  <Example
                    code={`
                    string_with_suggestions:
                      type: str
                      m_annotations:
                        eln:
                          component: EnumEditQuantity
                          suggestions: ['one', 'two', 'three']`}
                  >
                    <EnumEditQuantity {...createDefaultProps('string_with_suggestions')} suggestions={['one', 'two', 'three']}/>
                  </Example>
                </Grid>
                <Grid item>
                  <Example
                    code={`
                    autocomplete_enum:
                      type:
                        type_kind: enum
                        type_data: [
                          ...
                          ]
                      m_annotations:
                        eln:
                          component: AutocompleteEditQuantity`}
                  >
                    <AutocompleteEditQuantity {...createDefaultProps('autocomplete_enum', {
                      type: {
                        type_data: enumValues
                      }
                    })} />
                  </Example>
                </Grid>
                <Grid item>
                  <Example
                    code={`
                    slider:
                      type: int
                      m_annotations:
                        eln:
                          component: SliderEditQuantity,
                          minValue: 0
                          maxValue: 10`}
                  >
                    <SliderEditQuantity {...createDefaultProps('slider')} minValue={0} maxValue={10} />
                  </Example>
                </Grid>
                <Grid item>
                  <Example
                    code={`
                    date_time:
                      type: nomad.metainfo.Datetime
                      m_annotations:
                        eln:
                          component: DateTimeEditQuantity`}
                  >
                    <DateTimeEditQuantity {...createDefaultProps('date_time')} />
                  </Example>
                </Grid>
                <Grid item>
                  <Example
                    code={`
                    date:
                      type: nomad.metainfo.Datetime
                      m_annotations:
                        eln:
                          component: DateEditQuantity`}
                  >
                    <DateEditQuantity {...createDefaultProps('date')} />
                  </Example>
                </Grid>
                <Grid item>
                  <Example
                    code={`
                    time:
                      type: nomad.metainfo.Datetime
                      m_annotations:
                        eln:
                          component: TimeEditQuantity`}
                  >
                    <TimeEditQuantity {...createDefaultProps('date')} />
                  </Example>
                </Grid>
                <Grid item>
                  <Example
                    code={`
                    rich_text:
                      type: str
                      m_annotations:
                        eln:
                          component: RichTextEditQuantity`}
                  >
                    <RichTextEditQuantity {...createDefaultProps('rich_text')} />
                  </Example>
                </Grid>
                <Grid item>
                  <Example
                    code={`
                    list:
                      type: str
                      shape: ['*']
                      m_annotations:
                        eln:
                          component: StringEditQuantity`}
                  >
                    <ListEditQuantity
                      component={StringEditQuantity}
                      {...createDefaultProps('list', {shape: ['*']})}
                    />
                  </Example>
                </Grid>
                <Grid item>
                  <Example
                    code={`
                    list_fixed:
                      type: str
                      shape: [3]
                      m_annotations:
                        eln:
                          component: StringEditQuantity`}
                  >
                    <ListEditQuantity
                      component={StringEditQuantity}
                      {...createDefaultProps('list_fixed', {shape: [3]})}
                    />
                  </Example>
                </Grid>
                <Grid item>
                  <Example
                    code={`
                    author:
                      type: User
                      m_annotations:
                        eln:
                          component: AuthorEditQuantity`}
                  >
                    <AuthorEditQuantity {...createDefaultProps('User', {type: user})} />
                  </Example>
                </Grid>
                <Grid item>
                  <Example
                    code={`
                    author:
                      type: Author
                      m_annotations:
                        eln:
                          component: AuthorEditQuantity`}
                  >
                    <AuthorEditQuantity {...createDefaultProps('Author')} />
                  </Example>
                </Grid>
              </Grid>
            </Box>
          </CardContent>
        </Card>
      </Grid>
      <Grid item>
        <Box margin={2}>
          <Typography component="pre">
            {JSON.stringify(sectionRef.current, null, 4)}
          </Typography>
        </Box>
      </Grid>
    </Grid>
  </Box>
}
