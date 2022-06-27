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
import React, {useCallback, useEffect, useMemo, useState} from 'react'
import PropTypes from 'prop-types'
import {Card, Box, Grid, CardContent} from '@material-ui/core'
import { Code } from '../buttons/SourceDialogButton'
import { stripIndent } from '../../utils'
import {SectionPlots} from './ArchiveBrowser'
import {createGlobalMetainfo, resolveRef} from './metainfo'
import metainfoData from '../../metainfo.json'
import {JsonEditor} from './SectionEditor'

const section = {
  process_time: [-2099, -2026, -1953, -1880, -1807, -1734, -1661, -1588, -1515, -1442, -1369, -1296, -1223, -1150, -1077, -1004, -931, -858, -785, -712, -639, -566, -493, -419, -346, -273, -200, -127, -54, 19, 92, 165, 238, 311, 384, 457, 530, 603, 676, 749, 822, 895, 968, 1041, 1114, 1187, 1260, 1333, 1406, 1479, 1552, 1625, 1698, 1771, 1844, 1917, 1990, 2063, 2136, 2209, 2282, 2355, 2428, 2501, 2574, 2647, 2720, 2793, 2866, 2939, 3012, 3085, 3158, 3231, 3304, 3377, 3450, 3523, 3596, 3669, 3742, 3815, 3888, 3961, 4034, 4107, 4180, 4253, 4326, 4399, 4472, 4545, 4618, 4691, 4764, 4837, 4910, 4983, 5056, 5129, 5202, 5275, 5348, 5421, 5494, 5567, 5640, 5713, 5786, 5859, 5932, 6005, 6078, 6151, 6225, 6298, 6371, 6444, 6517, 6590, 6663, 6736, 6809, 6882, 6955, 7028, 7101, 7174, 7247, 7320, 7393, 7466, 7539],
  set_substrate_temperature: [323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15, 323.15],
  substrate_temperature: [305.563, 305.527, 305.518, 305.523, 305.506, 305.502, 305.506, 305.515, 305.563, 305.575, 305.591, 305.607, 305.609, 305.637, 305.683, 305.701, 305.74, 305.788, 305.813, 305.871, 305.913, 305.977, 306.003, 306.076, 306.143, 306.179, 306.267, 306.321, 306.374, 306.453, 306.635, 306.871, 307.098, 307.324, 307.549, 307.801, 308.03, 308.264, 308.499, 308.728, 308.971, 309.16, 309.374, 309.6, 309.779, 309.95, 310.16, 310.334, 310.501, 310.668, 310.84, 310.999, 311.145, 311.303, 311.45, 311.58, 311.709, 311.852, 311.991, 312.093, 312.243, 312.336, 312.464, 312.569, 312.664, 312.785, 312.887, 312.973, 313.099, 313.192, 313.27, 313.369, 313.47, 313.568, 313.638, 313.718, 314.497, 317.131, 316.92, 316.736, 316.561, 316.389, 316.204, 316.021, 315.853, 315.692, 315.552, 315.409, 315.284, 315.124, 314.984, 314.851, 314.719, 314.583, 314.466, 314.338, 314.226, 314.085, 313.986, 313.85, 313.765, 313.639, 313.536, 313.44, 313.315, 313.211, 313.123, 313.004, 312.9, 312.801, 312.717, 312.634, 312.553, 312.468, 312.371, 312.257, 312.172, 312.085, 312.004, 311.921, 311.833, 311.75, 311.658, 311.583, 310.238, 310.188, 310.204, 310.208, 310.185, 309.659, 309.457, 309.572, 309.602],
  chamber_pressure: [0.313, 0.328, 0.319, 0.313, 0.307, 0.306, 0.304, 0.303, 0.304, 0.303, 0.303, 0.303, 0.303, 0.303, 0.303, 0.303, 0.302, 0.303, 0.301, 0.303, 0.303, 0.302, 0.301, 0.303, 0.303, 0.303, 0.303, 0.304, 0.303, 0.303, 0.303, 0.304, 0.304, 0.304, 0.304, 0.305, 0.304, 0.305, 0.304, 0.305, 0.305, 0.307, 0.305, 0.307, 0.308, 0.307, 0.31, 0.309, 0.307, 0.308, 0.314, 0.309, 0.307, 0.31, 0.311, 0.31, 0.312, 0.312, 0.314, 0.313, 0.315, 0.314, 0.314, 0.317, 0.314, 0.315, 0.319, 0.318, 0.317, 0.319, 0.317, 0.318, 0.317, 0.318, 0.318, 0.319, 0.319, 0.318, 0.319, 0.319, 0.319, 0.319, 0.319, 0.319, 0.319, 0.319, 0.319, 0.32, 0.318, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.321, 0.32, 0.32, 0.32, 0.322, 0.32, 0.32, 0.321, 0.32, 0.32, 0.321, 0.321, 0.321, 0.321, 0.321, 0.321, 0.322, 0.322, 0.323, 0.94, 2.25, 3.49, 4.28, 2.21, 0.72, 0.291, 47700.0, 102000.0, 102000.0, 102000.0, 102000.0, 102000.0, 102000.0, 102000.0, 3220.0]}

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

const examples = [{
  a_plot: {
    label: 'T_substrate',
    x: 'process_time',
    y: './substrate_temperature'
  }
}, {
  a_plot: {
    label: 'Temperature',
    x: 'process_time',
    y: ['./set_substrate_temperature', './substrate_temperature'],
    config: {
      editable: true,
      scrollZoom: false
    }
  }
}, {
  a_plot: {
    x: 'process_time',
    y: ['./substrate_temperature', './chamber_pressure']
  }
}, {
  a_plot: [{
    label: 'Temperature and Pressure',
    x: 'process_time',
    y: ['./substrate_temperature', './chamber_pressure'],
    lines: [{
      mode: 'markers',
      marker: {
        color: 'rgb(40, 80, 130)'
      }}, {
      mode: 'lines',
      line: {
        color: 'rgb(100, 0, 0)'
      }}
    ]
  }, {
    label: 'Pressure of Chamber',
    x: 'process_time',
    y: 'chamber_pressure',
    layout: {
      xaxis: {title: 't (sec)'},
      yaxis: {title: 'P (GPa)', type: 'log'}
    }
  }]
}]

export function PlotExamples() {
  const [metainfo, setMetainfo] = useState()
  const [plots, setPlots] = useState(examples)

  useEffect(() => {
    createGlobalMetainfo(metainfoData).then(setMetainfo)
  }, [])

  const sectionDef = useMemo(() => metainfo ? resolveRef('/packages/7/section_definitions/4', metainfo._data) : undefined, [metainfo]) // Path for PVDEvaporation
  const exampleSectionDefs = useMemo(() => plots.map(example => ({...sectionDef, m_annotations: {plot: example.a_plot}})), [plots, sectionDef])

  const handleJsonChange = useCallback((data, index) => {
    const newPlots = [...plots]
    newPlots[index] = data
    setPlots(newPlots)
  }, [plots])

  return exampleSectionDefs.map((def, index) => (
    <Box key={index} display='flex' marginTop={3} width={'100%'} >
      <Box width={'50%'} marginLeft={6} marginRight={6}>
        <Card>
          <CardContent>
            <Grid item>
              <SectionPlots sectionDef={def} section={section}/>
            </Grid>
          </CardContent>
        </Card>
      </Box>
      <Box width={'50%'} marginLeft={6} marginRight={6}>
        <JsonEditor
          data={plots[index]}
          onChange={data => handleJsonChange(data, index)}
        />
      </Box>
    </Box>
  ))
}
