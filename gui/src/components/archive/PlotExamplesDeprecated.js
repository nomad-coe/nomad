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
import React, { useCallback, useMemo, useState } from 'react'
import PropTypes from 'prop-types'
import { Card, Box, CardContent } from '@material-ui/core'
import { Code } from '../buttons/SourceDialogButton'
import { stripIndent } from '../../utils'
import { SectionPlots } from './ArchiveBrowser'
import { useGlobalMetainfo } from './metainfo'
import InputConfig from '../search/input/InputConfig'
import Typography from '@material-ui/core/Typography'

const pvdSection = {
  name: 'Substrate 1',
  process_time: [-2099, -2026, -1953, -1880, -1807, -1734, -1661, -1588, -1515, -1442, -1369, -1296, -1223, -1150, -1077, -1004, -931, -858, -785, -712, -639, -566, -493, -419, -346, -273, -200, -127, -54, 19, 92, 165, 238, 311, 384, 457, 530, 603, 676, 749, 822, 895, 968, 1041, 1114, 1187, 1260, 1333, 1406, 1479, 1552, 1625, 1698, 1771, 1844, 1917, 1990, 2063, 2136, 2209, 2282, 2355, 2428, 2501, 2574, 2647, 2720, 2793, 2866, 2939, 3012, 3085, 3158, 3231, 3304, 3377, 3450, 3523, 3596, 3669, 3742, 3815, 3888, 3961, 4034, 4107, 4180, 4253, 4326, 4399, 4472, 4545, 4618, 4691, 4764, 4837, 4910, 4983, 5056, 5129, 5202, 5275, 5348, 5421, 5494, 5567, 5640, 5713, 5786, 5859, 5932, 6005, 6078, 6151, 6225, 6298, 6371, 6444, 6517, 6590, 6663, 6736, 6809, 6882, 6955, 7028, 7101, 7174, 7247, 7320, 7393, 7466, 7539],
  set_substrate_temperature: [309.356, 309.217, 309.084, 308.958, 308.838, 308.724, 308.616, 308.513, 308.415, 308.322, 308.233, 308.148, 308.068, 307.992, 307.919, 307.85, 307.784, 307.722, 307.662, 307.606, 307.552, 307.501, 307.452, 307.406, 307.362, 307.32, 307.28, 307.242, 307.206, 307.171, 307.139, 307.108, 307.078, 307.05, 307.023, 306.998, 306.974, 306.951, 306.929, 306.908, 306.888, 306.869, 306.851, 306.834, 306.818, 306.803, 306.788, 306.774, 306.761, 306.748, 306.736, 306.725, 306.714, 306.704, 306.694, 306.684, 306.676, 306.667, 306.659, 306.651, 306.644, 306.637, 306.631, 306.624, 306.618, 306.613, 306.607, 306.602, 306.597, 306.593, 306.588, 306.584, 306.58, 306.576, 306.573, 306.569, 306.566, 306.563, 306.56, 306.557, 306.554, 306.552, 306.549, 306.547, 306.545, 306.543, 306.541, 306.539, 306.537, 306.535, 306.534, 306.532, 306.531, 306.529, 306.528, 306.527, 306.525, 306.524, 306.523, 306.522, 306.521, 306.52, 306.519, 306.519, 306.518, 306.517, 306.516, 306.516, 306.515, 306.514, 306.514, 306.513, 306.513, 306.512, 306.512, 306.511, 306.511, 306.51, 306.51, 306.509, 306.509, 306.509, 306.508, 306.508, 306.508, 306.508, 306.507, 306.507, 306.507, 306.507, 306.506, 306.506, 306.506],
  substrate_temperature: [305.563, 305.527, 305.518, 305.523, 305.506, 305.502, 305.506, 305.515, 305.563, 305.575, 305.591, 305.607, 305.609, 305.637, 305.683, 305.701, 305.74, 305.788, 305.813, 305.871, 305.913, 305.977, 306.003, 306.076, 306.143, 306.179, 306.267, 306.321, 306.374, 306.453, 306.635, 306.871, 307.098, 307.324, 307.549, 307.801, 308.03, 308.264, 308.499, 308.728, 308.971, 309.16, 309.374, 309.6, 309.779, 309.95, 310.16, 310.334, 310.501, 310.668, 310.84, 310.999, 311.145, 311.303, 311.45, 311.58, 311.709, 311.852, 311.991, 312.093, 312.243, 312.336, 312.464, 312.569, 312.664, 312.785, 312.887, 312.973, 313.099, 313.192, 313.27, 313.369, 313.47, 313.568, 313.638, 313.718, 314.497, 317.131, 316.92, 316.736, 316.561, 316.389, 316.204, 316.021, 315.853, 315.692, 315.552, 315.409, 315.284, 315.124, 314.984, 314.851, 314.719, 314.583, 314.466, 314.338, 314.226, 314.085, 313.986, 313.85, 313.765, 313.639, 313.536, 313.44, 313.315, 313.211, 313.123, 313.004, 312.9, 312.801, 312.717, 312.634, 312.553, 312.468, 312.371, 312.257, 312.172, 312.085, 312.004, 311.921, 311.833, 311.75, 311.658, 311.583, 310.238, 310.188, 310.204, 310.208, 310.185, 309.659, 309.457, 309.572, 309.602],
  chamber_pressure: [0.313, 0.328, 0.319, 0.313, 0.307, 0.306, 0.304, 0.303, 0.304, 0.303, 0.303, 0.303, 0.303, 0.303, 0.303, 0.303, 0.302, 0.303, 0.301, 0.303, 0.303, 0.302, 0.301, 0.303, 0.303, 0.303, 0.303, 0.304, 0.303, 0.303, 0.303, 0.304, 0.304, 0.304, 0.304, 0.305, 0.304, 0.305, 0.304, 0.305, 0.305, 0.307, 0.305, 0.307, 0.308, 0.307, 0.31, 0.309, 0.307, 0.308, 0.314, 0.309, 0.307, 0.31, 0.311, 0.31, 0.312, 0.312, 0.314, 0.313, 0.315, 0.314, 0.314, 0.317, 0.314, 0.315, 0.319, 0.318, 0.317, 0.319, 0.317, 0.318, 0.317, 0.318, 0.318, 0.319, 0.319, 0.318, 0.319, 0.319, 0.319, 0.319, 0.319, 0.319, 0.319, 0.319, 0.319, 0.32, 0.318, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.321, 0.32, 0.32, 0.32, 0.322, 0.32, 0.32, 0.321, 0.32, 0.32, 0.321, 0.321, 0.321, 0.321, 0.321, 0.321, 0.322, 0.322, 0.323, 0.94, 2.25, 3.49, 4.28, 2.21, 0.72, 0.291, 47700.0, 102000.0, 102000.0, 102000.0, 102000.0, 102000.0, 102000.0, 102000.0, 3220.0]
}

const pvdSectionSecond = {
  name: 'Substrate 2',
  process_time: [-2099, -2026, -1953, -1880, -1807, -1734, -1661, -1588, -1515, -1442, -1369, -1296, -1223, -1150, -1077, -1004, -931, -858, -785, -712, -639, -566, -493, -419, -346, -273, -200, -127, -54, 19, 92, 165, 238, 311, 384, 457, 530, 603, 676, 749, 822, 895, 968, 1041, 1114, 1187, 1260, 1333, 1406, 1479, 1552, 1625, 1698, 1771, 1844, 1917, 1990, 2063, 2136, 2209, 2282, 2355, 2428, 2501, 2574, 2647, 2720, 2793, 2866, 2939, 3012, 3085, 3158, 3231, 3304, 3377, 3450, 3523, 3596, 3669, 3742, 3815, 3888, 3961, 4034, 4107, 4180, 4253, 4326, 4399, 4472, 4545, 4618, 4691, 4764, 4837, 4910, 4983, 5056, 5129, 5202, 5275, 5348, 5421, 5494, 5567, 5640, 5713, 5786, 5859, 5932, 6005, 6078, 6151, 6225, 6298, 6371, 6444, 6517, 6590, 6663, 6736, 6809, 6882, 6955, 7028, 7101, 7174, 7247, 7320, 7393, 7466, 7539],
  set_substrate_temperature: [309.356, 309.217, 309.084, 308.958, 308.838, 308.724, 308.616, 308.513, 308.415, 308.322, 308.233, 308.148, 308.068, 307.992, 307.919, 307.85, 307.784, 307.722, 307.662, 307.606, 307.552, 307.501, 307.452, 307.406, 307.362, 307.32, 307.28, 307.242, 307.206, 307.171, 307.139, 307.108, 307.078, 307.05, 307.023, 306.998, 306.974, 306.951, 306.929, 306.908, 306.888, 306.869, 306.851, 306.834, 306.818, 306.803, 306.788, 306.774, 306.761, 306.748, 306.736, 306.725, 306.714, 306.704, 306.694, 306.684, 306.676, 306.667, 306.659, 306.651, 306.644, 306.637, 306.631, 306.624, 306.618, 306.613, 306.607, 306.602, 306.597, 306.593, 306.588, 306.584, 306.58, 306.576, 306.573, 306.569, 306.566, 306.563, 306.56, 306.557, 306.554, 306.552, 306.549, 306.547, 306.545, 306.543, 306.541, 306.539, 306.537, 306.535, 306.534, 306.532, 306.531, 306.529, 306.528, 306.527, 306.525, 306.524, 306.523, 306.522, 306.521, 306.52, 306.519, 306.519, 306.518, 306.517, 306.516, 306.516, 306.515, 306.514, 306.514, 306.513, 306.513, 306.512, 306.512, 306.511, 306.511, 306.51, 306.51, 306.509, 306.509, 306.509, 306.508, 306.508, 306.508, 306.508, 306.507, 306.507, 306.507, 306.507, 306.506, 306.506, 306.506],
  substrate_temperature: [366.6756, 366.6324, 366.6216, 366.6276, 366.6072, 366.6024, 366.6072, 366.618, 366.6756, 366.69, 366.7092, 366.7284, 366.7308, 366.7644, 366.8196, 366.8412, 366.888, 366.9456, 366.9756, 367.0452, 367.0956, 367.1724, 367.2036, 367.2912, 367.3716, 367.4148, 367.5204, 367.5852, 367.6488, 367.7436, 367.962, 368.2452, 368.5176, 368.7888, 369.0588, 369.3612, 369.636, 369.9168, 370.1988, 370.4736, 370.7652, 370.992, 371.2488, 371.52, 371.7348, 371.94, 372.192, 372.4008, 372.6012, 372.8016, 373.008, 373.1988, 373.374, 373.5636, 373.74, 373.896, 374.0508, 374.2224, 374.3892, 374.5116, 374.6916, 374.8032, 374.9568, 375.0828, 375.1968, 375.342, 375.4644, 375.5676, 375.7188, 375.8304, 375.924, 376.0428, 376.164, 376.2816, 376.3656, 376.4616, 377.3964, 380.5572, 380.304, 380.0832, 379.8732, 379.6668, 379.4448, 379.2252, 379.0236, 378.8304, 378.6624, 378.4908, 378.3408, 378.1488, 377.9808, 377.8212, 377.6628, 377.4996, 377.3592, 377.2056, 377.0712, 376.902, 376.7832, 376.62, 376.518, 376.3668, 376.2432, 376.128, 375.978, 375.8532, 375.7476, 375.6048, 375.48, 375.3612, 375.2604, 375.1608, 375.0636, 374.9616, 374.8452, 374.7084, 374.6064, 374.502, 374.4048, 374.3052, 374.1996, 374.1, 373.9896, 373.8996, 372.2856, 372.2256, 372.2448, 372.2496, 372.222, 371.5908, 371.3484, 371.4864, 371.5224],
  chamber_pressure: [0.313, 0.328, 0.319, 0.313, 0.307, 0.306, 0.304, 0.303, 0.304, 0.303, 0.303, 0.303, 0.303, 0.303, 0.303, 0.303, 0.302, 0.303, 0.301, 0.303, 0.303, 0.302, 0.301, 0.303, 0.303, 0.303, 0.303, 0.304, 0.303, 0.303, 0.303, 0.304, 0.304, 0.304, 0.304, 0.305, 0.304, 0.305, 0.304, 0.305, 0.305, 0.307, 0.305, 0.307, 0.308, 0.307, 0.31, 0.309, 0.307, 0.308, 0.314, 0.309, 0.307, 0.31, 0.311, 0.31, 0.312, 0.312, 0.314, 0.313, 0.315, 0.314, 0.314, 0.317, 0.314, 0.315, 0.319, 0.318, 0.317, 0.319, 0.317, 0.318, 0.317, 0.318, 0.318, 0.319, 0.319, 0.318, 0.319, 0.319, 0.319, 0.319, 0.319, 0.319, 0.319, 0.319, 0.319, 0.32, 0.318, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.321, 0.32, 0.32, 0.32, 0.322, 0.32, 0.32, 0.321, 0.32, 0.32, 0.321, 0.321, 0.321, 0.321, 0.321, 0.321, 0.322, 0.322, 0.323, 0.94, 2.25, 3.49, 4.28, 2.21, 0.72, 0.291, 47700.0, 102000.0, 102000.0, 102000.0, 102000.0, 102000.0, 102000.0, 102000.0, 3220.0]
}

const processSection = {
  PVDEvaporation: [pvdSection, pvdSectionSecond]
}

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

const pvdExamples = [{
  a_plot: {
    label: 'T_substrate',
    x: 'process_time',
    y: './substrate_temperature'
  }
}, {
  a_plot: {
    label: 'Temperature',
    x: ['process_time', 'process_time'],
    y: ['./set_substrate_temperature', './substrate_temperature']
  }
}, {
  a_plot: {
    label: 'Temperature',
    x: 'process_time',
    y: ['./substrate_temperature', './chamber_pressure'],
    config: {
      editable: true,
      scrollZoom: false
    }
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

const processExamples = [{
    a_plot: {
      x: 'PVDEvaporation/0/process_time',
      y: 'PVDEvaporation/0/chamber_pressure'
    }
  }, {
    a_plot: {
      x: 'PVDEvaporation/1/process_time',
      y: 'PVDEvaporation/1/substrate_temperature'
    }
  }, {
    a_plot: {
      x: 'PVDEvaporation/0/process_time',
      y: ['PVDEvaporation/0/substrate_temperature', 'PVDEvaporation/2/substrate_temperature']
    }
  }, {
    a_plot: {
      x: 'PVDEvaporation/:2/process_time',
      y: 'PVDEvaporation/:2/substrate_temperature'
    }
  }

]

function getSection(metainfo, packageName, sectionName) {
  if (!metainfo) return undefined
  const packages = metainfo?._data?.packages
  const eln = packages?.find(section => section.name === packageName)
  const section_definitions = eln?.section_definitions
  return section_definitions?.find(section => section.name === sectionName)
}

export function PlotExamplesDeprecated() {
  const metainfo = useGlobalMetainfo()
  const [plots, setPlots] = useState(pvdExamples.concat(processExamples))
  const [keys, setKeys] = useState(pvdExamples.concat(processExamples).map((example, index) => `plot${index}-0`))

  const pvdSectionDef = useMemo(() => ({
    name: 'PVDEvaporation',
    _properties: {
      name: {
        m_def: 'nomad.metainfo.metainfo.Quantity'
      },
      process_time: {
        m_def: 'nomad.metainfo.metainfo.Quantity',
        unit: 'second',
        shape: ['*']
      },
      set_substrate_temperature: {
        m_def: 'nomad.metainfo.metainfo.Quantity',
        unit: 'kelvin',
        shape: ['*']
      },
      substrate_temperature: {
        m_def: 'nomad.metainfo.metainfo.Quantity',
        unit: 'kelvin',
        shape: ['*']
      },
      chamber_pressure: {
        m_def: 'nomad.metainfo.metainfo.Quantity',
        unit: 'pascal',
        shape: ['*']
      }
    }
  }), [])
  const archiveSectionDef = useMemo(() => metainfo ? getSection(metainfo, 'nomad.datamodel.data', 'ArchiveSection') : undefined, [metainfo])

  if (archiveSectionDef) {
    archiveSectionDef['sub_section'] = {
      m_def: 'nomad.metainfo.metainfo.Section',
      _properties: {
        PVDEvaporation: {
          m_def: 'nomad.metainfo.metainfo.SubSection',
          repeats: true,
          sub_section: pvdSectionDef
        }
      }
    }
  }

  const sectionDefs = useMemo(() => {
    return plots.map((example, index) => index < pvdExamples.length
      ? {...pvdSectionDef, m_annotations: {plot: example.a_plot}}
      : {...archiveSectionDef, m_annotations: {plot: example.a_plot}})
  }, [plots, pvdSectionDef, archiveSectionDef])

  const handleJsonChange = useCallback((data, index) => {
    const newPlots = [...plots]
    newPlots[index] = data
    setPlots(newPlots)
    const newKeys = [...keys]
    const [plotIndex, plotKey] = newKeys[index].split('-')
    newKeys[index] = `${plotIndex}-${Number(plotKey) + 1}`
    setKeys(newKeys)
  }, [keys, plots])

  return sectionDefs.map((def, index) => (
    <Box key={index}>
      {index === 0 && <Box margin={3} marginLeft={7}>
        <Typography variant="h6">
          Examples to show how to plot simple and multiline graphs and how to customize the line styles
        </Typography>
      </Box>}
      {index === pvdExamples.length && <Box margin={3} marginLeft={7}>
        <Typography variant="h6">
          Examples to plot data from an array of section
        </Typography>
      </Box>}
      {index === pvdExamples.length + 2 && <Box margin={3} marginLeft={7}>
        <Typography variant="h6">
          Non instantiated subsections are ignored
        </Typography>
      </Box>}
      {index === pvdExamples.length + 3 && <Box margin={3} marginLeft={7}>
        <Typography variant="h6">
          Subsection index can use python slice notation
        </Typography>
      </Box>}
      <Box display='flex' marginTop={3} width={'100%'} >
        <Box width={'50%'} marginLeft={6} marginRight={6}>
          <Card>
            <CardContent>
              <SectionPlots key={keys[index]} sectionDef={def} section={index < pvdExamples.length ? pvdSection : processSection}/>
            </CardContent>
          </Card>
        </Box>
        <Box width={'50%'} marginLeft={6} marginRight={6}>
          <InputConfig
            data={plots[index]}
            onChange={data => handleJsonChange(data, index)}
          />
        </Box>
      </Box>
    </Box>
  ))
}
