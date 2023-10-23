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
import { renderNoAPI, screen } from '../conftest.spec'
import XYPlot from './XYPlot'

const subSectionDef = {
  m_def: 'nomad.metainfo.metainfo.Section',
  more: {
    label_quantity: 'some_label_quantity'
  },
  _properties: {
    name: {
      type: {
        type_data: 'str',
        type_kind: 'python'
      }
    },
    some_label_quantity: {
      type: {
        type_data: 'str',
        type_kind: 'python'
      }
    },
    x: {
      type: {
        type_data: 'float64',
        type_kind: 'numpy'
      },
      shape: ['*']
    },
    y: {
      type: {
        type_data: 'float64',
        type_kind: 'numpy'
      },
      shape: ['*']
    }
  }
}

const middleSectionDef = {
  m_def: 'nomad.metainfo.metainfo.Section',
  _properties: {
    name: {
      type: {
        type_data: 'str',
        type_kind: 'python'
      }
    },
    my_sub_sections: {
      m_def: 'nomad.metainfo.metainfo.SubSection',
      repeats: true,
      sub_section: subSectionDef
    },
    x: {
      type: {
        type_data: 'float64',
        type_kind: 'numpy'
      },
      shape: []
    },
    y: {
      type: {
        type_data: 'float64',
        type_kind: 'numpy'
      },
      shape: []
    },
    y2: {
      type: {
        type_data: 'float64',
        type_kind: 'numpy'
      },
      shape: []
    }
  }
}

const topSectionDef = {
  m_def: 'nomad.metainfo.metainfo.Section',
  _properties: {
    my_middle_sections: {
      m_def: 'nomad.metainfo.metainfo.SubSection',
      repeats: true,
      sub_section: middleSectionDef
    }
  }
}

const topSectionNamed = {
  my_middle_sections: [
    {
      name: 'Middle 1',
      my_sub_sections: [
        {
          name: 'Sub 1',
          x: [0, 1],
          y: [1, 2]
        }, {
          name: 'Sub 2',
          x: [0, 0.5],
          y: [1, 2]
        }, {
          name: 'Sub 3',
          x: [0.2, 0.4],
          y: [1.5, 3]
        }
      ]
    }, {
      name: 'Middle 2',
      my_sub_sections: [
        {
          name: 'Sub 1',
          x: [0.5, 1],
          y: [3, 1]
        }, {
          name: 'Sub 2',
          x: [0.5, 1.5],
          y: [-1, 2]
        }, {
          name: 'Sub 3',
          x: [0.5, 1.5],
          y: [-1, 2]
        }
      ]
    }
  ]
}

const topSectionUnNamed = {
  my_middle_sections: [
    {
      my_sub_sections: [
        {
          x: [0, 1],
          y: [1, 2]
        }, {
          x: [0, 0.5],
          y: [1, 2]
        }, {
          x: [0.2, 0.4],
          y: [1.5, 3]
        }
      ]
    }, {
      my_sub_sections: [
        {
          x: [0.5, 1],
          y: [3, 1]
        }, {
          x: [0.5, 1.5],
          y: [-1, 2]
        }, {
          x: [0.5, 1.5],
          y: [-1, 2]
        }
      ]
    }
  ]
}

const topSectionLabeled = {
  my_middle_sections: [
    {
      name: 'Middle 1',
      my_sub_sections: [
        {
          some_label_quantity: 'Sub 1',
          x: [0, 1],
          y: [1, 2]
        }, {
          name: 'Sub 2',
          x: [0, 0.5],
          y: [1, 2]
        }, {
          x: [0.2, 0.4],
          y: [1.5, 3]
        }
      ]
    }
  ]
}

const topSectionScalar = {
  my_middle_sections: [
    {
      x: 0,
      y: 2,
      y2: 3
    },
    {
      x: 1,
      y: 3,
      y2: 4
    }
  ]
}

test.each([
  [
    'Plotting first two sub sections',
    {
      x: 'my_middle_sections/:/my_sub_sections/:2/x',
      y: 'my_middle_sections/:/my_sub_sections/:2/y'
    },
    topSectionNamed,
    topSectionDef,
    ['Middle 1, Sub 1, Y', 'Middle 1, Sub 2, Y', 'Middle 2, Sub 1, Y', 'Middle 2, Sub 2, Y']
  ],
  [
    'Plotting last two sub sections',
    {
      x: 'my_middle_sections/:/my_sub_sections/-2:/x',
      y: 'my_middle_sections/:/my_sub_sections/-2:/y'
    },
    topSectionNamed,
    topSectionDef,
    ['Middle 1, Sub 2, Y', 'Middle 1, Sub 3, Y', 'Middle 2, Sub 2, Y', 'Middle 2, Sub 3, Y']
  ],
  [
    'Plotting unnamed sub sections with slice',
    {
      x: 'my_middle_sections/:/my_sub_sections/1:2/x',
      y: 'my_middle_sections/:/my_sub_sections/1:2/y'
    },
    topSectionUnNamed,
    topSectionDef,
    ['My Middle Sections 0, My Sub Sections 1, Y', 'My Middle Sections 1, My Sub Sections 1, Y']
  ],
  [
    'Plotting only instantiated sub sections',
    {
      x: 'my_middle_sections/:/my_sub_sections/2:5/x',
      y: 'my_middle_sections/:/my_sub_sections/2:5/y'
    },
    topSectionNamed,
    topSectionDef,
    ['Middle 1, Sub 3, Y', 'Middle 2, Sub 3, Y']
  ],
  [
    'Plotting labeled, named and unlabeled',
    {
      x: 'my_middle_sections/:/my_sub_sections/:/x',
      y: 'my_middle_sections/:/my_sub_sections/:/y'
    },
    topSectionLabeled,
    topSectionDef,
    ['Middle 1, Sub 1, Y', 'Middle 1, Sub 2, Y', 'Middle 1, My Sub Sections 2, Y']
  ],
  [
    'Plotting scalars from repeating sub sections',
    {
      x: 'my_middle_sections/:/x',
      y: [
        'my_middle_sections/:/y',
        'my_middle_sections/:/y2'
      ]
    },
    topSectionScalar,
    topSectionDef,
    ['My Middle Sections, Y', 'My Middle Sections, Y2']
  ]
])('Render XYPlot: %s', async (title, plot, section, sectionDef, expected_texts) => {
  renderNoAPI(<XYPlot
      plot={plot}
      section={section}
      sectionDef={sectionDef}
      title={title}
    />)
  expected_texts.forEach(expected_text => {
    screen.getByText(expected_text)
  })
})
