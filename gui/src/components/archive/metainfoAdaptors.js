import React from 'react'
import PropTypes from 'prop-types'
import { useRecoilValue } from 'recoil'
import Adaptor from './adaptors'
import { Item, Content, Compartment, viewConfigState, filterConfigState } from './ArchiveBrowser'
import { Typography, Box } from '@material-ui/core'
import Markdown from '../Markdown'
import { metainfoDef, resolveRef } from './metainfo'

export function metainfoAdaptorFactory(obj) {
  if (obj.m_def === 'Section') {
    return new SectionDefAdaptor(obj)
  } else if (obj.m_def === 'SubSection') {
    return new SubSectionDefAdaptor(obj)
  } else if (obj.m_def === 'Quantity') {
    return new QuantityDefAdaptor(obj)
  } else {
    throw new Error('Unknown metainfo definition type')
  }
}

class MetainfoAdaptor extends Adaptor {
  itemAdaptor(key) {
    if (key === '_metainfo') {
      return metainfoAdaptorFactory(metainfoDef(this.e.m_def))
    } else if (this.e[key]) {
      return metainfoAdaptorFactory(resolveRef(this.e[key]))
    } else {
      throw new Error('Unknown item key')
    }
  }
}

export class SectionDefAdaptor extends MetainfoAdaptor {
  itemAdaptor(key) {
    const property = this.e._properties[key]
    if (!property) {
      return super.itemAdaptor(key)
    } else {
      return metainfoAdaptorFactory(property)
    }
  }
  render() {
    return <SectionDef def={this.e} />
  }
}

class QuantityDefAdaptor extends MetainfoAdaptor {
  render() {
    return <QuantityDef def={this.e} />
  }
}

class SubSectionDefAdaptor extends MetainfoAdaptor {
  render() {
    return <SubSectionDef def={this.e} />
  }
}

function SectionDef({def}) {
  const filterConfig = useRecoilValue(filterConfigState)
  const filter = filterConfig.showCodeSpecific ? def => true : def => !def.name.startsWith('x_')
  return <Content style={{backgroundColor: 'grey'}}>
    <Compartment>
      <Definition def={def} isDefinition/>
    </Compartment>
    <Compartment title="sub section definitions">
      {def.sub_sections.filter(filter)
        .map(subSectionDef => {
          const key = subSectionDef.name
          return <Item key={key} itemKey={key}>
            <Typography component="span">
              <Box fontWeight="bold" component="span">
                {subSectionDef.name}
              </Box>
            </Typography>
          </Item>
        })
      }
    </Compartment>
    <Compartment title="quantity definitions">
      {def.quantities.filter(filter)
        .map(quantityDef => {
          const key = quantityDef.name
          return <Item key={key} itemKey={key}>
            <Box component="span" whiteSpace="nowrap">
              <Typography component="span">
                <Box fontWeight="bold" component="span">
                  {quantityDef.name}
                </Box>
              </Typography>
            </Box>
          </Item>
        })
      }
    </Compartment>
  </Content>
}
SectionDef.propTypes = ({
  def: PropTypes.object
})

function SubSectionDef({def}) {
  return <Content>
    <Compartment>
      <Definition def={def} isDefinition/>
    </Compartment>
    <Item itemKey="sub_section">
      <Typography component="span">
        <Box fontWeight="bold" component="span">
          section
        </Box>
      </Typography>: {resolveRef(def.sub_section).name}
    </Item>
  </Content>
}
SubSectionDef.propTypes = ({
  def: PropTypes.object
})

function QuantityDef({def}) {
  return <Content>
    <Definition def={def} isDefinition/>
  </Content>
}
QuantityDef.propTypes = ({
  def: PropTypes.object
})

const definitionLabels = {
  'Section': 'section',
  'Quantity': 'quantity',
  'SubSection': 'sub section'
}
export function Definition({def, isDefinition}) {
  const viewConfig = useRecoilValue(viewConfigState)
  const color = isDefinition ? 'primary' : 'initial'
  if (viewConfig.showDefinitions) {
    return <Box>
      <Typography color={color} variant="h6">{def.name}</Typography>
      <Typography color={color} variant="caption">{definitionLabels[def.m_def]}{isDefinition ? ' definition' : ''}</Typography>
      {def.description && def.description !== 'None' && viewConfig.showDescriptions &&
        <Box marginTop={1}>
          <Markdown>{def.description}</Markdown>
        </Box>
      }
      <Box marginTop={1}>
        <Item itemKey="_metainfo">
          <Typography>
            metainfo definition
          </Typography>
        </Item>
      </Box>
    </Box>
  } else {
    return ''
  }
}
Definition.propTypes = ({
  def: PropTypes.object,
  isDefinition: PropTypes.bool
})
