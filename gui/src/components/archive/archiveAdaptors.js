import React, { useState } from 'react'
import Adaptor from './adaptors'
import { Item, Content, Compartment } from './ArchiveBrowser'
import { Typography, Box, IconButton } from '@material-ui/core'
import { jsonAdaptorFactory } from './jsonAdaptors'
import Markdown from '../Markdown'
import MoreVertIcon from '@material-ui/icons/MoreVert'
import metainfo, { metainfoDef, resolveRef, sectionDefs } from './metainfo'
import { Definition, metainfoAdaptorFactory } from './metainfoAdaptors'

export default function archiveAdaptorFactory(data, sectionDef) {
  return new SectionAdaptor(data, sectionDef || sectionDefs['EntryArchive'], {archive: data})
}

class ArchiveAdaptor extends Adaptor {
  constructor(obj, def, context) {
    super(obj)
    this.def = def
    this.context = context
  }

  adaptorFactory(obj, def, context) {
    if (def.m_def === 'Section') {
      return new SectionAdaptor(obj, def, context || this.context)
    } else if (def.m_def === 'SubSection') {
      return new SubSectionAdaptor(obj, def, context || this.context)
    } else if (def.m_def === 'Quantity') {
      return new QuantityAdaptor(obj, def, context || this.context)
    }
  }

  itemAdaptor(key) {
    if (key === '_metainfo') {
      return metainfoAdaptorFactory(this.def)
    } else {
      throw new Error('Unknown item key')
    }
  }
}

export class SectionAdaptor extends ArchiveAdaptor {
  itemAdaptor(key) {
    const property = this.def._properties[key]
    const value = this.e[key]
    if (!property) {
      return super.itemAdaptor(key)
    } else if (property.m_def === 'SubSection') {
      const sectionDef = resolveRef(property.sub_section)
      if (Array.isArray(value)) {
        if (value.length === 1) {
          return this.adaptorFactory(value[0], sectionDef)
        } else {
          return this.adaptorFactory(value, property)
        }
      } else {
        return this.adaptorFactory(value, sectionDef)
      }
    } else if (property.m_def === 'Quantity') {
      if (property.type.type_kind === 'reference' && property.shape.length === 0) {
        return this.adaptorFactory(resolveRef(value, this.context.archive), resolveRef(property.type.type_data))
      }
      return this.adaptorFactory(value, property)
    } else {
      throw new Error('Unknown metainfo meta definition')
    }
  }
  render() {
    return <Section section={this.e} def={this.def} />
  }
}

class QuantityAdaptor extends ArchiveAdaptor {
  render() {
    return <Quantity value={this.e} def={this.def} />
  }
}

class SubSectionAdaptor extends ArchiveAdaptor {
  itemAdaptor(key) {
    const sectionDef = resolveRef(this.def.sub_section)
    return this.adaptorFactory(this.e[key], sectionDef)
  }
  render() {
    return <SubSection sections={this.e} def={this.def} />
  }
}

function QuantityItemPreview({value, def}) {
  if (def.type.type_kind === 'reference') {
    return <Box component="span" fontStyle="italic">
      <Typography component="span">reference ...</Typography>
    </Box>
  }
  if (def.shape.length > 0) {
    const dimensions = []
    let current = value
    for(let i = 0; i < def.shape.length; i++) {
      dimensions.push(current.length)
      current = current[0]
    }
    let typeLabel
    if (def.type.type_kind === 'python') {
      typeLabel = 'list'
    } else {
      if (dimensions.length === 1) {
        typeLabel = 'vector'
      } else if (dimensions.length === 2) {
        typeLabel = 'matrix'
      } else {
        typeLabel = 'tensor'
      }
    }
    return <Box component="span" whiteSpace="nowrap" fontStyle="italic">
        <Typography component="span">{`[${dimensions.join(', ')}] ${typeLabel}`}</Typography>
      </Box>
  } else {
    return <Box component="span" whiteSpace="nowarp">
      <Typography component="span">{String(value)}</Typography>
    </Box>
  }
}

function QuantityValue({value, def}) {
  return <Box
    marginTop={2} marginBottom={2} textAlign="center" fontWeight="bold"
  >
    <Typography>
      {String(value)}
    </Typography>
  </Box>
}

function Section({section, def}) {
  return <Content>
    <Compartment>
      <Definition def={def} />
    </Compartment>
    <Compartment title="sub sections">
      {def.sub_sections
        .filter(subSectionDef => section[subSectionDef.name])
        .map(subSectionDef => {
          const sectionDef = resolveRef(subSectionDef.sub_section)
          const key = subSectionDef.name
          return <Item key={key} itemKey={key}>
            <Typography component="span">
              <Box component="span" fontWeight="bold" component="span">
                {subSectionDef.name}
              </Box>
              {subSectionDef.repeats ? ` (${section[subSectionDef.name].length})` : ''}
            </Typography>
          </Item>
        })
      }
    </Compartment>
    <Compartment title="quantities">
      {def.quantities
        .filter(quantityDef => section[quantityDef.name])
        .map(quantityDef => {
          const key = quantityDef.name
          return <Item key={key} itemKey={key}>
            <Box component="span" whiteSpace="nowrap">
              <Typography component="span">
                <Box fontWeight="bold"  component="span">
                  {quantityDef.name}
                </Box>
              </Typography> = <QuantityItemPreview value={section[quantityDef.name]} def={quantityDef} />
            </Box>
          </Item>
        })
      }
    </Compartment>
  </Content>
}

function SubSection({sections, def}) {
  const [showAll, setShowAll] = useState(false)
  const length = sections.length

  const renderItem = (section, index) => (
    <Item key={index} itemKey={index}>
      <Typography><Box fontWeight="bold">{index}</Box></Typography>
    </Item>
  )

  if (length <= 5 || showAll) {
    return <Content>{sections.map(renderItem)}</Content>
  } else {
    return <Content>
      {sections.slice(0, 3).map(renderItem)}
      <Box marginLeft={3} marginTop={1} marginBottom={1}>
        <IconButton onClick={() => setShowAll(true)}>
          <MoreVertIcon />
        </IconButton>
      </Box>
      {sections.slice(length - 2, length).map((section, index) => renderItem(section, index + length - 2))}
    </Content>
  }
}

function Quantity({value, def}) {
  return <Content>
    <Definition def={def} />
    <QuantityValue value={value} def={def} />
  </Content>
}
