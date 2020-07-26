import React from 'react'
import Adaptor from './adaptors'
import { Item, Content } from './ArchiveBrowser'
import { Typography } from '@material-ui/core'

export function jsonAdaptorFactory(child) {
  if (Array.isArray(child)) {
    if (child.length === 1) {
      return new ObjectAdaptor(child[0])
    }
    return new ArrayAdaptor(child)
  } else if (typeof child === 'string') {
    return new ValueAdaptor(child)
  } else if (typeof child === 'object') {
    return new ObjectAdaptor(child)
  } else {
    return new ValueAdaptor(child)
  }
}

class ObjectAdaptor extends Adaptor {
  itemAdaptor(key) {
    return jsonAdaptorFactory(this.e[key])
  }
  render() {
    return <React.Fragment>
      {Object.keys(this.e).map(key => (
        <Item key={key} itemKey={key}>
          <Typography>
            {key}
          </Typography>
        </Item>
      ))}
    </React.Fragment>
  }
}

class ValueAdaptor extends Adaptor {
  render() {
    return <Content>
      <Typography>{String(this.e)}</Typography>
    </Content>
  }
}

class ArrayAdaptor extends Adaptor {
  itemAdaptor(index) {
    return jsonAdaptorFactory(this.e[index])
  }
  render() {
    return <React.Fragment>
      {this.e.map((_, index) => (
        <Item key={index} itemKey={index}>
          <Typography>
            {index + 1}
          </Typography>
        </Item>
      ))}
    </React.Fragment>
  }
}
