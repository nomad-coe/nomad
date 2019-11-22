import React from 'react'
import PropTypes from 'prop-types'
import { schema } from '../MetaInfoRepository'
import ValueCard from './ValueCard'
import Section from './Section'

export class Definition extends React.Component {
  static propTypes = {
    definition: PropTypes.object.isRequired,
    state: PropTypes.arrayOf(PropTypes.object).isRequired
  }

  render() {
    const {definition} = this.props
    if (schema.isSection(definition)) {
      return <Section {...this.props} section={definition}/>
    } else if (schema.isProperty(definition)) {
      return <ValueCard {...this.props}/>
    } else {
      return <div>{definition.mType}</div>
    }
  }
}

export class Definitions extends React.Component {
  static propTypes = {
    state: PropTypes.arrayOf(PropTypes.object).isRequired
  }
  render() {
    return (
      <div>
        {this.props.state.map(({definition, state}) => (
          <Definition key={definition.name} {...this.props} definition={definition} state={state}/>
        ))}
      </div>
    )
  }
}
