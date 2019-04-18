import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import grey from '@material-ui/core/colors/grey'
import blue from '@material-ui/core/colors/blue'
import { AfterRender } from './util/after-render'
import { Definitions } from './Definition'
import { updateList } from './util/data'
import { schema } from '../MetaInfoRepository'

class ViewerUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    rootElement: PropTypes.object.isRequired,
    packages: PropTypes.arrayOf(PropTypes.object).isRequired,
    visiblePackages: PropTypes.arrayOf(PropTypes.string).isRequired
  }

  static styles = (theme) => ({
    root: {

    },
    packages: {
      padding: theme.spacing.unit
    },
    container: {
      position: 'relative'
    },
    canvas: {
      display: 'flex',
      flexDirection: 'row',
      flexWrap: 'nowrap',
      justifyContent: 'flex-start',
      alignItems: 'center',
      padding: theme.spacing.unit,
      zIndex: 1
    },
    sankey: {
      position: 'absolute',
      left: 0,
      top: 0,
      pointerEvents: 'none',
      zIndex: 0
    },
    references: {
      fill: blue[200],
      fillOpacity: 0.8
    },
    containments: {
      fill: grey[200],
      fillOpacity: 0.8
    }
  })

  constructor(props) {
    super(props)
    this.canvasRef = React.createRef()
    this.renderSankey = this.renderSankey.bind(this)
    this.toggleDefinition = this.toggleDefinition.bind(this)
    this.isVisible = this.isVisible.bind(this)

    this.state = {
      packages: this.props.packages.map(pkg => ({
        name: pkg.name,
        visible: this.props.visiblePackages.find(name => pkg.name === name) !== undefined
      })),
      definitions: [{
        definition: this.props.rootElement,
        state: []
      }]
    }
  }

  isVisible(definition) {
    let state = { state: this.state.definitions }
    while (state && state.definition !== definition) {
      state = state.state.find(state => state.definition === definition || schema.isParent(definition, state.definition))
    }
    return state
  }

  toggleDefinition(toggledDefinition, toggledVisible) {
    const isParentOfCurrent = schema.isParent(this.props.rootElement, toggledDefinition)
    if (isParentOfCurrent && !toggledVisible) {
      // a special case

      // determine the overall root of the toggledDefinition/rootElement (its the same)
      const rootIndex = this.state.definitions.findIndex(state => schema.isParent(this.props.rootElement, state.definition))
      const rootState = this.state.definitions[rootIndex]

      // go down from the root until, we find find the toggledDefinition
      let newRootState = rootState
      while (newRootState.definition !== toggledDefinition) {
        newRootState = newRootState.state.find(state => state.definition === toggledDefinition || schema.isParent(toggledDefinition, state.definition))
      }

      // find the child that is rootElement or a parent of it
      newRootState = newRootState.state.find(state => state.definition === this.props.rootElement || schema.isParent(this.props.rootElement, state.definition))

      // replace the old with the new root
      this.setState({ definitions: updateList(this.state.definitions, rootIndex, () => newRootState) })
      return
    }

    let state = this.state.definitions
    // traverse the state to create a copy, find its parent (or itself), and possibly add it
    const handle = (state) => {
      const index = state.findIndex(({definition}) => definition === toggledDefinition)
      if (toggledVisible) {
        if (index === -1) {
          return [...state, {definition: toggledDefinition, state: []}]
        }
      } else {
        if (index !== -1) {
          return [...state.slice(0, index), ...state.slice(index + 1)]
        }
      }
      return state
    }

    let wasHandled = false
    const clone = (definition, state) => {
      if (definition === toggledDefinition.parent) {
        state = handle(state)
        wasHandled = true
        return state
      } else {
        return state.map(child => ({
          ...child,
          state: clone(child.definition, child.state)
        }))
      }
    }

    state = clone(null, state)
    if (!wasHandled) {
      const childIndex = state.findIndex(({definition}) => definition.parent === toggledDefinition)
      if (childIndex === -1) {
        // add, remove it to the top-level if parent not found in the tree
        state = handle(state)
      } else {
        // replace the root definition with the new top-most parent
        state = updateList(state, childIndex, (oldRoot) => ({
          definition: toggledDefinition,
          state: [oldRoot]
        }))
      }
    }

    this.setState({definitions: state})
  }

  renderSankey(measures) {
    // create connections from sections back to their features
    let containments = []
    let references = []
    Object.values(measures).forEach(({key, box, data}) => {
      const splitKey = key.split('.')
      const mType = splitKey[0]
      const name = splitKey[1]
      if (mType === 'section' || mType === 'property') {
        let featureMeasure = measures[`feature.${name}`]
        if (!featureMeasure) {
          if (data.parent) {
            if (mType === 'property') {
              featureMeasure = measures[`properties.${data.parent.name}`]
            } else if (mType === 'section') {
              featureMeasure = measures[`sections.${data.parent.name}`]
            }
          }
        }
        if (featureMeasure) {
          containments.push({
            section: box,
            feature: featureMeasure.box
          })
        }
      }
      if (data && mType === 'property' && data.mType === 'reference' && data.referencedSection) {
        const referencedMeasure = measures[`section.${data.referencedSection.name}`]
        if (referencedMeasure) {
          references.push({
            referencedSection: referencedMeasure.box,
            reference: box
          })
        }
      }
    })

    // create polygons
    containments = containments.map(({feature, section}) => ([
      feature.right, feature.top,
      feature.right, feature.bottom,
      section.left, section.bottom,
      section.left, section.top
    ]))
    references = references.map(({referencedSection, reference}) => ([
      referencedSection.right, referencedSection.top,
      referencedSection.right, referencedSection.bottom,
      reference.left, reference.bottom,
      reference.left, reference.top
    ]))

    // render the connections
    const {classes} = this.props
    return (
      <svg width={this.canvasRef.current.offsetWidth} height={this.canvasRef.current.offsetHeight}>
        <g className={classes.references}>
          {references.map((points, index) => <polygon key={index} points={points}/>)}
        </g>
        <g className={classes.containments}>
          {containments.map((points, index) => <polygon key={index} points={points}/>)}
        </g>
      </svg>
    )
  }

  render() {
    const {classes} = this.props
    return (
      <div className={classes.root}>
        {/* <div className={classes.packages}>
          {this.state.packages.map(({name, visible}, index) => (
            <CheckChip key={index} label={name} checked={visible}/>
          ))}
        </div> */}
        <div ref={this.canvasRef} className={classes.container}>
          <AfterRender classes={{content: classes.canvas, afterRender: classes.sankey}} afterRender={this.renderSankey}>
            <Definitions current={this.props.rootElement}
              state={this.state.definitions}
              isVisible={this.isVisible}
              toggleDefinition={this.toggleDefinition}/>
          </AfterRender>
        </div>
      </div>
    )
  }
}

export default withStyles(ViewerUnstyled.styles)(ViewerUnstyled)
