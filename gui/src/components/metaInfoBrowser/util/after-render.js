import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core'

/**
 * Component that calls onMutation callback if something in its children's DOM changes. Uses
 * `window.MuationObserver`.
 */
class MutationObserver extends React.Component {
  static propTypes = {
    children: PropTypes.any,
    onMutation: PropTypes.func.isRequired
  }

  constructor(props) {
    super(props)
    this.ref = React.createRef()
    this.observer = new window.MutationObserver(mutations => {
      mutations.forEach(mutation => this.onMutation(mutation))
    })
  }

  onMutation(mutation) {
    this.props.onMutation(mutation)
  }

  componentDidMount() {
    const config = { attributes: true, childList: true, characterData: true, subtree: true }
    this.observer.observe(this.ref.current, config)
  }

  componentWillUnmount() {
    this.observer.disconnect()
  }

  render() {
    return <div ref={this.ref}>{this.props.children}</div>
  }
}

const AfterRenderContext = React.createContext(() => null)

/**
 * Component that measures and transfers measures of child components to an AfterRender component.
 * It uses the AfterRenderContext provided by the AfterRender component to communicate
 * with the AfterRender instance.
 */
export class AfterRenderMeasure extends React.Component {
  static propTypes = {
    children: PropTypes.any,
    measureId: PropTypes.oneOfType([PropTypes.number, PropTypes.string]).isRequired,
    measureData: PropTypes.object
  }

  constructor(props) {
    super(props)
    this.measureRef = React.createRef()
  }

  render() {
    return (
      <AfterRenderContext.Consumer>
        {(measureRefCallback) => {
          measureRefCallback(this.measureRef, this.props.measureId, this.props.measureData)
          return (
            <div ref={this.measureRef}>
              {this.props.children}
            </div>
          )
        }}
      </AfterRenderContext.Consumer>
    )
  }
}

class AfterRenderUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    children: PropTypes.any,
    afterRender: PropTypes.func.isRequired
  }

  static styles = {
    container: {},
    content: {},
    afterRender: {}
  }

  constructor(props) {
    super(props)
    this.state = {afterRender: false, needsRender: true}
    this.containerRef = React.createRef()
    this.measureRefs = {}
    this.measureRefCallback = this.measureRefCallback.bind(this)
    this.onMutation = this.onMutation.bind(this)
  }

  update() {
    if (this.state.needsRender) {
      this.setState({afterRender: true, needsRender: false})
    }
  }

  componentDidMount() {
    this.update()
  }

  componentDidUpdate() {
    this.update()
  }

  onMutation() {
    if (!this.state.needsRender) {
      this.setState({needsRender: true, afterRender: false})
    }
  }

  measureRefCallback(ref, key, data) {
    this.measureRefs[key] = { ref: ref, key: key, data: data }
  }

  afterRender() {
    // a function to calculate the bounding box relative to the whole document
    const getMeasure = (elem) => {
      const box = elem.getBoundingClientRect()

      const body = document.body
      var rootNode = document.documentElement

      const scrollTop = window.pageYOffset || rootNode.scrollTop || body.scrollTop
      const scrollLeft = window.pageXOffset || rootNode.scrollLeft || body.scrollLeft

      const clientTop = rootNode.clientTop || body.clientTop || 0
      const clientLeft = rootNode.clientLeft || body.clientLeft || 0

      const top = box.top + scrollTop - clientTop
      const left = box.left + scrollLeft - clientLeft

      return {
        top: Math.round(top),
        left: Math.round(left),
        right: Math.round(left + box.width),
        bottom: Math.round(top + box.height)
      }
    }

    // a function to calculate the bounding box relative to the container
    const containerMeasure = getMeasure(this.containerRef.current)
    const getRelativeMeasure = (elem) => {
      const elemMeasure = getMeasure(elem)
      return {
        top: elemMeasure.top - containerMeasure.top,
        left: elemMeasure.left - containerMeasure.left,
        right: elemMeasure.right - containerMeasure.left,
        bottom: elemMeasure.bottom - containerMeasure.top
      }
    }

    const relativeMeasures = {}
    Object.keys(this.measureRefs).filter(measureKey => this.measureRefs[measureKey].ref.current).forEach(measureKey => {
      relativeMeasures[measureKey] = {
        ...this.measureRefs[measureKey],
        box: getRelativeMeasure(this.measureRefs[measureKey].ref.current)
      }
    })
    return this.props.afterRender(relativeMeasures)
  }

  render() {
    const {classes} = this.props
    return (
      <AfterRenderContext.Provider value={this.measureRefCallback}>
        <div ref={this.containerRef} className={classes.container}>
          <MutationObserver onMutation={this.onMutation}>
            <div className={classes.content}>
              {this.props.children}
            </div>
          </MutationObserver>
          {this.state.afterRender
            ? <div className={classes.afterRender}>
              {this.afterRender()}
            </div>
            : ''}
        </div>
      </AfterRenderContext.Provider>
    )
  }
}

export const AfterRender = withStyles(AfterRenderUnstyled.styles)(AfterRenderUnstyled)

class Button extends React.Component {
  constructor(props) {
    super(props)
    this.state = {a: 0}
  }

  render() {
    return <div onClick={() => { this.setState({a: 1}) }}>{this.state.a}</div>
  }
}

const afterRenderTestStyles = {
  afterRender: {
    color: 'red'
  }
}

class AfterRenderTestUnstyled extends React.Component {
  constructor(props) {
    super(props)
    this.measures = {}
    this.afterRender = this.afterRender.bind(this)
  }

  afterRender(measures) {
    console.log(measures)
    return ('..empty..' + measures[1])
  }

  render() {
    const {classes} = this.props
    return (
      <div style={{padding: 10}}>
        <AfterRender afterRender={this.afterRender} classes={{afterRender: classes.afterRender}}>
          <AfterRenderMeasure measureId={1}>Hello</AfterRenderMeasure>
          <AfterRenderMeasure measureId={2}>Hello World</AfterRenderMeasure>
          <Button>click</Button>
        </AfterRender>
      </div>
    )
  }
}

AfterRenderTestUnstyled.propTypes = {
  classes: PropTypes.object.isRequired
}

export const AfterRenderTest = withStyles(afterRenderTestStyles)(AfterRenderTestUnstyled)
