import React from 'react'
import PropTypes from 'prop-types'

/**
 * This is a kinda-HOC that allows to keep a component alive while not being visible.
 */
export default class KeepState extends React.Component {
  static propTypes = {
    visible: PropTypes.bool,
    render: PropTypes.func.isRequired
  }

  state = {
    props: null
  }

  update() {
    const { visible, render, ...props } = this.props
    if (this.props.visible) {
      this.setState({props: props})
    }
  }

  componentDidMount() {
    this.update()
  }

  componentDidUpdate(prevProps) {
    if (prevProps.visible !== this.props.visible && this.props.visible) {
      this.update()
    }
  }

  render() {
    const { visible, render, ...other } = this.props
    const props = visible ? other : this.state.props

    if (props) {
      return <div style={{display: visible ? 'block' : 'none'}}>
        {render(props)}
      </div>
    } else {
      return ''
    }
  }
}