import React from 'react'
import PropTypes from 'prop-types'
import { withRouter } from 'react-router-dom'

/**
 * This is a kinda-HOC that allows to keep a component alive while not being visible.
 */
class KeepState extends React.Component {
  static propTypes = {
    visible: PropTypes.bool,
    render: PropTypes.func.isRequired
  }

  state = {
    props: null,
    update: 0
  }

  update(becameVisible) {
    const { visible, render, ...props } = this.props
    if (visible) {
      this.setState({props: props, update: this.state.update + (this.state.props ? 1 : 0)})
    }
  }

  componentDidMount() {
    if (this.props.visible) {
      this.update()
    }
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
      props.update = this.state.update
      return <div style={{display: visible ? 'block' : 'none'}}>
        {render(props)}
      </div>
    } else {
      return ''
    }
  }
}

const KeepStateWithRouter = withRouter(KeepState)
export default KeepStateWithRouter
