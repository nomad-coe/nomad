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
