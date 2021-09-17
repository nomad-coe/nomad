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
import { withStyles, IconButton, Tooltip, Link } from '@material-ui/core'
import ClipboardIcon from '@material-ui/icons/Assignment'
import { CopyToClipboard } from 'react-copy-to-clipboard'

class DOIUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    doi: PropTypes.string.isRequired,
    parentheses: PropTypes.bool
  }

  static defaultProps = {
    parentheses: false
  }

  static styles = theme => ({
    root: {
      display: 'inline-flex',
      alignItems: 'center',
      flexDirection: 'row',
      flexWrap: 'nowrap'
    }
  })

  render() {
    const {classes, doi, parentheses} = this.props
    const url = `https://dx.doi.org/${doi}`
    return <span className={classes.root}>
      {parentheses && <div style={{marginRight: 0}}>(</div>}
      <Link href={url}>{doi}</Link>
      <CopyToClipboard
        text={url} onCopy={() => null}
      >
        <Tooltip title={`Copy DOI to clipboard`}>
          <IconButton style={{margin: 3, marginRight: 0, padding: 4}}>
            <ClipboardIcon style={{fontSize: 16}} />
          </IconButton>
        </Tooltip>
      </CopyToClipboard>
      {parentheses && <div style={{marginRight: 0}}>)</div>}
    </span>
  }
}

export const DOI = withStyles(DOIUnstyled.styles)(DOIUnstyled)
