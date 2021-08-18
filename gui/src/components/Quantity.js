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
import { withStyles, Typography, Tooltip, IconButton } from '@material-ui/core'
import ClipboardIcon from '@material-ui/icons/Assignment'
import { CopyToClipboard } from 'react-copy-to-clipboard'
import _ from 'lodash'
import searchQuantities from '../searchQuantities'

class Quantity extends React.Component {
  static propTypes = {
    classes: PropTypes.object,
    children: PropTypes.node,
    label: PropTypes.string,
    typography: PropTypes.string,
    loading: PropTypes.bool,
    placeholder: PropTypes.string,
    noWrap: PropTypes.bool,
    row: PropTypes.bool,
    column: PropTypes.bool,
    flex: PropTypes.bool,
    data: PropTypes.object,
    quantity: PropTypes.oneOfType([
      PropTypes.string,
      PropTypes.func
    ]),
    withClipboard: PropTypes.bool,
    ellipsisFront: PropTypes.bool,
    hideIfUnavailable: PropTypes.bool,
    description: PropTypes.string
  }

  static styles = theme => ({
    root: {
      maxWidth: '100%'
    },
    valueContainer: {
      display: 'flex',
      alignItems: 'center',
      flexDirection: 'row',
      maxWidth: '100%'
    },
    value: {
      flexGrow: 1
    },
    ellipsis: {
      whiteSpace: 'nowrap',
      overflow: 'hidden',
      textOverflow: 'ellipsis'
    },
    ellipsisFront: {
      direction: 'rtl',
      textAlign: 'left'
    },
    valueAction: {},
    valueActionButton: {
      padding: 4
    },
    valueActionIcon: {
      fontSize: 16
    },
    row: {
      display: 'flex',
      flexWrap: 'wrap',
      flexDirection: 'row',
      '& > :not(:last-child)': {
        marginRight: theme.spacing(3)
      }
    },
    column: {
      display: 'flex',
      flexDirection: 'column',
      '& > :not(:first-child)': {
        marginTop: theme.spacing(1)
      }
    },
    flex: {
      display: 'flex',
      flexDirection: 'row',
      flexWrap: 'wrap',
      alignContent: 'flex-start',
      '& div': {
        marginRight: theme.spacing(1)
      }
    },
    label: {
      color: 'rgba(0, 0, 0, 0.54)',
      fontSize: '0.75rem',
      fontWeight: 500
    },
    quantityList: {
      display: 'flex',
      flexDirection: 'column'
    }
  })

  render() {
    const {
      classes, children, label, typography, loading, placeholder, noWrap, row, column, flex,
      quantity, data, withClipboard, ellipsisFront, hideIfUnavailable, description
    } = this.props
    let content = null
    let clipboardContent = null

    let valueClassName = classes.value
    if (noWrap && ellipsisFront) {
      valueClassName = `${valueClassName} ${classes.ellipsisFront}`
    }

    let value
    if (!loading) {
      if (typeof quantity === 'string') {
        value = data && quantity && _.get(data, quantity)
      } else if (children) {
      } else {
        try {
          value = quantity(data)
        } catch {
          value = undefined
        }
      }

      if (value === 'not processed') {
        value = 'unavailable'
      }

      if (value === 'unavailable') {
        value = ''
      }

      if ((!value && !children) && hideIfUnavailable) {
        return ''
      }

      if (children && children.length !== 0) {
        content = children
      } else if (value || value === 0) {
        if (Array.isArray(value)) {
          value = value.join(', ')
        }
        clipboardContent = value
        content = <Typography noWrap={noWrap} variant={typography} className={valueClassName}>
          {value}
        </Typography>
      } else {
        content = <Typography noWrap={noWrap} variant={typography} className={valueClassName}>
          <i>{placeholder || 'unavailable'}</i>
        </Typography>
      }
    }

    const useLabel = label || (typeof quantity === 'string' ? quantity : 'MISSING LABEL')

    if (row || column || flex) {
      return <div className={row ? classes.row : (column ? classes.column : classes.flex)}>{children}</div>
    } else {
      return (
        <Tooltip title={description || (searchQuantities[quantity] && searchQuantities[quantity].description) || ''}>
          <div className={classes.root}>
            <Typography noWrap classes={{root: classes.label}} variant="caption">{useLabel}</Typography>
            <div className={classes.valueContainer}>
              {loading
                ? <Typography noWrap={noWrap} variant={typography} className={valueClassName}>
                  <i>loading ...</i>
                </Typography> : content}
              {withClipboard
                ? <CopyToClipboard className={classes.valueAction} text={clipboardContent} onCopy={() => null}>
                  <Tooltip title={`Copy ${useLabel} to clipboard`}>
                    <div>
                      <IconButton disabled={!clipboardContent} classes={{root: classes.valueActionButton}} >
                        <ClipboardIcon classes={{root: classes.valueActionIcon}}/>
                      </IconButton>
                    </div>
                  </Tooltip>
                </CopyToClipboard> : ''}
            </div>
          </div>
        </Tooltip>
      )
    }
  }
}

export default withStyles(Quantity.styles)(Quantity)
