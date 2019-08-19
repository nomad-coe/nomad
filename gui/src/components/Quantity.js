import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Typography, Tooltip, IconButton } from '@material-ui/core'
import ClipboardIcon from '@material-ui/icons/Assignment'
import { CopyToClipboard } from 'react-copy-to-clipboard'

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
    data: PropTypes.object,
    quantity: PropTypes.string,
    withClipboard: PropTypes.bool
  }

  static styles = theme => ({
    root: {},
    valueContainer: {
      display: 'flex',
      alignItems: 'center',
      flexDirection: 'row'
    },
    value: {
      flexGrow: 1
    },
    valueAction: {},
    valueActionButton: {
      padding: 2,
      paddingButtom: 3
    },
    valueActionIcon: {
      fontSize: 16
    },
    row: {
      display: 'flex',
      flexDirection: 'row',
      '& > :not(:first-child)': {
        marginLeft: theme.spacing.unit * 3
      }
    },
    column: {
      display: 'flex',
      flexDirection: 'column',
      '& > :not(:first-child)': {
        marginTop: theme.spacing.unit * 1
      }
    }
  })

  render() {
    const {classes, children, label, typography, loading, placeholder, noWrap, row, column, quantity, data, withClipboard} = this.props
    let content = null
    let clipboardContent = null
    if (!loading) {
      if (!(data && quantity && !data[quantity])) {
        if (!children || children.length === 0) {
          const value = data && quantity ? data[quantity] : null
          if (value) {
            clipboardContent = value
            content = <Typography noWrap={noWrap} variant={typography} className={classes.value}>
              {value}
            </Typography>
          } else {
            content = <Typography noWrap={noWrap} variant={typography} className={classes.value}>
              <i>{placeholder || 'unavailable'}</i>
            </Typography>
          }
        } else {
          content = children
        }
      } else {
        content = <Typography noWrap={noWrap} variant={typography} className={classes.value}>
          <i>{placeholder || 'unavailable'}</i>
        </Typography>
      }
    }

    if (row || column) {
      return <div className={row ? classes.row : classes.column}>{children}</div>
    } else {
      return (
        <div className={classes.root}>
          <Typography noWrap variant="caption">{label || quantity}</Typography>
          <div className={classes.valueContainer}>
            {loading
              ? <Typography noWrap={noWrap} variant={typography} className={classes.value}>
                <i>loading ...</i>
              </Typography> : content}
            {withClipboard
              ? <CopyToClipboard className={classes.valueAction} text={clipboardContent} onCopy={() => null}>
                <Tooltip title={`Copy ${label || quantity} to clipboard`}>
                  <div>
                    <IconButton disabled={!clipboardContent} classes={{root: classes.valueActionButton}} >
                      <ClipboardIcon classes={{root: classes.valueActionIcon}}/>
                    </IconButton>
                  </div>
                </Tooltip>
              </CopyToClipboard> : ''}
          </div>
        </div>
      )
    }
  }
}

export default withStyles(Quantity.styles)(Quantity)
