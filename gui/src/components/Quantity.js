import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Typography } from '@material-ui/core'

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
    quantity: PropTypes.string
  }

  static styles = theme => ({
    root: {},
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
    const {classes, children, label, typography, loading, placeholder, noWrap, row, column, quantity, data} = this.props
    let content = null
    if (!loading) {
      if (!(data && quantity && !data[quantity])) {
        if (!children || children.length === 0) {
          const value = data && quantity ? data[quantity] : null
          if (value) {
            content = <Typography noWrap={noWrap} variant={typography}>
              {value}
            </Typography>
          } else {
            content = <Typography noWrap={noWrap} variant={typography}>
              <i>{placeholder || 'unavailable'}</i>
            </Typography>
          }
        } else {
          content = children
        }
      } else {
        content = <Typography noWrap={noWrap} variant={typography}>
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
          {loading ? <Typography noWrap={noWrap} variant={typography}>
            <i>loading ...</i>
          </Typography> : content}
        </div>
      )
    }
  }
}

export default withStyles(Quantity.styles)(Quantity)
