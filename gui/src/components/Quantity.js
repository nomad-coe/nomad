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
    noWrap: PropTypes.bool
  }

  static styles = theme => ({
    root: {
      margin: '8px 24px 0px 0'
    }
  })

  render() {
    const {classes, children, label, typography, loading, placeholder, noWrap} = this.props
    const content = (!children || children.length === 0) ? null : children

    return (
      <div className={classes.root}>
        <Typography variant="caption">{label}</Typography>
        <Typography noWrap={noWrap} variant={typography || 'body1'}>{content || <i>{loading ? 'loading...' : placeholder || 'unavailable'}</i>}</Typography>
      </div>
    )
  }
}

export default withStyles(Quantity.styles)(Quantity)
