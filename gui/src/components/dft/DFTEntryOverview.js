import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core'
import Quantity from '../Quantity'

class DFTEntryOverview extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    data: PropTypes.object.isRequired,
    loading: PropTypes.bool
  }

  static styles = theme => ({
    quantityColumn: {
      display: 'flex',
      flexDirection: 'column'
    },
    quantityRow: {
      display: 'flex',
      flexDirection: 'row',
      flexWrap: 'wrap',
      marginBottom: theme.spacing.unit
    }
  })

  render() {
    const { classes, data, loading } = this.props

    return (
      <div className={classes.quantityColumn}>
        <div className={classes.quantityRow}>
          <Quantity label='dft code' loading={loading}>
            {data.code_name}
          </Quantity>
          <Quantity label='dft code version' loading={loading}>
            {data.code_version}
          </Quantity>
        </div>
        <div className={classes.quantityRow}>
          <Quantity label='basis set' loading={loading}>
            {data.basis_set}
          </Quantity>
          <Quantity label='xc functional' loading={loading}>
            {data.xc_functional}
          </Quantity>
        </div>
        <div className={classes.quantityRow}>
          <Quantity label='system type' loading={loading}>
            {data.system}
          </Quantity>
          <Quantity label='crystal system' loading={loading}>
            {data.crystal_system}
          </Quantity>
          <Quantity label='spacegroup' loading={loading}>
            {data.spacegroup_symbol} ({data.spacegroup})
          </Quantity>
        </div>
      </div>
    )
  }
}

export default withStyles(DFTEntryOverview.styles)(DFTEntryOverview)
