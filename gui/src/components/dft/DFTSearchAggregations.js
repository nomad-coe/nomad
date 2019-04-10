import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Grid } from '@material-ui/core'
import PeriodicTable from '../search/PeriodicTable'
import QuantityHistogram from '../search/QuantityHistogram'

class DFTSearchAggregations extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    aggregations: PropTypes.object.isRequired,
    metric: PropTypes.string.isRequired,
    searchValues: PropTypes.object.isRequired,
    onChange: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {},
    quantityGrid: {
      minWidth: 524,
      maxWidth: 924,
      margin: 'auto',
      width: '100%'
    },
    quantity: {
      marginTop: theme.spacing.unit * 2
    }
  })

  handleAtomsChanged(atoms) {
    const searchValues = {...this.props.searchValues}
    searchValues.atoms = atoms
    if (searchValues.atoms.length === 0) {
      delete searchValues.atoms
    }
    this.props.onChange({searchValues: searchValues})
  }

  handleQuantityChanged(quantity, selection) {
    const searchValues = {...this.props.searchValues}
    if (selection) {
      searchValues[quantity] = selection
    } else {
      delete searchValues[quantity]
    }
    this.props.onChange({searchValues: searchValues})
  }

  render() {
    const { classes, aggregations, metric, searchValues } = this.props

    const quantity = (key, title) => (<QuantityHistogram
      classes={{root: classes.quantity}} title={title || key} width={300}
      data={aggregations[key]} metric={metric}
      value={searchValues[key]}
      onChanged={(selection) => this.handleQuantityChanged(key, selection)}/>)

    return (
      <div className={classes.root}>
        <PeriodicTable
          aggregations={aggregations.atoms} metric={metric}
          values={searchValues.atoms || []}
          onChanged={(selection) => this.handleAtomsChanged(selection)}
        />

        <Grid container spacing={24} className={classes.quantityGrid}>

          <Grid item xs={4}>
            {quantity('system', 'System')}
            {quantity('crystal_system', 'Crystal system')}
          </Grid>
          <Grid item xs={4}>
            {quantity('basis_set', 'Basis set')}
            {quantity('xc_functional', 'XC functionals')}
          </Grid>
          <Grid item xs={4}>
            {quantity('code_name', 'Code')}
          </Grid>
        </Grid>
      </div>
    )
  }
}

export default withStyles(DFTSearchAggregations.styles)(DFTSearchAggregations)
