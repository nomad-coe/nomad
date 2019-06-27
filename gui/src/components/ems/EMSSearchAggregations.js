import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Grid, Card, CardContent } from '@material-ui/core'
import PeriodicTable from '../search/PeriodicTable'
import QuantityHistogram from '../search/QuantityHistogram'

class EMSSearchAggregations extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    aggregations: PropTypes.object.isRequired,
    metric: PropTypes.string.isRequired,
    searchValues: PropTypes.object.isRequired,
    onChange: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {},
    quantity: {
      marginTop: theme.spacing.unit * 2
    },
    quantityGrid: {
      marginBottom: theme.spacing.unit * 2
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
        <Card>
          <CardContent>
            <PeriodicTable
              aggregations={aggregations.atoms} metric={metric}
              values={searchValues.atoms || []}
              onChanged={(selection) => this.handleAtomsChanged(selection)}
            />
          </CardContent>
        </Card>

        <Grid container spacing={24} className={classes.quantityGrid}>
          <Grid item xs={6}>
            {quantity('method', 'Method')}
            {quantity('probing_method', 'Probing')}
          </Grid>
          <Grid item xs={6}>
            {quantity('sample_microstructure', 'Sample structure')}
            {quantity('sample_constituents', 'Sample constituents')}
          </Grid>
        </Grid>
      </div>
    )
  }
}

export default withStyles(EMSSearchAggregations.styles)(EMSSearchAggregations)
