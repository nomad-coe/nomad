import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Grid, Card, CardContent } from '@material-ui/core'
import PeriodicTable from '../search/PeriodicTable'
import QuantityHistogram from '../search/QuantityHistogram'
import { compose } from 'recompose'
import { withApi } from '../api'

class DFTSearchAggregations extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    quantities: PropTypes.object.isRequired,
    metric: PropTypes.string.isRequired,
    searchValues: PropTypes.object.isRequired,
    onChange: PropTypes.func.isRequired,
    info: PropTypes.object
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

  constructor(props) {
    super(props)
    this.handleExclusiveChanged = this.handleExclusiveChanged.bind(this)
  }

  state = {
    exclusive: false
  }

  handleExclusiveChanged() {
    const { searchValues } = this.props
    const value = !this.state.exclusive
    this.setState({exclusive: value})
    if (value) {
      searchValues.only_atoms = searchValues.only_atoms || searchValues.atoms
      delete searchValues.atoms
    } else {
      searchValues.atoms = searchValues.only_atoms || searchValues.atoms
      delete searchValues.only_atoms
    }
    this.props.onChange({searchValues: searchValues})
  }

  handleAtomsChanged(atoms) {
    if (this.state.exclusive) {
      this.setState({exclusive: false})
    }
    const searchValues = {...this.props.searchValues}
    searchValues.atoms = atoms
    if (searchValues.atoms.length === 0) {
      delete searchValues.atoms
    }
    delete searchValues.only_atoms
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
    const { classes, quantities, metric, searchValues, info } = this.props

    if (quantities.code_name && info) {
      // filter based on known codes, since elastic search might return 0 aggregations on
      // obsolete code names
      const filteredCodeNames = {}
      const defaultValue = {
        code_runs: 0
      }
      defaultValue[metric] = 0
      info.codes.forEach(key => {
        filteredCodeNames[key] = quantities.code_name[key] || defaultValue
      })
      quantities.code_name = filteredCodeNames
    }

    const quantity = (key, title, scale) => (<QuantityHistogram
      classes={{root: classes.quantity}} title={title || key} width={300}
      data={quantities[key]} metric={metric}
      value={searchValues[key]} defaultScale={scale}
      onChanged={(selection) => this.handleQuantityChanged(key, selection)}/>)

    return (
      <div className={classes.root}>
        <Card>
          <CardContent>
            <PeriodicTable
              aggregations={quantities.atoms} metric={metric}
              exclusive={this.state.exclusive}
              values={searchValues.atoms || searchValues.only_atoms || []}
              onChanged={(selection) => this.handleAtomsChanged(selection)}
              onExclusiveChanged={this.handleExclusiveChanged}
            />
          </CardContent>
        </Card>

        <Grid container spacing={24} className={classes.quantityGrid}>
          <Grid item xs={4}>
            {quantity('code_name', 'Code', 0.25)}
          </Grid>
          <Grid item xs={4}>
            {quantity('system', 'System type', 0.25)}
            {quantity('crystal_system', 'Crystal system', 1)}
          </Grid>
          <Grid item xs={4}>
            {quantity('basis_set', 'Basis set', 0.25)}
            {quantity('xc_functional', 'XC functionals', 0.5)}
          </Grid>
        </Grid>
      </div>
    )
  }
}

export default compose(withApi(false), withStyles(DFTSearchAggregations.styles))(DFTSearchAggregations)
