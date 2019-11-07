import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Grid, Card, CardContent } from '@material-ui/core'
import PeriodicTable from '../search/PeriodicTable'
import QuantityHistogram from '../search/QuantityHistogram'
import { compose } from 'recompose'
import { withApi } from '../api'
import SearchContext from '../search/SearchContext'


class QuantityUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    quantity: PropTypes.string.isRequired,
    metric: PropTypes.string.isRequired,
    title: PropTypes.string,
    scale: PropTypes.number
  }
  static styles = theme => ({
    root: {
      marginTop: theme.spacing.unit * 2
    }
  })

  static contextType = SearchContext.type

  render() {
    const {classes, scale, quantity, title, ...props} = this.props
    const {state: {response, query}, setQuery} = this.context

    return <QuantityHistogram
      classes={{root: classes.root}}
      width={300}
      defaultScale={scale || 1}
      title={title || quantity}
      data={response.statistics[quantity]}
      value={query[quantity]}
      onChanged={selection => setQuery({...query, [quantity]: selection})}
      {...props} />
  }
}

const Quantity = withStyles(QuantityUnstyled.styles)(QuantityUnstyled)


class DFTSearchAggregations extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    metric: PropTypes.string.isRequired,
    info: PropTypes.object
  }

  static styles = theme => ({
    root: {},
    quantityGrid: {
      marginBottom: theme.spacing.unit * 2
    }
  })

  constructor(props) {
    super(props)
    this.handleExclusiveChanged = this.handleExclusiveChanged.bind(this)
    this.handleAtomsChanged = this.handleAtomsChanged.bind(this)
  }

  state = {
    exclusive: false
  }

  handleExclusiveChanged() {
    this.setState({exclusive: !this.state.exclusive}, () => {
      const {state: {query}, setQuery} = this.context
      if (this.state.exclusive) {
        setQuery({...query, only_atoms: query.atoms, atoms: []})
      } else {
        setQuery({...query, atoms: query.only_atoms, only_atoms: []})
      }
    })
  }

  handleAtomsChanged(atoms) {
    if (this.state.exclusive) {
      this.setState({exclusive: false})
    }

    const {state: {query}, setQuery} = this.context
    setQuery({...query, atoms: atoms, only_atoms: []})
  }

  componentDidMount() {
    const {state: {query}, setQuery} = this.context
    setQuery({...query, atoms: [], only_atoms: []})
  }

  static contextType = SearchContext.type

  render() {
    const {classes, info, metric} = this.props
    const {state: {response: {statistics}, query}} = this.context

    if (statistics.code_name && info) {
      // filter based on known codes, since elastic search might return 0 aggregations on
      // obsolete code names
      const filteredCodeNames = {}
      const defaultValue = {
        code_runs: 0
      }
      defaultValue[metric] = 0
      info.codes.forEach(key => {
        filteredCodeNames[key] = statistics.code_name[key] || defaultValue
      })
      statistics.code_name = filteredCodeNames
    }

    return (
      <div className={classes.root}>
        <Card>
          <CardContent>
            <PeriodicTable
              aggregations={statistics.atoms}
              metric={metric}
              exclusive={this.state.exclusive}
              values={[...(query.atoms || []), ...(query.only_atoms || [])]}
              onChanged={this.handleAtomsChanged}
              onExclusiveChanged={this.handleExclusiveChanged}
            />
          </CardContent>
        </Card>

        <Grid container spacing={24} className={classes.quantityGrid}>
          <Grid item xs={4}>
            <Quantity quantity="code_name" title="Code" scale={0.25} metric={metric} />
          </Grid>
          <Grid item xs={4}>
            <Quantity quantity="system" title="System type" scale={0.25} metric={metric} />
            <Quantity quantity="crystal_system" title="Crystal system" scale={1} metric={metric} />
          </Grid>
          <Grid item xs={4}>
            <Quantity quantity="basis_set" title="Basis set" scale={0.25} metric={metric} />
            <Quantity quantity="xc_functional" title="XC functionals" scale={0.5} metric={metric} />
          </Grid>
        </Grid>
      </div>
    )
  }
}

export default compose(withApi(false), withStyles(DFTSearchAggregations.styles))(DFTSearchAggregations)
