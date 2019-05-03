import React from 'react'
import PropTypes from 'prop-types'
import periodicTableData from './PeriodicTableData'
import { withStyles, Typography, Button, Tooltip } from '@material-ui/core'
import chroma from 'chroma-js'

const elements = []
for (var i = 0; i < 10; i++) {
  elements[i] = Array.apply(null, Array(18))
}
periodicTableData.elements.forEach(element => {
  elements[element.ypos - 1][element.xpos - 1] = element
  element.category = element.category.replace(' ', '')
})

class ElementUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    element: PropTypes.object.isRequired,
    onClick: PropTypes.func,
    selected: PropTypes.bool,
    count: PropTypes.number.isRequired,
    heatmapScale: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {
      position: 'relative'
    },
    button: {
      border: '1px solid',
      borderColor: '#555',
      paddingTop: theme.spacing.unit,
      paddingBottom: theme.spacing.unit,
      paddingLeft: 0,
      paddingRight: 0,
      width: '100%',
      textAlign: 'center',
      fontSize: '1rem',
      fontWeight: 700,
      textTransform: 'none',
      minWidth: 0,
      minHeight: 0,
      borderRadius: 0
    },
    containedPrimary: {
      backgroundColor: theme.palette.primary.main
    },
    number: {
      position: 'absolute',
      top: 2,
      left: 2,
      margin: 0,
      padding: 0,
      fontSize: 8
    },
    count: {
      position: 'absolute',
      bottom: 2,
      right: 2,
      margin: 0,
      padding: 0,
      fontSize: 8
    }
  })

  render() {
    const {classes, element, selected, count, heatmapScale} = this.props
    const buttonClasses = {
      root: classes.button,
      containedPrimary: classes.containedPrimary
    }
    const disabled = count <= 0

    return (
      <div className={classes.root}>
        <Tooltip title={element.name}>
          <div>
            <Button
              disabled={disabled}
              classes={buttonClasses}
              style={{backgroundColor: count > 0 && !selected ? heatmapScale(count).hex() : undefined}}
              onClick={this.props.onClick} variant="contained"
              color={selected ? 'primary' : 'default'}
            >
              {element.symbol}
            </Button>
          </div>
        </Tooltip>
        <Typography
          classes={{root: classes.number}} variant="caption"
          style={selected ? {color: 'white'} : disabled ? {color: '#BDBDBD'} : {}}>
          {element.number}
        </Typography>
        {count >= 0
          ? <Typography
            classes={{root: classes.count}} variant="caption"
            style={selected ? {color: 'white'} : disabled ? {color: '#BDBDBD'} : {}}>
            {count}
          </Typography> : ''
        }
      </div>
    )
  }
}

const Element = withStyles(ElementUnstyled.styles)(ElementUnstyled)

class PeriodicTable extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    aggregations: PropTypes.object,
    metric: PropTypes.string.isRequired,
    values: PropTypes.array.isRequired,
    onChanged: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {
      overflowX: 'scroll'
    },
    table: {
      margin: 'auto',
      width: '100%',
      minWidth: 500,
      maxWidth: 900,
      tableLayout: 'fixed',
      borderSpacing: theme.spacing.unit * 0.5
    }
  })

  onElementClicked(element) {
    const index = this.props.values.indexOf(element)
    const isClicked = index >= 0
    let selected
    if (isClicked) {
      selected = [...this.props.values]
      selected.splice(index, 1)
    } else {
      selected = [element, ...this.props.values]
    }

    this.props.onChanged(selected)
  }

  unSelectedAggregations() {
    const { aggregations, metric, values } = this.props
    return Object.keys(aggregations)
      .filter(key => values.indexOf(key) === -1)
      .map(key => aggregations[key][metric])
  }

  render() {
    const {classes, aggregations, metric, values} = this.props
    const max = aggregations ? Math.max(...this.unSelectedAggregations()) || 1 : 1
    const heatmapScale = chroma.scale(['#ffcdd2', '#d50000']).domain([1, max], 10, 'log')
    return (
      <div className={classes.root}>
        <table className={classes.table}>
          <tbody>
            {elements.map((row, i) => (
              <tr key={i}>
                {row.map((element, j) => (
                  <td key={j}>
                    {element
                      ? <Element
                        element={element}
                        count={aggregations ? (aggregations[element.symbol] || {})[metric] || 0 : 0}
                        heatmapScale={heatmapScale}
                        relativeCount={aggregations ? ((aggregations[element.symbol] || {})[metric] || 0) / max : 0}
                        onClick={() => this.onElementClicked(element.symbol)}
                        selected={values.indexOf(element.symbol) >= 0}
                      /> : ''}
                  </td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    )
  }
}

export default withStyles(PeriodicTable.styles)(PeriodicTable)
