import React from 'react'
import PropTypes from 'prop-types'
import periodicTableData from './PeriodicTableData'
import { withStyles, Paper, Typography, Button, Tooltip } from '@material-ui/core'

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
    selected: PropTypes.bool
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
    contained: {
      backgroundColor: 'white'
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
    actinide: {
      backgroundColor: '#E8EAF6'
    },
    alkalimetal: {
      backgroundColor: '#E3F2FD'
    },
    alkalineearthmetal: {
      backgroundColor: '#EDE7F6'
    },
    diatomicnonmetal: {
      backgroundColor: '#F3E5F5'
    },
    lanthanide: {
      backgroundColor: '#FCE4EC'
    },
    metalloid: {
      backgroundColor: '#FFEBEE'
    },
    noblegas: {
      backgroundColor: '#E0F7FA'
    },
    polyatomicnonmetal: {
      backgroundColor: '#E0F2F1'
    },
    'post-transitionmetal': {
      backgroundColor: '#E1F5FE'
    },
    transitionmetal: {
      backgroundColor: '#F9FBE7'
    }
  })

  render() {
    const {classes, element, selected} = this.props
    const buttonClasses = {
      root: classes.button,
      contained: (!selected ? classes[element.category] : null) || classes.contained,
      containedPrimary: classes.containedPrimary
    }

    return (
      <div className={classes.root}>
        <Tooltip title={element.name}>
          <Button
            classes={buttonClasses}
            onClick={this.props.onClick} variant="contained"
            color={selected ? 'primary' : 'default'}
          >
            {element.symbol}
          </Button>
        </Tooltip>
        <Typography
          classes={{root: classes.number}} variant="caption"
          style={selected ? {color: 'white'} : {}}>
          {element.number}
        </Typography>
      </div>
    )
  }
}

const Element = withStyles(ElementUnstyled.styles)(ElementUnstyled)

class PeriodicTable extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired
  }

  static styles = theme => ({
    root: {
      padding: theme.spacing.unit * 2,
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

  state = {
    selected: []
  }

  onElementClicked(element) {
    const index = this.state.selected.indexOf(element)
    const isClicked = index >= 0
    if (isClicked) {
      const selected = [...this.state.selected]
      selected.slice(index, 1)
      this.setState({selected: selected})
    } else {
      this.setState({selected: [element, ...this.state.selected]})
    }
  }

  render() {
    const {classes} = this.props
    return (
      <Paper className={classes.root}>
        <table className={classes.table}>
          <tbody>
            {elements.map((row, i) => (
              <tr key={i}>
                {row.map((element, j) => (
                  <td key={j}>
                    {element
                      ? <Element
                        element={element}
                        onClick={() => this.onElementClicked(element)}
                        selected={this.state.selected.indexOf(element) >= 0}
                      /> : ''}
                  </td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
      </Paper>
    )
  }
}

export default withStyles(PeriodicTable.styles)(PeriodicTable)
