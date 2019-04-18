import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import { Typography, IconButton, Icon, Popover } from '@material-ui/core'
import grey from '@material-ui/core/colors/grey'
import { updateListState } from './data'

const ToggleCardCompartmentContext = React.createContext()

export class CardButton extends React.Component {
  static propTypes = {
    position: PropTypes.string
  }
  render() {
    return (
      <MIButton size="tiny" {...{...this.props, position: undefined}}/>
    )
  }
}

class CardUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    children: PropTypes.any
  }
  static styles = (theme) => ({
    root: {
      width: 300,
      position: 'relative'
    },
    type: {

    },
    header: {
      position: 'relative'
    },
    left: {
      position: 'absolute',
      top: 0,
      left: 0
    },
    center: {
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'center'
    },
    right: {
      position: 'absolute',
      top: 0,
      right: 0
    }
  })

  constructor(props) {
    super(props)
    this.handleMaximize = this.handleMaximize.bind(this)
    this.handleMinimize = this.handleMinimize.bind(this)
    this.handleToggleCompartment = this.handleToggleCompartment.bind(this)
    this.compartmentIsVisible = this.compartmentIsVisible.bind(this)

    this.state = {
      compartments: []
    }
  }

  handleMinimize() {
    this.setState({
      compartments: this.state.compartments.map(state => ({
        ...state,
        visible: false
      }))
    })
  }

  handleMaximize() {
    this.setState({
      compartments: this.state.compartments.map(state => ({
        ...state,
        visible: true
      }))
    })
  }

  handleToggleCompartment(compartment) {
    const state = updateListState(this.state.compartments,
      state => state.compartment.props.compId === compartment.props.compId,
      state => ({...state, visible: !state.visible}),
      ({compartment: compartment, visible: !compartment.props.folded}))
    this.setState({compartments: state})
  }

  compartmentIsVisible(compartment, visibleDefault) {
    const state = this.state.compartments.find(state => {
      return state.compartment.props.compId === compartment.props.compId
    })
    if (state) {
      return state.visible
    } else {
      const state = this.state // this is a deliberate hack
      state.compartments = [...this.state.compartments, {compartment: compartment, visible: visibleDefault}]
      return visibleDefault
    }
  }

  getButtons(position, isDefault) {
    const components = [CardButton, PopoverCardButton]
    const children = React.Children.map(this.props.children, child => child)
      .filter(child => {
        return components.find(component => child.type === component)
      })
    if (position) {
      return children.filter(child => {
        if (child.props.position) {
          return child.props.position === position
        } else {
          return isDefault
        }
      })
    } else {
      return children
    }
  }

  getChildren() {
    const components = [CardButton, PopoverCardButton]
    return React.Children.map(this.props.children, child => child)
      .filter(child => {
        return !components.find(component => child.type === component)
      })
  }

  render() {
    const { classes } = this.props
    return (
      <div className={classes.card}>
        <div className={classes.header}>
          <div className={classes.left}>
            {this.getButtons('left')}
          </div>
          <div className={classes.center}>
            {this.getButtons('center', true)}
          </div>
          <div className={classes.right}>
            <MIButton icon="minimize" size="tiny" onClick={this.handleMinimize}/>
            <MIButton icon="add" size="tiny" onClick={this.handleMaximize}/>
            {this.getButtons('right')}
          </div>
        </div>
        <ToggleCardCompartmentContext.Provider value={{
          onToggleCompartment: this.handleToggleCompartment,
          isVisible: this.compartmentIsVisible}}>
          {this.getChildren()}
        </ToggleCardCompartmentContext.Provider>
      </div>
    )
  }
}

export const Card = withStyles(CardUnstyled.styles)(CardUnstyled)

class CardCompartmentUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    compId: PropTypes.string.isRequired,
    label: PropTypes.string,
    padded: PropTypes.bool,
    borderless: PropTypes.bool,
    foldable: PropTypes.bool,
    folded: PropTypes.bool,
    children: PropTypes.any
  }
  static styles = (theme) => ({
    root: {
      borderTop: `1px solid ${grey[400]}`,
      overflowX: 'hidden'
    },
    borderless: {
      overflowX: 'hidden'
    },
    header: {
      position: 'relative'
    },
    minHeader: {
      position: 'relative',
      minHeight: 16 + theme.spacing.unit
    },
    foldButton: {
      position: 'absolute',
      top: 0,
      right: 0,
      margin: theme.spacing.unit * 0.5
    },
    content: {

    },
    padded: {
      padding: theme.spacing.unit
    },
    label: {
      fontSize: 10,
      textTransform: 'uppercase',
      paddingTop: theme.spacing.unit * 0.5,
      paddingLeft: theme.spacing.unit * 0.5,
      color: grey[600]
    }
  })

  render() {
    return (
      <ToggleCardCompartmentContext.Consumer>
        {context => this.renderCompartment(
          () => context.onToggleCompartment(this),
          context.isVisible(this, !this.props.folded)
        )}
      </ToggleCardCompartmentContext.Consumer>
    )
  }

  renderCompartment(onToggleCompartment, visible) {
    const {classes, foldable, label, borderless, children, padded} = this.props
    return (
      <div className={borderless ? classes.borderless : classes.root}>
        {label || foldable
          ? <div className={visible ? classes.header : classes.minHeader}>
            {this.props.label
              ? <Typography classes={{root: classes.label}}>{label}</Typography>
              : ''}
            {foldable
              ? <MIButton size="tiniest" classes={{tiniest: classes.foldButton}}
                onClick={onToggleCompartment} icon={visible ? 'minimize' : 'add'}/>
              : ''}
          </div>
          : ''}
        {visible || !foldable
          ? <div className={padded ? classes.padded : classes.content}>
            {children}
          </div> : ''
        }
      </div>
    )
  }
}

export const CardCompartment = withStyles(CardCompartmentUnstyled.styles)(CardCompartmentUnstyled)

class PopoverCardButtonUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    children: PropTypes.object,
    label: PropTypes.string,
    icon: PropTypes.string,
    size: PropTypes.string
  }

  static styles = theme => ({
    root: {
      display: 'inline'
    },
    content: {
      padding: theme.spacing.unit
    },
    iconButton: {

    }
  })

  state = {
    anchorEl: null
  };

  handleClick = event => {
    this.setState({
      anchorEl: event.currentTarget
    })
  };

  handleClose = () => {
    this.setState({
      anchorEl: null
    })
  };

  render() {
    const { classes, size } = this.props
    const { anchorEl } = this.state

    return (
      <div className={classes.root}>
        <CardButton {...this.props} variant={'contained'} onClick={this.handleClick} classes={{[size || 'small']: classes.iconButton}}/>
        <Popover
          classes={{paper: classes.content}}
          open={Boolean(anchorEl)}
          anchorEl={anchorEl}
          onClose={this.handleClose}
          anchorOrigin={{
            vertical: 'bottom',
            horizontal: 'center'
          }}
          transformOrigin={{
            vertical: 'top',
            horizontal: 'center'
          }}
        >
          {this.props.children}
        </Popover>
      </div>
    )
  }
}

export const PopoverCardButton = withStyles(PopoverCardButtonUnstyled.styles)(PopoverCardButtonUnstyled)

class MIButtonUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    size: PropTypes.string,
    icon: PropTypes.string.isRequired,
    onClick: PropTypes.func
  }
  static styles = theme => ({
    small: {
      padding: 3
    },
    tiny: {
      padding: 3,
      '& .material-icons': {
        fontSize: 18
      }
    },
    tiniest: {
      padding: 3,
      '& .material-icons': {
        fontSize: 12
      }
    }
  })
  render() {
    const {classes, size, onClick, icon} = this.props
    let buttonClass = classes.small
    if (!size) {
      buttonClass = classes.small
    } else if (size === 'tiny') {
      buttonClass = classes.tiny
    } else if (size === 'tiniest') {
      buttonClass = classes.tiniest
    }
    return (
      <IconButton size="small" {...this.props} classes={{root: buttonClass}}
        onClick={onClick}>
        <Icon>{icon}</Icon>
      </IconButton>)
  }
}

const MIButton = withStyles(MIButtonUnstyled.styles)(MIButtonUnstyled)
