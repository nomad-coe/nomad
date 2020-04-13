import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import grey from '@material-ui/core/colors/grey'
import { Chip } from '@material-ui/core'

class LimitedTextUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    text: PropTypes.string.isRequired,
    limit: PropTypes.number.isRequired,
    render: PropTypes.func
  }

  static style = (theme) => ({
    root: {
    },
    withLimit: {
      cursor: 'pointer'
    }
  })

  constructor(props) {
    super(props)
    this.handleMoreClick = this.handleMoreClick.bind(this)
    this.limittedText = this.limittedText()
    this.state = {
      shownText: this.limittedText
    }
  }

  isLimited() {
    return this.state.shownText !== this.props.text
  }

  limittedText() {
    if (this.props.text.length > this.props.limit) {
      this.withLimit = true
      const words = this.props.text.split(' ')
      let result = ''
      words.forEach(word => {
        if (result.length + word.length <= this.props.limit) {
          result += ' ' + word
        }
      })
      return result
    } else {
      this.withLimit = false
      return this.props.text
    }
  }

  handleMoreClick() {
    if (this.isLimited()) {
      this.setState({shownText: this.props.text})
    } else {
      this.setState({shownText: this.limittedText})
    }
  }

  renderText() {
    if (this.props.render) {
      return this.props.render(this.state.shownText)
    } else {
      return this.state.shownText
    }
  }

  render() {
    const {classes} = this.props
    return (
      <div className={this.withLimit ? classes.withLimit : classes.root} onClick={this.handleMoreClick}>
        {this.renderText()}
        {this.isLimited()
          ? <span>...</span> : ''
        }
      </div>
    )
  }
}

export const LimitedText = withStyles(LimitedTextUnstyled.style)(LimitedTextUnstyled)

class CheckChipUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    checked: PropTypes.bool,
    onChange: PropTypes.func,
    label: PropTypes.string
  }
  static styles = (theme) => ({
    chip: {
      height: 20,
      borderRadius: 10,
      marginRight: theme.spacing(1),
      '&:focus': {
        background: grey[300]
      }
    },
    checkedChip: {
      height: 20,
      borderRadius: 10,
      marginRight: theme.spacing(1),
      color: 'white',
      backgroundColor: theme.palette.primary.main,
      '&:focus': {
        background: theme.palette.primary.main
      },
      '&:hover': {
        background: theme.palette.primary.light
      }
    },
    label: {
      fontSize: 12,
      paddingLeft: 8,
      paddingRight: 8
    }
  })

  constructor(props) {
    super(props)
    this.state = {
      checked: this.props.checked
    }
    this.onClick = this.onClick.bind(this)
  }

  componentDidUpdate(prevProps) {
    if (this.props.checked !== prevProps.checked && this.props.checked !== this.state.checked) {
      this.setState({checked: this.props.checked})
    }
  }

  onClick() {
    const newState = !this.state.checked
    this.setState({checked: newState})
    if (this.props.onChange) {
      this.props.onChange(newState)
    }
  }

  render() {
    const {classes} = this.props
    const {checked} = this.state
    return (
      <Chip classes={{root: checked ? classes.checkedChip : classes.chip, label: classes.label}}
        onClick={this.onClick} label={this.props.label}/>
    )
  }
}

export const CheckChip = withStyles(CheckChipUnstyled.styles)(CheckChipUnstyled)
