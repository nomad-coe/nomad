import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import Typography from '@material-ui/core/Typography'
import { ListItem, ListItemIcon, Icon } from '@material-ui/core'
import grey from '@material-ui/core/colors/grey'
import { schema } from '../MetaInfoRepository'
import { renderName } from './DefinitionCard'

class SectionFeature extends React.Component {
  static styles = (theme) => ({
    root: {

    },
    visible: {
      extends: 'root',
      backgroundColor: grey[300]
    },
    gutters: {
      paddingLeft: theme.spacing(1),
      paddingRight: theme.spacing(1)
    },
    stateIcon: {
      paddingLeft: theme.spacing(1)
    }
  })

  static propTypes = {
    classes: PropTypes.object.isRequired,
    feature: PropTypes.object.isRequired,
    visible: PropTypes.bool,
    onClick: PropTypes.func
  }

  renderIcon() {
    const {feature} = this.props
    if (schema.isSection(feature)) {
      return <Icon>code</Icon>
    } else if (schema.isValue(feature)) {
      return <Icon>fiber_manual_record</Icon>
    } else if (schema.isReference(feature)) {
      return <Icon>arrow_right_alt</Icon>
    } else {
      return <Icon>help</Icon>
    }
  }

  render() {
    const {feature, classes} = this.props
    return (
      <ListItem
        classes={{gutters: classes.gutters}}
        className={this.props.visible ? classes.visible : classes.root}
        dense button onClick={this.props.onClick}>
        <ListItemIcon>
          {this.renderIcon()}
        </ListItemIcon>
        <Typography>{renderName(feature)}</Typography>
      </ListItem>
    )
  }
}

export default withStyles(SectionFeature.styles)(SectionFeature)
