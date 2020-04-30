import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import { schema } from '../MetaInfoRepository'
import SectionCard from './SectionCard'
import { Definition } from './Definition'

class SectionUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    section: PropTypes.object.isRequired,
    state: PropTypes.arrayOf(PropTypes.object).isRequired,
    toggleDefinition: PropTypes.func,
    onClose: PropTypes.func
  }

  static styles = (theme) => ({
    container: {
      display: 'flex',
      flexDirection: 'row',
      flexWrap: 'nowrap',
      justifyContent: 'flex-start',
      alignItems: 'center',
      padding: theme.spacing(1)
    },
    children: {
      paddingLeft: theme.spacing(5),
      display: 'flex',
      flexDirection: 'column',
      flexWrap: 'nowrap',
      justifyContent: 'center'
    }
  })

  constructor(props) {
    super(props)
    this.handleHideFeatures = this.handleHideFeatures.bind(this)
    this.handleFeatureClick = this.handleFeatureClick.bind(this)

    this.features = this.computeFeatures()
  }

  computeFeatures() {
    const mTypeOrder = [schema.section, schema.value, schema.reference]
    const mTypeCmp = (t1, t2) => {
      const i1 = mTypeOrder.indexOf(t1)
      const i2 = mTypeOrder.indexOf(t2)
      if (i1 === i2) {
        return 0
      } else if (i1 < i2) {
        return -1
      } else {
        return 1
      }
    }
    const sort = (list) => {
      if (!list) {
        return []
      } else {
        list.sort((e1, e2) => mTypeCmp(e1.mType, e2.mType) || e1.name.localeCompare(e2.name))
        return list
      }
    }
    return sort(this.props.section.features || [])
  }

  visible(feature) {
    return this.props.state.findIndex(({definition}) => definition.name === feature.name) !== -1
  }

  getState(feature) {
    return this.props.state.find(({definition}) => definition === feature).state
  }

  handleFeatureClick(feature) {
    this.props.toggleDefinition(feature, !this.visible(feature))
  }

  handleHideFeatures(featuresToHide) {
    featuresToHide.forEach(feature => {
      this.props.toggleDefinition(feature, false)
    })
  }

  render() {
    const features = this.features.map(feature => ({feature, visible: this.visible(feature)}))
    const {classes} = this.props
    const inherited = {...this.props, classes: undefined}
    return (
      <div className={classes.container}>
        <SectionCard {...inherited} features={features} onHideFeatures={this.handleHideFeatures} onFeatureClick={this.handleFeatureClick}/>
        <div className={classes.children}>
          {features.map(({feature, visible}) => (
            visible
              ? <Definition {...inherited}
                key={feature.name}
                definition={feature}
                state={this.getState(feature)}
                onClose={() => this.handleHideFeatures([feature])}/>
              : ''
          ))}
        </div>
      </div>
    )
  }
}

const Section = withStyles(SectionUnstyled.styles)(SectionUnstyled)

export default Section
