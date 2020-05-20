import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import { List } from '@material-ui/core'
import SectionFeature from './SectionFeature'
import { AfterRenderMeasure } from './util/after-render'
import { schema } from '../MetaInfoRepository'
import { CheckChip } from './ui'
import { DefinitionCard } from './DefinitionCard'
import { updateList } from './util/data'
import { CardCompartment } from './util/cards'

class SectionCardUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    section: PropTypes.object.isRequired,
    /* called if user closes the whole section */
    onClose: PropTypes.func,
    /* called if user clicks feature */
    onFeatureClick: PropTypes.func,
    /* called if features must become invisible due to hidden categories */
    onHideFeatures: PropTypes.func,
    /* all features as {feature, visibilty} of the section in the order they should appear */
    features: PropTypes.arrayOf(PropTypes.object).isRequired,

    toggleDefinition: PropTypes.func.isRequired
  }

  static styles = (theme) => ({
    categories: {
      paddingTop: theme.spacing(1),
      paddingLeft: theme.spacing(1),
      paddingRight: theme.spacing(1),
      lineHeight: 1.5
    }
  })

  constructor(props) {
    super(props)
    const {features} = this.props

    const categories = {
      'OTHERS': []
    }

    features
      .map(({feature}) => feature)
      .filter(schema.isProperty)
      .forEach(property => {
        if (property.categories && property.categories.length > 0) {
          const visit = (category) => {
            if (category.super && category.super.length > 0) {
              category.super.forEach(visit)
            } else {
              categories[category.name] = categories[category.name] || []
              categories[category.name].push(property)
            }
          }
          property.categories.forEach(visit)
        } else {
          categories['OTHERS'].push(property)
        }
      })

    this.state = {
      categories: Object.keys(categories).map(name => ({
        name: name,
        features: categories[name],
        visible: false
      }))
    }
  }

  onToggleCategory(name, visible) {
    const index = this.state.categories.findIndex(state => state.name === name)
    this.setState({
      categories: updateList(this.state.categories, index, category => ({
        ...category,
        visible: visible
      }))
    })

    if (!visible) {
      // fire toggle visible=false events for features of the new hidden category
      this.props.onHideFeatures(this.state.categories[index].features)
    }
  }

  allCategoriesSelected() {
    return this.state.categories.find(({name, visible}) => !visible) === false
  }

  onToggleAllCategories(visible) {
    this.setState({categories: this.state.categories.map(category => ({...category, visible: visible}))})

    if (!visible) {
      // fire toggle visible=false events for features of the new hidden category
      const allFeaturesToHide = this.props.features.filter(({feature}) => schema.isProperty(feature)).map(({feature}) => feature)
      this.props.onHideFeatures(allFeaturesToHide)
    }
  }

  sections() {
    return this.props.features.filter(({feature}) => schema.isSection(feature))
  }

  properties() {
    return this.props.features
      .filter(({feature}) => this.state.categories.find(({features, visible}) => visible && features.indexOf(feature) !== -1))
  }

  renderFeatures(features) {
    const renderedFeatures = features.map(({feature, visible}, index) => (
      <AfterRenderMeasure measureId={`feature.${feature.name}`} measureData={feature} key={index}>
        <SectionFeature feature={feature}
          visible={visible}
          onClick={() => this.props.onFeatureClick(feature)}/>
      </AfterRenderMeasure>
    ))

    return (
      <List dense={true}>{renderedFeatures}</List>
    )
  }

  render() {
    const { classes, section } = this.props
    const inherited = {...this.props, classes: undefined}
    return (
      <AfterRenderMeasure measureId={`section.${section.name}`} measureData={section}>
        <DefinitionCard {...inherited} definition={section}>
          {section.features.find(schema.isSection)
            ? <AfterRenderMeasure measureId={`sections.${section.name}`}>
              <CardCompartment label="sections" foldable compId="sections">
                {this.renderFeatures(this.sections())}
              </CardCompartment>
            </AfterRenderMeasure> : ''
          }
          {section.features.find(schema.isProperty)
            ? <AfterRenderMeasure measureId={`properties.${section.name}`}>
              <CardCompartment label={'quantities'} foldable compId="features">
                <div className={classes.categories}>
                  <CheckChip label={`ALL (${section.features.filter(schema.isProperty).length})`} checked={this.allCategoriesSelected()}
                    onChange={(value) => this.onToggleAllCategories(value)}/>
                  {this.state.categories.length > 1
                    ? this.state.categories.map(({name, features, visible}, index) => (
                      <CheckChip key={index} label={`${name} (${features.length})`} checked={visible}
                        onChange={(value) => this.onToggleCategory(name, value)}/>
                    )) : ''}
                </div>
                <div>
                  {this.renderFeatures(this.properties())}
                </div>
              </CardCompartment>
            </AfterRenderMeasure> : ''
          }
        </DefinitionCard>
      </AfterRenderMeasure>
    )
  }
}

const SectionCard = withStyles(SectionCardUnstyled.styles)(SectionCardUnstyled)

export default SectionCard
