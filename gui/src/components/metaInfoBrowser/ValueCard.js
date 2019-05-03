import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import { AfterRenderMeasure } from './util/after-render'
import { Typography } from '@material-ui/core'
import { DefinitionCard, renderName } from './DefinitionCard'
import { schema } from '../MetaInfoRepository'
import { CardCompartment, CardButton } from './util/cards'

export const MetaAttribute = (props) => {
  return (
    <Typography component={'p'} variant={'body1'}>
      <span>{props.label}</span>: <span style={{fontWeight: 'bold'}}>{props.value}</span>
    </Typography>
  )
}

MetaAttribute.propTypes = {
  label: PropTypes.string.isRequired,
  value: PropTypes.any
}

export const ValueAttributes = props => {
  const { definition } = props
  return <div>
    <MetaAttribute label={'repeats'} value={'' + (!!definition.miJson.repeats)}/>
    <MetaAttribute label={'shape'} value={`[${definition.miJson.shape.join(', ')}]`}/>
    <MetaAttribute label={'type'} value={definition.type}/>
    {definition.miJson.units ? <MetaAttribute label={'units'} value={definition.miJson.units}/> : ''}
  </div>
}

ValueAttributes.propTypes = {
  definition: PropTypes.object.isRequired
}

class ValueUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    definition: PropTypes.object.isRequired,
    onClose: PropTypes.func,
    toggleDefinition: PropTypes.func
  }
  static styles = (theme) => ({
    root: {
      padding: theme.spacing.unit
    }
  })
  render() {
    const {classes, definition, toggleDefinition} = this.props
    const inherited = {...this.props, classes: undefined}
    return (
      <div className={classes.root}>
        <AfterRenderMeasure measureId={`property.${definition.name}`} measureData={definition}>
          <DefinitionCard {...inherited}>
            <CardCompartment compId="type_info" padded label={'type info'} foldable>
              <ValueAttributes definition={definition} />
              {schema.isReference(definition) && definition.referencedSection
                ? <MetaAttribute label={'referenced'} value={<span>
                  {renderName(definition.referencedSection)} <CardButton size="tiny" icon="arrow_right_alt"
                    onClick={() => toggleDefinition(definition.referencedSection)}/>
                </span>}/> : ''
              }
            </CardCompartment>
          </DefinitionCard>
        </AfterRenderMeasure>
      </div>
    )
  }
}

const Value = withStyles(ValueUnstyled.styles)(ValueUnstyled)

export default Value
