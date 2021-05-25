/*
 * Copyright The NOMAD Authors.
 *
 * This file is part of NOMAD. See https://nomad-lab.eu for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
import React, {useCallback} from 'react'
import PropTypes from 'prop-types'
import periodicTableData from './PeriodicTableData'
import { makeStyles, Typography, Button, Tooltip, FormControlLabel, Checkbox } from '@material-ui/core'
import chroma from 'chroma-js'
import { nomadSecondaryColor } from '../../config.js'

const elements = []
for (var i = 0; i < 10; i++) {
  elements[i] = Array.apply(null, Array(18))
}
periodicTableData.elements.forEach(element => {
  elements[element.ypos - 1][element.xpos - 1] = element
  element.category = element.category.replace(' ', '')
})

const useElementStyles = makeStyles(theme => ({
  root: {
    position: 'relative'
  },
  button: {
    border: '1px solid',
    paddingTop: theme.spacing(1),
    paddingBottom: theme.spacing(1),
    paddingLeft: 0,
    paddingRight: 0,
    width: '100%',
    textAlign: 'center',
    fontSize: '1rem',
    fontWeight: 700,
    textTransform: 'none',
    minWidth: 0,
    minHeight: 0,
    borderRadius: 0,
    boxShadow: 'none'
  },
  containedPrimary: {
    backgroundColor: theme.palette.primary.dark,
    color: 'white'
  },
  number: {
    position: 'absolute',
    top: 2,
    left: 2,
    margin: 0,
    padding: 0,
    fontSize: 8,
    pointerEvents: 'none'
  },
  count: {
    position: 'absolute',
    bottom: 2,
    right: 2,
    margin: 0,
    padding: 0,
    fontSize: 8,
    pointerEvents: 'none'
  }
}))

const Element = React.memo(({
  element,
  selected,
  count,
  heatmapScale,
  onClick
}) => {
  const styles = useElementStyles()
  const buttonClasses = {
    root: styles.button,
    containedPrimary: styles.containedPrimary
  }
  const disabled = count <= 0

  const style = (count > 0) ? {
    backgroundColor: !selected ? heatmapScale(count).hex() : undefined,
    borderColor: '#555'
  } : undefined

  return (
    <div className={styles.root}>
      <Tooltip title={element.name}>
        <div>
          <Button
            disabled={disabled}
            classes={buttonClasses}
            style={style}
            onClick={onClick} variant="contained"
            color={selected ? 'primary' : 'default'}
          >
            {element.symbol}
          </Button>
        </div>
      </Tooltip>
      <Typography
        classes={{root: styles.number}} variant="caption"
        style={selected ? {color: 'white'} : disabled ? {color: '#BDBDBD'} : {}}>
        {element.number}
      </Typography>
      {count >= 0
        ? <Typography
          classes={{root: styles.count}} variant="caption"
          style={selected ? {color: 'white'} : disabled ? {color: '#BDBDBD'} : {}}>
          {count.toLocaleString()}
        </Typography> : ''
      }
    </div>
  )
})

Element.propTypes = {
  element: PropTypes.object.isRequired,
  onClick: PropTypes.func,
  selected: PropTypes.bool,
  count: PropTypes.number.isRequired,
  heatmapScale: PropTypes.func.isRequired
}

const useTableStyles = makeStyles(theme => ({
  root: {
    position: 'relative'
  },
  table: {
    margin: 'auto',
    width: '100%',
    minWidth: 500,
    maxWidth: 900,
    tableLayout: 'fixed',
    borderSpacing: theme.spacing(0.5)
  },
  formContainer: {
    position: 'absolute',
    top: theme.spacing(0),
    left: '10%',
    textAlign: 'center'
  }
}))

const NewPeriodicTable = React.memo(({
  aggregations,
  metric,
  values,
  exclusive,
  onExclusiveChanged,
  onChanged
}) => {
  const styles = useTableStyles()

  const onElementClicked = useCallback((element) => {
    if (values) {
      const isSelected = values?.has(element)
      isSelected ? values.delete(element) : values.add(element)
    } else {
      values = new Set()
      values.add(element)
    }
    onChanged(values)
  }, [values, onChanged])

  const unSelectedAggregations = useCallback(() => {
    return Object.keys(aggregations)
      .filter(key => !values?.has(key))
      .map(key => aggregations[key][metric])
  })

  const max = aggregations ? Math.max(...unSelectedAggregations()) || 1 : 1
  const heatmapScale = chroma.scale([nomadSecondaryColor.veryLight, nomadSecondaryColor.main]).domain([1, max], 10, 'log')

  return (
    <div className={styles.root}>
      <table className={styles.table}>
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
                      onClick={() => onElementClicked(element.symbol)}
                      selected={values?.has(element.symbol)}
                    /> : ''}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
      <div className={styles.formContainer}>
        <Tooltip title={
          'Search for entries with compositions that only (exclusively) contain the ' +
          'selected atoms. The default is to return all entries that have at least ' +
          '(inclusively) the selected atoms.'}>
          <FormControlLabel
            control={<Checkbox checked={exclusive} onChange={onExclusiveChanged} />}
            label={'only composition that exclusively contain these atoms'}
          />
        </Tooltip>
      </div>
    </div>
  )
})

NewPeriodicTable.propTypes = {
  aggregations: PropTypes.object,
  metric: PropTypes.string.isRequired,
  values: PropTypes.object,
  onChanged: PropTypes.func.isRequired,
  exclusive: PropTypes.bool,
  onExclusiveChanged: PropTypes.func.isRequired
}

export default NewPeriodicTable
