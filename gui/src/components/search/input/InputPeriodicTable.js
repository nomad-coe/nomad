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
import React, { useCallback } from 'react'
import PropTypes from 'prop-types'
import periodicTableData from './PeriodicTableData'
import {
  Typography,
  Button,
  Tooltip
} from '@material-ui/core'
import InputCheckbox from './InputCheckbox'
import { useTheme, makeStyles } from '@material-ui/core/styles'

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
    width: '100%',
    height: '100%',
    border: '1px solid',
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
  containerOuter: {
    width: '100%',
    paddingBottom: '100%',
    position: 'relative'
  },
  containerInner: {
    position: 'absolute',
    width: '100%',
    height: '100%'
  },
  buttonRoot: {
    padding: 16,
    paddingTop: 17,
    paddingBottom: 15,
    width: '100%',
    height: '100%',
    boxSizing: 'border-box',
    minWidth: 0,
    border: '1px solid',
    textAlign: 'center',
    fontSize: '1rem',
    fontWeight: 600,
    textTransform: 'none',
    minHeight: 0,
    borderRadius: 0,
    boxShadow: 'none'
  },
  number: {
    position: 'absolute',
    top: 0,
    left: 2,
    margin: 0,
    padding: 0,
    fontSize: 8,
    pointerEvents: 'none'
  }
}))

const Element = React.memo(({
  element,
  selected,
  disabled,
  onClick
}) => {
  const styles = useElementStyles()
  const buttonClasses = {
    root: styles.buttonRoot,
    containedPrimary: styles.containedPrimary
  }
  const theme = useTheme()

  const style = (!disabled) ? {
    backgroundColor: !selected ? theme.palette.secondary.light : undefined,
    borderColor: '#555'
  } : undefined

  return (
    <div className={styles.root}>
      <Tooltip title={element.name}>
        <div className={styles.containerOuter}>
          <div className={styles.containerInner}>
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
        </div>
      </Tooltip>
      <Typography
        classes={{root: styles.number}} variant="caption"
        style={selected ? {color: 'white'} : disabled ? {color: '#BDBDBD'} : {}}>
        {element.number}
      </Typography>
      {!disabled >= 0
        ? <Typography
          classes={{root: styles.count}} variant="caption"
          style={selected ? {color: 'white'} : disabled ? {color: '#BDBDBD'} : {}}>
        </Typography> : ''
      }
    </div>
  )
})

Element.propTypes = {
  element: PropTypes.object.isRequired,
  onClick: PropTypes.func,
  selected: PropTypes.bool,
  disabled: PropTypes.bool
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

const InputPeriodicTable = React.memo(({
  availableValues,
  values,
  onChanged
}) => {
  const styles = useTableStyles()

  const onElementClicked = useCallback((element) => {
    let newValues
    if (values) {
      const isSelected = values?.has(element)
      isSelected ? values.delete(element) : values.add(element)
      newValues = new Set(values)
    } else {
      newValues = new Set()
      newValues.add(element)
    }
    onChanged(newValues)
  }, [values, onChanged])

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
                      disabled={!availableValues[element.symbol] && !values?.has(element.symbol)}
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
        <InputCheckbox
          quantity="exclusive"
          label="only composition that exclusively contain these atoms"
          description="Search for entries with compositions that only (exclusively) contain the selected atoms. The default is to return all entries that have at least (inclusively) the selected atoms."
          initialValue={false}
        ></InputCheckbox>
        {/* <Tooltip title={
          'Search for entries with compositions that only (exclusively) contain the ' +
          'selected atoms. The default is to return all entries that have at least ' +
          '(inclusively) the selected atoms.'}>
          <FormControlLabel
            control={<Checkbox
              checked={exclusive}
              onChange={(event) => { onExclusiveChanged(event.target.checked) }}
            />}
            label={'only composition that exclusively contain these atoms'}
          />
        </Tooltip> */}
      </div>
    </div>
  )
})

InputPeriodicTable.propTypes = {
  availableValues: PropTypes.object,
  values: PropTypes.object,
  onChanged: PropTypes.func.isRequired
}

export default InputPeriodicTable
