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
import React, { useState, useCallback, useMemo } from 'react'
import PropTypes from 'prop-types'
import { string, bool } from 'yup'
import clsx from 'clsx'
import { isNil } from 'lodash'
import {
  TextField,
  Typography,
  MenuItem,
  FormControlLabel,
  Checkbox,
  makeStyles
} from '@material-ui/core'
import { useResizeDetector } from 'react-resize-detector'
import { useSearchContext } from '../SearchContext'
import { InputMetainfo } from '../input/InputMetainfo'
import { InputTextQuantity, InputTextField } from '../input/InputText'
import InputItem, { inputItemHeight } from '../input/InputItem'
import InputUnavailable from '../input/InputUnavailable'
import InputTooltip from '../input/InputTooltip'
import Placeholder from '../../visualization/Placeholder'
import { Widget, schemaWidget } from './Widget'
import { ActionSelect } from '../../Actions'
import { WidgetEditDialog, WidgetEditGroup, WidgetEditOption } from './WidgetEdit'
import { DType, pluralize } from '../../../utils'
import { scales } from '../../plotting/common'

// Predefined in order to not break memoization
const dtypes = new Set([DType.String, DType.Enum])

/**
 * Displays a terms widget.
 */
const useStyles = makeStyles(theme => ({
  root: {
  },
  outerContainer: {
    position: 'relative',
    width: '100%',
    height: '100%'
  },
  innerContainer: {
    position: 'absolute',
    top: 0,
    bottom: 0,
    right: 0,
    left: 0,
    display: 'flex',
    flexDirection: 'column'
  },
  menuItem: {
    height: inputItemHeight
  },
  menuItemSelected: {
    '&.Mui-selected': {
      backgroundColor: 'transparent'
    }
  },
  chips: {
    display: 'flex',
    flexWrap: 'wrap'
  },
  icon: {
    right: theme.spacing(1)
  },
  spacer: {
    flex: 1,
    minHeight: 0
  },
  textField: {
    marginBottom: theme.spacing(1)
  },
  count: {
    height: '1.2rem',
    marginTop: theme.spacing(0.25),
    marginRight: theme.spacing(1),
    display: 'flex',
    justifyContent: 'flex-end',
    alignItems: 'center',
    color: theme.palette.text.disabled
  }
}))
export const WidgetTerms = React.memo((
{
  id,
  title,
  description,
  quantity,
  scale,
  showinput,
  className,
  'data-testid': testID
}) => {
  const {useAgg, useFilterState, filterData} = useSearchContext()
  const styles = useStyles()
  const [filter, setFilter] = useFilterState(quantity)
  const { height, ref } = useResizeDetector()
  const { useSetWidget } = useSearchContext()
  const setWidget = useSetWidget(id)

  // The terms aggregations need to request at least 1 item or an API error is thrown
  const aggSize = useMemo(() => Math.max(Math.floor(height / inputItemHeight), 1), [height])
  const aggConfig = useMemo(() => {
    const config = {type: 'terms', size: aggSize}
    // If a fixed list of options is used, we must restrict the aggregation
    // return values with 'include'. Otherwise the returned results may not
    // contain the correct values.
    const options = filterData[quantity]?.options
    if (options) config.include = Object.keys(options)
    return config
  }, [aggSize, filterData, quantity])
  const agg = useAgg(quantity, !isNil(height), id, aggConfig)
  const max = agg ? Math.max(...agg.data.map(option => option.count)) : 0

  const handleChange = useCallback((event, key, selected) => {
    setFilter(old => {
      if (!old) return new Set([key])
      const newValue = new Set(old)
      selected ? newValue.add(key) : newValue.delete(key)
      return newValue
    })
  }, [setFilter])

  const handleEdit = useCallback(() => {
    setWidget(old => { return {...old, editing: true } })
  }, [setWidget])

  const handleChangeScale = useCallback((value) => {
    setWidget(old => { return {...old, scale: value} })
  }, [setWidget])

  const [aggComp, nShown] = useMemo(() => {
    let aggComp
    let nShown = 0
    if (!agg) {
      aggComp = <Placeholder
        variant="rect"
        data-testid={`${testID}-placeholder`}
        margin={0}
      />
      nShown = 0
    } else {
      if (agg?.data && agg.data.length > 0) {
        aggComp = []
        const maxSize = Math.min(aggConfig.size, agg.data.length)
        for (let i = 0; i < maxSize; ++i) {
          const option = agg.data[i]
          if (option.count > 0 && nShown < maxSize) {
            aggComp.push(<InputItem
              key={option.value}
              value={option.value}
              selected={filter ? filter.has(option.value) : false}
              max={max}
              onChange={handleChange}
              variant="checkbox"
              count={option.count}
              scale={scale}
            />)
            ++nShown
          }
        }
      } else {
        aggComp = <InputUnavailable/>
        nShown = 0
      }
    }
    return [aggComp, nShown]
  }, [agg, aggConfig, filter, handleChange, max, scale, testID])

  const count = pluralize('item', nShown, true)

  return <Widget
    id={id}
    quantity={quantity}
    title={title}
    description={description}
    onEdit={handleEdit}
    className={clsx(className)}
    actions={
      <ActionSelect
        value={scale}
        options={Object.keys(scales)}
        tooltip="Statistics scaling"
        onChange={handleChangeScale}
      />
    }
  >
    <InputTooltip>
      <div className={clsx(styles.outerContainer)}>
        <div className={clsx(styles.innerContainer)}>
          {showinput
            ? <InputTextQuantity
                className={styles.textField}
                quantity={quantity}
                disabled={false}
                disableSuggestions={false}
                fullWidth
              />
            : null
          }
          <div ref={ref} className={styles.spacer}>
            {aggComp}
          </div>
          {nShown !== 0 && <div className={styles.count}>
            <Typography variant="overline">
              {nShown < aggConfig.size
                ? nShown === 1
                  ? 'Showing the only item'
                  : `Showing all ${count}`
                : `Showing top ${count}`
              }
            </Typography>
          </div>}
        </div>
      </div>
    </InputTooltip>
  </Widget>
})

WidgetTerms.propTypes = {
  id: PropTypes.string.isRequired,
  title: PropTypes.string,
  description: PropTypes.string,
  quantity: PropTypes.string,
  nbins: PropTypes.number,
  scale: PropTypes.string,
  autorange: PropTypes.bool,
  showinput: PropTypes.bool,
  className: PropTypes.string,
  'data-testid': PropTypes.string
}

WidgetTerms.defaultProps = {
  'data-testid': 'widgetterms'
}

/**
 * A dialog that is used to configure a scatter plot widget.
 */
export const WidgetTermsEdit = React.memo((props) => {
    const {id, editing, visible} = props
    const { useSetWidget } = useSearchContext()
    const [settings, setSettings] = useState(props)
    const [errors, setErrors] = useState({})
    const setWidget = useSetWidget(id)
    const hasError = useMemo(() => {
      return Object.values(errors).some((d) => !!d) || !schemaWidgetTerms.isValidSync(settings)
    }, [errors, settings])

    const handleSubmit = useCallback((settings) => {
      setWidget(old => ({...old, ...settings}))
    }, [setWidget])

    const handleChange = useCallback((key, value) => {
      setSettings(old => ({...old, [key]: value}))
    }, [setSettings])

    const handleError = useCallback((key, value) => {
      setErrors(old => ({...old, [key]: value}))
    }, [setErrors])

    const handleClose = useCallback(() => {
      setWidget(old => ({...old, editing: false}))
    }, [setWidget])

    const handleAccept = useCallback((key, value) => {
      try {
        schemaWidgetTerms.validateSyncAt(key, {[key]: value})
      } catch (e) {
        handleError(key, e.message)
        return
      }
      setErrors(old => ({...old, [key]: undefined}))
      setSettings(old => ({...old, [key]: value}))
    }, [handleError, setSettings])

    const handleEditAccept = useCallback(() => {
      handleSubmit({...settings, editing: false, visible: true})
    }, [handleSubmit, settings])

    return <WidgetEditDialog
        id={id}
        open={editing}
        visible={visible}
        title="Edit terms widget"
        onClose={handleClose}
        onAccept={handleEditAccept}
        error={hasError}
      >
      <WidgetEditGroup title="x axis">
        <WidgetEditOption>
          <InputMetainfo
            label="quantity"
            value={settings.quantity}
            error={errors.quantity}
            onChange={(value) => handleChange('quantity', value)}
            onSelect={(value) => handleAccept('quantity', value)}
            onError={(value) => handleError('quantity', value)}
            dtypes={dtypes}
            dtypesRepeatable={dtypes}
            disableNonAggregatable
          />
        </WidgetEditOption>
        <WidgetEditOption>
          <TextField
            select
            fullWidth
            label="Statistics scaling"
            variant="filled"
            value={settings.scale}
            onChange={(event) => { handleChange('scale', event.target.value) }}
          >
            {Object.keys(scales).map((key) =>
              <MenuItem value={key} key={key}>{key}</MenuItem>
            )}
          </TextField>
        </WidgetEditOption>
      </WidgetEditGroup>
      <WidgetEditGroup title="general">
        <WidgetEditOption>
          <InputTextField
            label="title"
            fullWidth
            value={settings?.title}
            onChange={(event) => handleChange('title', event.target.value)}
          />
        </WidgetEditOption>
        <WidgetEditOption>
          <FormControlLabel
            control={<Checkbox checked={settings.showinput} onChange={(event, value) => handleChange('showinput', value)}/>}
            label='Show input field'
          />
        </WidgetEditOption>
      </WidgetEditGroup>
    </WidgetEditDialog>
})

WidgetTermsEdit.propTypes = {
  id: PropTypes.string.isRequired,
  editing: PropTypes.bool,
  visible: PropTypes.bool,
  quantity: PropTypes.string,
  scale: PropTypes.string,
  nbins: PropTypes.number,
  autorange: PropTypes.bool,
  showinput: PropTypes.bool,
  onClose: PropTypes.func
}

export const schemaWidgetTerms = schemaWidget.shape({
  quantity: string().required('Quantity is required.'),
  scale: string().required('Scale is required.'),
  showinput: bool()
})
