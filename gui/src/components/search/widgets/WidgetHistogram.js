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
import { string, number, bool } from 'yup'
import {
  TextField,
  MenuItem,
  Checkbox,
  FormControlLabel
} from '@material-ui/core'
import { useSearchContext } from '../SearchContext'
import { InputMetainfo } from '../input/InputMetainfo'
import { InputTextField } from '../input/InputText'
import { Widget, schemaWidget } from './Widget'
import { ActionCheckbox, ActionSelect } from '../../Actions'
import { WidgetEditDialog, WidgetEditGroup, WidgetEditOption } from './WidgetEdit'
import { Range } from '../input/InputRange'
import { DType } from '../../../utils'
import { scales } from '../../plotting/common'

// Predefined in order to not break memoization
const dtypes = new Set([DType.Float, DType.Int, DType.Timestamp])

/**
 * Displays a histogram widget.
 */
export const autorangeDescription = 'Automatically center the view on the data'
export const WidgetHistogram = React.memo((
{
  id,
  title,
  description,
  quantity,
  nbins,
  scale,
  autorange,
  showinput,
  className
}) => {
  const { useSetWidget } = useSearchContext()
  const setWidget = useSetWidget(id)

  const handleEdit = useCallback(() => {
    setWidget(old => { return {...old, editing: true } })
  }, [setWidget])

  const handleChangeScale = useCallback((value) => {
    setWidget(old => { return {...old, scale: value} })
  }, [setWidget])

  return <Widget
    id={id}
    quantity={quantity}
    title={title}
    description={description}
    onEdit={handleEdit}
    className={className}
    actions={<>
      <ActionCheckbox
        tooltip={autorangeDescription}
        label="autorange"
        value={autorange}
        onChange={(value) => setWidget(old => ({...old, autorange: value}))}
      />
      <ActionSelect
        value={scale}
        options={Object.keys(scales)}
        tooltip="Statistics scaling"
        onChange={handleChangeScale}
      />
    </>}
  >
    <Range
      quantity={quantity}
      visible={true}
      nBins={nbins}
      scale={scale}
      anchored={true}
      autorange={autorange}
      showinput={showinput}
      disableHistogram={false}
      aggId={id}
    />
  </Widget>
})

WidgetHistogram.propTypes = {
  id: PropTypes.string.isRequired,
  title: PropTypes.string,
  description: PropTypes.string,
  quantity: PropTypes.string,
  nbins: PropTypes.number,
  scale: PropTypes.string,
  autorange: PropTypes.bool,
  showinput: PropTypes.bool,
  className: PropTypes.string
}

/**
 * A dialog that is used to configure a scatter plot widget.
 */
export const WidgetHistogramEdit = React.memo((props) => {
    const {id, editing, visible} = props
    const { useSetWidget } = useSearchContext()
    const [settings, setSettings] = useState(props)
    const [errors, setErrors] = useState({})
    const setWidget = useSetWidget(id)
    const hasError = useMemo(() => {
      return Object.values(errors).some((d) => !!d) || !schemaWidgetHistogram.isValidSync(settings)
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

    const handleAccept = useCallback((key, value) => {
      try {
        schemaWidgetHistogram.validateSyncAt(key, {[key]: value})
      } catch (e) {
        handleError(key, e.message)
        return
      }
      setErrors(old => ({...old, [key]: undefined}))
      setSettings(old => ({...old, [key]: value}))
    }, [handleError, setSettings])

    const handleClose = useCallback(() => {
      setWidget(old => ({...old, editing: false}))
    }, [setWidget])

    const handleEditAccept = useCallback(() => {
      handleSubmit({...settings, editing: false, visible: true})
    }, [handleSubmit, settings])

    return <WidgetEditDialog
        id={id}
        open={editing}
        visible={visible}
        title="Edit histogram widget"
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
        <WidgetEditOption>
          <TextField
            select
            fullWidth
            label="Maximum number of bins"
            variant="filled"
            value={settings.nbins}
            onChange={(event) => { handleChange('nbins', event.target.value) }}
          >
            <MenuItem value={10}>10</MenuItem>
            <MenuItem value={20}>20</MenuItem>
            <MenuItem value={30}>30</MenuItem>
            <MenuItem value={40}>40</MenuItem>
            <MenuItem value={50}>50</MenuItem>
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
            control={<Checkbox checked={settings.autorange} onChange={(event, value) => handleChange('autorange', value)}/>}
            label={autorangeDescription}
          />
        </WidgetEditOption>
        <WidgetEditOption>
          <FormControlLabel
            control={<Checkbox checked={settings.showinput} onChange={(event, value) => handleChange('showinput', value)}/>}
            label='Show input fields'
          />
        </WidgetEditOption>
      </WidgetEditGroup>
    </WidgetEditDialog>
})

WidgetHistogramEdit.propTypes = {
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

export const schemaWidgetHistogram = schemaWidget.shape({
  quantity: string().required('Quantity is required.'),
  scale: string().required('Scale is required.'),
  nbins: number().integer().required(),
  autorange: bool(),
  showinput: bool()
})
