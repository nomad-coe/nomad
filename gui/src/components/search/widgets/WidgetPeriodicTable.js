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
import { string } from 'yup'
import { TextField, MenuItem } from '@material-ui/core'
import { useSearchContext } from '../SearchContext'
import { Widget, schemaWidget } from './Widget'
import { ActionSelect } from '../../Actions'
import { WidgetEditDialog, WidgetEditGroup, WidgetEditOption } from './WidgetEdit'
import { InputTextField } from '../input/InputText'
import { PeriodicTable } from '../input/InputPeriodicTable'
import { scales } from '../../plotting/common'

/**
 * Displays a periodic table as a widget.
 */
export const WidgetPeriodicTable = React.memo((
{
  id,
  title,
  description,
  quantity,
  scale,
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
    actions={
      <ActionSelect
        value={scale}
        options={Object.keys(scales)}
        tooltip="Statistics scaling"
        onChange={handleChangeScale}
      />
    }
  >
    <PeriodicTable
      quantity={quantity}
      scale={scale}
      anchored={true}
      disableStatistics={true}
      visible={true}
      // Can't use the same aggregation identifier as the periodic table in the
      // filter menu: due to rendering order the aggregation may otherwise get
      // turned off when menu is not open.
      aggId={'widget'}
    />
  </Widget>
})

WidgetPeriodicTable.propTypes = {
  id: PropTypes.string.isRequired,
  title: PropTypes.string,
  description: PropTypes.string,
  quantity: PropTypes.string,
  scale: PropTypes.string,
  className: PropTypes.string
}

/**
 * A dialog that is used to configure a scatter plot widget.
 */
export const WidgetPeriodicTableEdit = React.memo((props) => {
    const {id, editing, visible} = props
    const { useSetWidget } = useSearchContext()
    const [settings, setSettings] = useState(props)
    const [errors, setErrors] = useState({})
    const setWidget = useSetWidget(id)
    const hasError = useMemo(() => {
      return Object.values(errors).some((d) => !!d) || !schemaWidgetPeriodicTable.isValidSync(settings)
    }, [errors, settings])

    const handleSubmit = useCallback((settings) => {
      setWidget(old => ({...old, ...settings}))
    }, [setWidget])

    const handleChange = useCallback((key, value) => {
      setSettings(old => ({...old, [key]: value}))
    }, [setSettings])

    const handleClose = useCallback(() => {
      setWidget(old => ({...old, editing: false}))
    }, [setWidget])

    const handleError = useCallback((key, value) => {
      setErrors(old => ({...old, [key]: value}))
    }, [setErrors])

    const handleAccept = useCallback((key, value) => {
      try {
        schemaWidgetPeriodicTable.validateSyncAt(key, {[key]: value})
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
        title="Edit periodic table widget"
        onClose={handleClose}
        onAccept={handleEditAccept}
        error={hasError}
      >
      <WidgetEditGroup title="Elements">
        <WidgetEditOption>
          <TextField
            select
            fullWidth
            label="Statistics scaling"
            variant="filled"
            value={settings.scale}
            onChange={(event) => {
              handleChange('scale', event.target.value)
              handleAccept('scale', event.target.value)
            }}
          >
            {Object.keys(scales).map((key) =>
              <MenuItem value={key} key={key}>{key}</MenuItem>
            )}
          </TextField>
        </WidgetEditOption>
      </WidgetEditGroup>
      <WidgetEditGroup title="General">
        <WidgetEditOption>
          <InputTextField
            label="title"
            fullWidth
            value={settings?.title}
            onChange={(event) => handleChange('title', event.target.value)}
          />
        </WidgetEditOption>
      </WidgetEditGroup>
    </WidgetEditDialog>
})

WidgetPeriodicTableEdit.propTypes = {
  id: PropTypes.string.isRequired,
  editing: PropTypes.bool,
  visible: PropTypes.bool,
  quantity: PropTypes.string,
  scale: PropTypes.string,
  onClose: PropTypes.func
}

export const schemaWidgetPeriodicTable = schemaWidget.shape({
  quantity: string().required('Quantity is required.'),
  scale: string().required('Scale is required.')
})
