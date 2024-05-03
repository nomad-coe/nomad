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
import React, { useCallback, useMemo } from 'react'
import clsx from 'clsx'
import PropTypes from 'prop-types'
import { cloneDeep, isNil } from 'lodash'
import { object, string, number } from 'yup'
import { makeStyles } from '@material-ui/core/styles'
import { Edit } from '@material-ui/icons'
import { Action } from '../../Actions'
import WidgetHeader from './WidgetHeader'

/**
 * A thin wrapper for displaying widgets.
 */
const useStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
    flexDirection: 'column'
  }
}))
export const Widget = React.memo(({
  id,
  quantity,
  title,
  description,
  onEdit,
  className,
  children,
  actions
}) => {
  const styles = useStyles()

  const handleEdit = useCallback(() => {
    onEdit && onEdit()
  }, [onEdit])

  const actionsFinal = useMemo(() => {
    return <>
      {!!actions && actions}
      <Action tooltip='Edit' onClick={handleEdit}>
        <Edit fontSize="small"/>
      </Action>
    </>
  }, [handleEdit, actions])

  return <div className={clsx(styles.root, className)}>
    <WidgetHeader
      id={id}
      quantity={quantity}
      label={title}
      description={description}
      actions={actionsFinal}
      anchored
    />
    {children}
  </div>
})

Widget.propTypes = {
  id: PropTypes.string.isRequired,
  quantity: PropTypes.string,
  title: PropTypes.string,
  description: PropTypes.string,
  onEdit: PropTypes.func,
  className: PropTypes.string,
  children: PropTypes.node,
  actions: PropTypes.node
}

export const schemaLayout = object({
  x: number().test(
    'is-integer-or-infinity',
    // eslint-disable-next-line no-template-curly-in-string
    '${path} is not a valid integer number',
    (value, context) => Number.isInteger(value) || value === Infinity || isNil(value)
  ),
  y: number().integer(),
  w: number().integer(),
  h: number().integer(),
  minW: number().integer(),
  minH: number().integer()
})
export const schemaWidget = object({
  type: string().required(),
  title: string(),
  layout: object({
    sm: schemaLayout,
    md: schemaLayout,
    lg: schemaLayout,
    xl: schemaLayout,
    xxl: schemaLayout
  }),
  id: string().strip(),
  editing: string().strip(),
  visible: string().strip()
})
export const schemaAxis = object({
  quantity: string().required(),
  unit: string().nullable(),
  title: string().nullable()
})
export const schemaAxisOptional = object({
  quantity: string().nullable(),
  unit: string().nullable(),
  title: string().nullable()
})
export const schemaMarkers = object({
  color: schemaAxisOptional
})

/**
 * Transforms a list of widgets into the internal object representation.
 * @param {array} widgetList Array of widget configs.
 */
export function getWidgetsObject(widgetList) {
  const widgets = cloneDeep(widgetList)
  return widgets
    ? Object.fromEntries(widgets.map((widget, index) => {
      return [
        index.toString(),
        {...widget, id: index.toString(), editing: false, visible: true}
      ]
    }))
    : {}
}
