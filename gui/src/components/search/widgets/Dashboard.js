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
import React, {useCallback, useState, useEffect} from 'react'
import PropTypes from 'prop-types'
import {
  Accordion,
  AccordionSummary,
  AccordionDetails,
  Typography,
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  withStyles
} from '@material-ui/core'
import { isArray, cloneDeep } from 'lodash'
import CodeIcon from '@material-ui/icons/Code'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ExpandLessIcon from '@material-ui/icons/ExpandLess'
import AddCircleIcon from '@material-ui/icons/AddCircle'
import ReplayIcon from '@material-ui/icons/Replay'
import WidgetGrid from './WidgetGrid'
import { Actions, Action } from '../../Actions'
import { useSearchContext } from '../SearchContext'
import { WidgetScatterPlotEdit, schemaWidgetScatterPlot } from './WidgetScatterPlotEdit'
import { WidgetHistogramEdit, schemaWidgetHistogram } from './WidgetHistogram'
import { WidgetTermsEdit, schemaWidgetTerms } from './WidgetTerms'
import { WidgetPeriodicTableEdit, schemaWidgetPeriodicTable } from './WidgetPeriodicTable'
import InputConfig from '../input/InputConfig'
import { cleanse } from '../../../utils'
import Markdown from '../../Markdown'
import { getWidgetsObject } from './Widget'
import { ContentButton } from '../../buttons/SourceDialogButton'

const DashboardAccordion = withStyles({
  root: {
    backgroundColor: 'transparent'
  },
  expanded: {}
})(Accordion)

const summaryMinHeight = 40
const DashboardAccordionSummary = withStyles({
  root: {
    padding: 0,
    margin: 0,
    minHeight: summaryMinHeight,
    '&$expanded': {
      minHeight: summaryMinHeight
    }
  },
  content: {
    margin: 0,
    '&$expanded': {
      margin: '0'
    }
  },
  expanded: {}
})(AccordionSummary)

const DashboardAccordionDetails = withStyles((theme) => ({
  root: {
    padding: 0,
    display: 'block'
  }
}))(AccordionDetails)

/**
 * A collapsible area that shows the currently active widgets in a grid.
 */
const Dashboard = React.memo(() => {
  const [expanded, setExpanded] = useState(true)
  const { useWidgetsValue, useAddWidget, useResetWidgets } = useSearchContext()
  const resetWidgets = useResetWidgets()
  const widgets = useWidgetsValue()
  const addWidget = useAddWidget()

  // Adds a new scatter plot widget
  const handleScatterplot = useCallback(() => {
    const index = Date.now()
    const id = index.toString()
    const layout = {x: Infinity, y: 0, w: 9, h: 6}
    const value = {
      id: id,
      editing: true,
      visible: false,
      layout: {
        sm: {...layout},
        md: {...layout},
        lg: {...layout},
        xl: {...layout},
        xxl: {...layout}
      },
      size: 1000,
      autorange: true,
      type: 'scatterplot'
    }
    addWidget(id, value)
  }, [addWidget])

  // Adds a new periodic table widget
  const handlePeriodicTable = useCallback(() => {
    const index = Date.now()
    const id = index.toString()
    const layout = {x: Infinity, y: 0, w: 12, h: 9}
    const value = {
      id: id,
      editing: true,
      visible: false,
      layout: {
        sm: {...layout},
        md: {...layout},
        lg: {...layout},
        xl: {...layout},
        xxl: {...layout}
      },
      scale: 'linear',
      quantity: 'results.material.elements',
      type: 'periodictable'
    }
    addWidget(id, value)
  }, [addWidget])

  // Adds a new histogram widget
  const handleHistogram = useCallback(() => {
    const index = Date.now()
    const id = index.toString()
    const layout = {x: Infinity, y: 0, w: 8, h: 3}
    const value = {
      id: id,
      editing: true,
      visible: false,
      layout: {
        sm: {...layout},
        md: {...layout},
        lg: {...layout},
        xl: {...layout},
        xxl: {...layout}
      },
      autorange: false,
      nbins: 30,
      scale: 'linear',
      type: 'histogram'
    }
    addWidget(id, value)
  }, [addWidget])

  // Adds a new histogram widget
  const handleTerms = useCallback(() => {
    const index = Date.now()
    const id = index.toString()
    const layout = {x: Infinity, y: 0, w: 6, h: 9}
    const value = {
      id: id,
      editing: true,
      visible: false,
      layout: {
        sm: {...layout},
        md: {...layout},
        lg: {...layout},
        xl: {...layout},
        xxl: {...layout}
      },
      scale: 'linear',
      type: 'terms'
    }
    addWidget(id, value)
  }, [addWidget])

  const handleExpand = useCallback(() => {
    setExpanded(old => !old)
  }, [])

  const handleReset = useCallback(() => {
    resetWidgets()
  }, [resetWidgets])

  return <DashboardAccordion expanded={expanded} elevation={0}>
    <DashboardAccordionSummary
      aria-controls="panel1a-content"
      id="panel1a-header"
    >
      <Actions justifyContent="flex-start">
        <DashboardAction
          title="Terms"
          onClick={handleTerms}
          tooltip="Add a terms widget to the dashboard"
        />
        <DashboardAction
          title="Histogram"
          onClick={handleHistogram}
          tooltip="Add a histogram widget to the dashboard"
        />
        <DashboardAction
          title="Scatter plot"
          onClick={handleScatterplot}
          tooltip="Add a scatter plot widget to the dashboard"
        />
        <DashboardAction
          title="Periodic table"
          onClick={handlePeriodicTable}
          tooltip="Add a periodic table widget to the dashboard"
        />
        <Action
          onClick={handleReset}
          tooltip="Reset the default dashboard configuration"
        >
            <ReplayIcon fontSize="small"/>
        </Action>
        <Action
          tooltip=""
          ButtonComponent={DashboardExportButton}
          ButtonProps={{
            tooltip: "Dashboard export/import",
            title: "Dashboard export/import",
            DialogProps: {
              maxWidth: "md",
              fullWidth: true
            },
            ButtonProps: {
              size: "small"
            }
          }}>
            <CodeIcon fontSize="small"/>
        </Action>
        <Action
          onClick={handleExpand}
          tooltip={expanded ? 'Hide dashboard' : 'Show dashboard'}
        >
            {expanded ? <ExpandLessIcon fontSize="small"/> : <ExpandMoreIcon fontSize="small"/>}
        </Action>
      </Actions>
    </DashboardAccordionSummary>
    <DashboardAccordionDetails>
      <WidgetGrid />
    </DashboardAccordionDetails>
    {widgets && Object.entries(widgets).map(([id, value]) => {
      if (!value.editing) return null
      const comp = {
        scatterplot: <WidgetScatterPlotEdit key={id} widget={value}/>,
        periodictable: <WidgetPeriodicTableEdit key={id} {...value}/>,
        histogram: <WidgetHistogramEdit key={id} {...value}/>,
        terms: <WidgetTermsEdit key={id} {...value}/>
      }[value.type]
      return comp || null
    })}
  </DashboardAccordion>
})

Dashboard.propTypes = {
}

const buttonProps = {
  variant: 'outlined',
  startIcon: <AddCircleIcon color="primary"/>
}
const DashboardAction = React.memo(({title, ...rest}) => {
  return <Action {...rest} ButtonComponent={Button} ButtonProps={buttonProps}>
    <Typography variant="caption">{title}</Typography>
  </Action>
})
DashboardAction.propTypes = {
  title: PropTypes.string,
  children: PropTypes.node
}

const schemas = {
  scatterplot: schemaWidgetScatterPlot,
  periodictable: schemaWidgetPeriodicTable,
  histogram: schemaWidgetHistogram,
  terms: schemaWidgetTerms
}

/**
 * A button that displays a dialog that can be used to modify the current
 * dashboard setup.
 */
export const DashboardExportButton = React.memo(({tooltip, title, DialogProps, ButtonProps}) => {
  const {useWidgetsState} = useSearchContext()
  const [widgets, setWidgets] = useWidgetsState()
  const [widgetExport, setWidgetExport] = useState()
  const [widgetImport, setWidgetImport] = useState()
  const [widgetImportError, setWidgetImportError] = useState(null)
  const [open, setOpen] = useState(false)

  // The widget export data. Only reacts if the menu is open.
  useEffect(() => {
    if (!open) return
    // Work on a copy: data in Recoil is unmutable
    const widgetsExport = cloneDeep(widgets)

    const exp = Object
      .values(widgetsExport)
      .map((widget) => {
        const schema = schemas[widget.type]
        const casted = schema?.cast(widget, {stripUnknown: true})

        // Remove undefined values. YUP cannot do this, and the YAML
        // serialization will otherwise include these.
        cleanse(casted)

        // Perform custom sort: type first, layout last
        const sorted = {}
        if (casted['type']) sorted['type'] = casted['type']
        const fields = new Set(['type', 'layout'])
        for (const [key, value] of Object.entries(casted)) {
          if (!fields.has(key)) {
            sorted[key] = value
          }
        }
        if (casted['layout']) sorted['layout'] = casted['layout']
        return sorted
      })
      setWidgetExport(exp)
  }, [open, widgets])

  // Validate new widget config
  const handleChange = useCallback((data) => {
    if (!isArray(data)) return
    let error
    for (const widget of data) {
      const type = widget.type
      const schema = schemas[type]
      if (!schema) return
      try {
        schema.validateSync(widget)
      } catch (e) {
        setWidgetImportError(e.message)
        error = e
        break
      }
    }
    if (!error) {
      setWidgetImport(data)
      setWidgetImportError(null)
    }
  }, [])

  // Save the dashboard setup
  const handleSave = useCallback(() => {
    if (widgetImport) {
      if (widgetImportError) return
      const widgetsInternal = getWidgetsObject(widgetImport)
      setWidgets(widgetsInternal)
    }
    setOpen(false)
  }, [widgetImportError, widgetImport, setWidgets])

  return <ContentButton
    tooltip={tooltip}
    buttonContent={<CodeIcon fontSize="small"/>}
    ButtonProps={{...ButtonProps, onClick: () => { setOpen(true) }}}
  >
    <Dialog {...DialogProps} open={open}>
      {title && <DialogTitle>{title}</DialogTitle>}
      <DialogContent>
        <Markdown text="The current dashboard configuration in YAML. You can modify it here directly and save it to update the dashboard."/>
        <InputConfig
          data={widgetExport}
          format='YAML'
          maxRows={30}
          minRows={30}
          onChange={handleChange}
          error={widgetImportError}
          onError={setWidgetImportError}
        />
      </DialogContent>
      <DialogActions>
        <Button onClick={() => setOpen(false)}>
          Cancel
        </Button>
        <Button disabled={!!widgetImportError} onClick={handleSave}>
          Save
        </Button>
      </DialogActions>
    </Dialog>
  </ContentButton>
})
DashboardExportButton.propTypes = {
  label: PropTypes.string,
  title: PropTypes.string,
  tooltip: PropTypes.string,
  DialogProps: PropTypes.object,
  ButtonProps: PropTypes.object,
  children: PropTypes.node
}

export default Dashboard
