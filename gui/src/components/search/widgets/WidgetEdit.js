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
import React, {useCallback, useState} from 'react'
import PropTypes from 'prop-types'
import {
  List,
  ListItem,
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  makeStyles,
  withStyles,
  Typography,
  Accordion as MuiAccordion,
  AccordionDetails as MuiAccordionDetails,
  AccordionSummary as MuiAccordionSummary
} from '@material-ui/core'
import { useSearchContext } from '../SearchContext'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'

/**
 * A dialog that is used to configure a widget.
 */
export const WidgetEditDialog = React.memo(({id, title, open, visible, onAccept, onClose, error, children}) => {
    const { useRemoveWidget } = useSearchContext()
  const removeWidget = useRemoveWidget()
    const handleClose = useCallback(() => {
      // If the widget has not bee visualized, then closing the dialog deletes
      // the widget completely.
      onClose && onClose()
      if (!visible) {
        removeWidget(id)
      }
    }, [onClose, removeWidget, id, visible])

    const handleAccept = useCallback(() => {
      onAccept && onAccept()
    }, [onAccept])

    return <Dialog
      fullWidth={true}
      maxWidth="sm"
      open={open}
      onClose={handleClose}
    >
      <DialogTitle>{title || ''}</DialogTitle>
      <DialogContent>
        {children}
      </DialogContent>
      <DialogActions>
        <Button onClick={handleClose} color="primary">
          Cancel
        </Button>
        <Button disabled={!!error} onClick={handleAccept} color="primary">
          Done
        </Button>
      </DialogActions>
    </Dialog>
})

WidgetEditDialog.propTypes = {
  id: PropTypes.string,
  title: PropTypes.string,
  open: PropTypes.bool,
  visible: PropTypes.bool,
  onClose: PropTypes.func,
  onAccept: PropTypes.func,
  error: PropTypes.bool,
  children: PropTypes.node
}

/**
 * A group of options in an edit dialog.
 */

const Accordion = withStyles((theme) => ({
  root: {
    boxShadow: 'none',
    margin: 0,
    '&$expanded': {
      margin: 0
    }
  },
  expanded: {
    margin: 0
  }
}))(MuiAccordion)

const AccordionDetails = withStyles((theme) => ({
  root: {
    padding: theme.spacing(0)
  }
}))(MuiAccordionDetails)

const AccordionSummary = withStyles((theme) => ({
  root: {
    minHeight: 24,
    padding: 0,
    '&$expanded': {
      minHeight: 24
    }
  },
  content: {
    '&$expanded': {
      margin: '0px 0'
    }
  },
  expanded: {}
}))(MuiAccordionSummary)

const useEditGroupStyles = makeStyles((theme) => ({
  heading: {
    color: theme.palette.primary.main
  },
  list: {
    width: '100%'
  }
}))
export const WidgetEditGroup = React.memo(({title, children}) => {
  const styles = useEditGroupStyles()
  const [expanded, setExpanded] = useState(true)

  return <Accordion expanded={expanded} onChange={() => { setExpanded(old => !old) }}>
    <AccordionSummary
        expandIcon={<ExpandMoreIcon />}
        IconButtonProps={{
          size: "small"
        }}
    >
      <Typography variant="button" className={styles.heading}>{title}</Typography>
    </AccordionSummary>
    <AccordionDetails>
      <List dense className={styles.list}>
        {children}
      </List>
    </AccordionDetails>
  </Accordion>
})

WidgetEditGroup.propTypes = {
  title: PropTypes.string,
  children: PropTypes.node
}

/**
 * An option in an edit dialog.
 */
export const WidgetEditOption = React.memo(({children}) => {
  return <ListItem>
    {children}
  </ListItem>
})

WidgetEditOption.propTypes = {
  children: PropTypes.node
}
