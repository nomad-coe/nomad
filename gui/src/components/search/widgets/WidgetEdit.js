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
import {
  List,
  ListItem,
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  makeStyles,
  MenuItem,
  TextField,
  Typography
} from '@material-ui/core'
import { useSearchContext } from '../SearchContext'

/**
 * A dialog that is used to configure a widget.
 */
export const useStyles = makeStyles(theme => ({
  width: {
    maxWidth: "700px"
  }
}))
export const WidgetEditDialog = React.memo(({id, title, open, visible, onAccept, onClose, error, children}) => {
  const styles = useStyles()
  const { useRemoveWidget } = useSearchContext()
  const removeWidget = useRemoveWidget()

  const handleClose = useCallback((event, reason) => {
    // Do not close dialog on backdrop click: user may lose lot of filled info
    // accidentally
    if (reason === 'backdropClick') {
      return
    }

    // If the widget has not bee visualized, then closing the dialog deletes
    // the widget completely.
    onClose && onClose()
    if (!visible) {
      removeWidget(id)
    }
  }, [id, onClose, removeWidget, visible])

  const handleAccept = useCallback(() => {
    onAccept && onAccept()
  }, [onAccept])

  return <Dialog
    fullWidth={true}
    maxWidth="sm"
    open={open}
    classes={{paperWidthSm: styles.width}}
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

  return <>
    <Typography variant="button" className={styles.heading}>{title}</Typography>
    <List dense className={styles.list}>
      {children}
    </List>
  </>
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

/**
 * Select (=dropdown) component for widget edit.
 */
export const WidgetEditSelect = React.memo(({label, disabled, options, value, onChange}) => {
  return <TextField
      select
      fullWidth
      label={label}
      variant="filled"
      value={value || ''}
      onChange={onChange}
      disabled={disabled}
    >
      {Object.entries(options).map(([key, value]) =>
        <MenuItem value={value} key={key}>{key}</MenuItem>
      )}
    </TextField>
})

WidgetEditSelect.propTypes = {
  label: PropTypes.string,
  disabled: PropTypes.bool,
  options: PropTypes.object,
  value: PropTypes.oneOfType([PropTypes.string, PropTypes.number]),
  onChange: PropTypes.func
}
