import React, { useState } from 'react'
import PropTypes from 'prop-types'
import { IconButton, Dialog, DialogTitle, DialogContent, DialogActions, Button, Tooltip } from '@material-ui/core'
import CodeIcon from '@material-ui/icons/Code'
import ReactJson from 'react-json-view'
import { CopyToClipboard } from 'react-copy-to-clipboard'

export default function CodeDialogButton({component, buttonProps, dialogProps, children, tooltip, title, actions, icon, ...moreProps}) {
  const [showDialog, setShowDialog] = useState(false)

  let button
  if (component) {
    button = component({onClick: () => setShowDialog(true)})
  } else {
    button = <Tooltip title={tooltip || 'Show code'}>
      <IconButton onClick={() => setShowDialog(true)} {...buttonProps}>
        {icon || <CodeIcon />}
      </IconButton>
    </Tooltip>
  }

  return (
    <React.Fragment>
      {button}
      <Dialog maxWidth="lg" fullWidth {...dialogProps} {...moreProps} open={showDialog}>
        <DialogTitle>{title || 'Code'}</DialogTitle>
        <DialogContent>
          {children}
        </DialogContent>
        <DialogActions>
          {actions}
          <Button onClick={() => setShowDialog(false)}>
            Close
          </Button>
        </DialogActions>
      </Dialog>
    </React.Fragment>
  )
}
CodeDialogButton.propTypes = {
  title: PropTypes.string,
  tooltip: PropTypes.string,
  component: PropTypes.func,
  buttonProps: PropTypes.object,
  dialogProps: PropTypes.object,
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ]).isRequired,
  actions: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ]),
  icon: PropTypes.node
}

export function JsonCodeDialogButton({json, jsonProps, ...props}) {
  const jsonToRender = typeof json === 'object' ? json : {data: json}
  // serialized json has "private" values that might contain circles removed
  const code = JSON.stringify(json, (name, value) => name.startsWith('_') ? undefined : value, 2)
  const actions = <CopyToClipboard text={code} onCopy={() => null}>
    <Button>Copy to clipboard</Button>
  </CopyToClipboard>

  return <CodeDialogButton actions={actions} {...props}>
    <ReactJson
      src={jsonToRender}
      enableClipboard={false}
      collapsed={2}
      displayObjectSize={false}
      {...jsonProps}
    />
  </CodeDialogButton>
}
JsonCodeDialogButton.propTypes = {
  json: PropTypes.any,
  jsonProps: PropTypes.object
}
