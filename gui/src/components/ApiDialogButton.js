import React, { useState } from 'react'
import PropTypes from 'prop-types'
import { IconButton, Dialog, DialogTitle, DialogContent, DialogActions, Button, Tooltip, Typography, makeStyles } from '@material-ui/core'
import CodeIcon from '@material-ui/icons/Code'
import ReactJson from 'react-json-view'
import Markdown from './Markdown'
import { CopyToClipboard } from 'react-copy-to-clipboard'
import ClipboardIcon from '@material-ui/icons/Assignment'

const useApiDialogStyles = makeStyles(theme => ({
  content: {
    paddingBottom: 0
  },
  json: {
    marginTop: theme.spacing(2),
    marginBottom: theme.spacing(2)
  },
  codeContainer: {
    display: 'flex',
    flexDirection: 'row',
    alignItems: 'flex-start'
  },
  code: {
    flexGrow: 1,
    marginRight: theme.spacing(1),
    overflow: 'hidden'
  },
  codeActions: {
    marginTop: theme.spacing(3)
  }
}))

export function ApiDialog({title, data, onClose, ...dialogProps}) {
  const classes = useApiDialogStyles()

  const renderCode = (title, code) => {
    return <React.Fragment>
      <Typography>{title}</Typography>
      <div className={classes.codeContainer}>
        <div className={classes.code}>
          <Markdown text={'```\n' + code + '\n```'} />
        </div>
        <div className={classes.codeActions}>
          <CopyToClipboard text={code} onCopy={() => null}>
            <Tooltip title="Copy to clipboard">
              <IconButton>
                <ClipboardIcon />
              </IconButton>
            </Tooltip>
          </CopyToClipboard>
        </div>
      </div>
    </React.Fragment>
  }

  return (
    <Dialog maxWidth="lg" fullWidth {...dialogProps}>
      <DialogTitle>{title || 'API Code'}</DialogTitle>

      <DialogContent classes={{root: classes.content}}>
        { data.code && data.code.repo_url &&
          renderCode(<span>URL to this query on the repository API:</span>, data.code.repo_url)}
        { data.code && data.code.curl &&
          renderCode(<span>Access the archive as JSON via <i>curl</i>:</span>, data.code.curl)}
        { data.code && data.code.python &&
          renderCode(<span>Access the archive in <i>python</i>:</span>, data.code.python)}
        { data.code && data.code.clientlib &&
          renderCode(<span>Access the archive with the <i>NOMAD client library</i>:</span>, data.code.clientlib)}

        <Typography>The repository API response as JSON:</Typography>
        <div className={classes.codeContainer}>
          <div className={classes.code}>
            <div className={classes.json}>
              <ReactJson
                src={data}
                enableClipboard={false}
                collapsed={2}
                displayObjectSize={false}
              />
            </div>
          </div>
          <div className={classes.codeActions}>
            <CopyToClipboard text={data} onCopy={() => null}>
              <Tooltip title="Copy to clipboard">
                <IconButton>
                  <ClipboardIcon />
                </IconButton>
              </Tooltip>
            </CopyToClipboard>
          </div>
        </div>
      </DialogContent>

      <DialogActions>
        <Button onClick={onClose}>
          Close
        </Button>
      </DialogActions>
    </Dialog>
  )
}
ApiDialog.propTypes = {
  data: PropTypes.any.isRequired,
  title: PropTypes.string,
  onClose: PropTypes.func
}

export default function ApiDialogButton({component, ...dialogProps}) {
  const [showDialog, setShowDialog] = useState(false)

  return (
    <React.Fragment>
      {component ? component({onClick: () => setShowDialog(true)}) : <Tooltip title="Show API code">
        <IconButton onClick={() => setShowDialog(true)}>
          <CodeIcon />
        </IconButton>
      </Tooltip>
      }
      <ApiDialog
        {...dialogProps} open={showDialog}
        onClose={() => setShowDialog(false)}
      />
    </React.Fragment>
  )
}
ApiDialogButton.propTypes = {
  data: PropTypes.any.isRequired,
  title: PropTypes.string,
  component: PropTypes.func
}
