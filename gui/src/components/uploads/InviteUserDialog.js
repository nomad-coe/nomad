import {
  CircularProgress,
  Dialog,
  DialogContent,
  DialogTitle, TextField, makeStyles
} from '@material-ui/core'
import Button from '@material-ui/core/Button'
import DialogActions from '@material-ui/core/DialogActions'
import DialogContentText from '@material-ui/core/DialogContentText'
import React, { useCallback, useState } from 'react'
import { useApi } from '../api'

export const useInviteUserDialogStyles = makeStyles(theme => ({
  button: {
    marginLeft: theme.spacing(1)
  },
  dialog: {
    width: '100%'
  },
  submitWrapper: {
    margin: theme.spacing(1),
    position: 'relative'
  },
  submitProgress: {
    position: 'absolute',
    top: '50%',
    left: '50%',
    marginTop: -12,
    marginLeft: -12
  }
}))

export const InviteUserDialog = React.memo(function InviteUserDialog(props) {
  const classes = useInviteUserDialogStyles()
  const [open, setOpen] = useState(false)
  const [submitting, setSubmitting] = useState(false)
  const [canSubmit, setCanSubmit] = useState(false)
  const [error, setError] = useState(null)
  const [data, setData] = useState({
    first_name: '',
    last_name: '',
    email: '',
    affiliation: ''
  })
  const { api } = useApi()

  const handleClose = useCallback((event, reason) => {
    if (reason !== 'backdropClick') {
      setOpen(false)
    }
  }, [setOpen])

  const handleSubmit = useCallback(() => {
    setSubmitting(true)

    api.inviteUser(data).then(() => {
      setSubmitting(false)
      setOpen(false)
    }).catch(error => {
      let message = '' + error
      try {
        message = JSON.parse(error.request.responseText).detail || '' + error
      } catch (e) { }
      setError(message)
      setSubmitting(false)
      setCanSubmit(false)
    })
  }, [data, setSubmitting, setOpen, setError, api])

  const handleChange = useCallback((key, value) => {
    const valid = value && !Object.keys(data).find(dataKey => !(key === dataKey || data[dataKey]))
    setData({ ...data, [key]: value })
    setCanSubmit(valid)
  }, [setData, data, setCanSubmit])

  const handleOpen = useCallback(() => {
    setOpen(true)
  }, [setOpen])

  const input = (key, label) => <TextField
    variant="filled"
    label={label}
    value={data[key]}
    onChange={event => handleChange(key, event.target.value)}
    margin="normal"
    fullWidth />
  return <React.Fragment>
    <Button className={classes.button}
      onClick={handleOpen}
      color="secondary" disabled={submitting}
    >
      Invite new user
    </Button>
    <Dialog
      classes={{ paper: classes.dialog }}
      open={open}
      onClose={handleClose} disableEscapeKeyDown>
      <DialogTitle>Invite a new user to NOMAD</DialogTitle>
      <DialogContent>
        <DialogContentText>
          If you want to add a user as co-author or share your data with someone that
          is not already a NOMAD user, you can invite this person here. We need just a few
          details about this person. After your invite, the new user will receive an
          Email that allows her to set a password and further details. Anyhow, you will
          be able to add the user as co-author or someone to share with immediately after the
          invite.
        </DialogContentText>
        {error && <DialogContentText color="error">
          {error}
        </DialogContentText>}
        {input('email', 'Email')}
        {input('first_name', 'First name')}
        {input('last_name', 'Last name')}
        {input('affiliation', 'Affiliation')}
      </DialogContent>
      <DialogActions>
        <Button onClick={handleClose} disabled={submitting}>
          Cancel
        </Button>
        <div className={classes.submitWrapper}>
          <Button onClick={handleSubmit} color="primary" disabled={!canSubmit}>
            Submit
          </Button>
          {submitting && <CircularProgress size={24} className={classes.submitProgress} />}
        </div>
      </DialogActions>
    </Dialog>
  </React.Fragment>
})
InviteUserDialog.propTypes = {}
