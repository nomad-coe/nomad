
import React, { useEffect, useMemo, useState } from 'react'
import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  FormControlLabel,
  FormGroup,
  Switch
} from '@material-ui/core'
import PropTypes from 'prop-types'
import Markdown from '../Markdown'
import { matomo } from '../App'
import { useCookies } from 'react-cookie'
import { guiBase, consent } from '../../config'

/**
 * Consent form.
 */
function Consent({open, onAccept}) {
  const [cookies, setCookie] = useCookies()
  const [optOut, setOptOut] = useState(cookies['tracking-enabled'] === 'false')
  const cookieOptions = useMemo(() => ({
    expires: new Date(2147483647 * 1000),
    path: '/' + guiBase.split('/').slice(1).join('/')
  }), [])

  useEffect(() => {
    if (!optOut) {
      matomo.push(['setConsentGiven'])
    } else {
      matomo.push(['requireConsent'])
    }
  })

  // Write again on initial render to push forwards Safari's hard-coded 7 days
  // ITP window
  useEffect(() => {
    setCookie('terms-accepted', cookies['terms-accepted'], cookieOptions)
    setCookie('tracking-enabled', cookies['tracking-enabled'], cookieOptions)
  },
  // eslint-disable-next-line react-hooks/exhaustive-deps
  [])

  // Write cookies when form accepted.
  const handleAccept = () => {
    setCookie('terms-accepted', true, cookieOptions)
    setCookie('tracking-enabled', !optOut, cookieOptions)
    onAccept()
  }

  // When the form is opened, the cookies are by default rejected until the user
  // accepts the terms.
  useEffect(() => {
    if (open) {
      setCookie('terms-accepted', false, cookieOptions)
      setCookie('tracking-enabled', false, cookieOptions)
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [cookieOptions, open])

  return <Dialog
    disableBackdropClick
    disableEscapeKeyDown
    open={open}
  >
    <DialogTitle>Terms of Use</DialogTitle>
    <DialogContent>
      <Markdown>{consent}</Markdown>
      <FormGroup>
        <FormControlLabel
          control={<Switch
            checked={optOut}
            onChange={(e) => {
              setOptOut(old => !old)
            }}
            color="primary"
          />}
          label="Do not provide information about your use of NOMAD (opt-out)."
        />
      </FormGroup>
    </DialogContent>
    <DialogActions>
      <Button onClick={handleAccept} color="primary">
        Accept
      </Button>
    </DialogActions>
  </Dialog>
}

Consent.propTypes = {
  open: PropTypes.bool,
  onAccept: PropTypes.func
}

export default Consent
