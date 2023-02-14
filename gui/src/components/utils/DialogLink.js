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
import React, { useCallback, useState } from 'react'
import PropTypes from 'prop-types'
import { Button, Dialog, DialogActions, DialogContent, DialogTitle, Link } from '@material-ui/core'

/** A link that opens a simple help/explanation dialog */
const DialogLink = React.memo(function DialogLink(props) {
  const {title, text, linkProps, children, ...dialogProps} = props
  const [open, setOpen] = useState(false)

  const handleClick = useCallback(() => {
    setOpen(true)
  }, [setOpen])

  const handleClose = useCallback(() => {
    setOpen(false)
  }, [setOpen])

  return (
    <React.Fragment>
      <Dialog open={open} onClose={handleClose} {...dialogProps}>
        <DialogTitle>{title}</DialogTitle>
        <DialogContent>
          {children}
        </DialogContent>
        <DialogActions>
          <Button color="primary" onClick={handleClose}>close</Button>
        </DialogActions>
      </Dialog>
      <Link onClick={handleClick} href="#" {...linkProps}>
        {text}
      </Link>
    </React.Fragment>
  )
})
DialogLink.propTypes = {
  /** The dialog title. */
  title: PropTypes.string.isRequired,
  /** The text to render the link. */
  text: PropTypes.string.isRequired,
  /** The dialog contents. */
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ]),
  /** Additional optional dialog props. */
  linkProps: PropTypes.object
  /** Other props are forwarded to the dialog compotnent. */
}

export default DialogLink
