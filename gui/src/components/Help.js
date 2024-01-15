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
import React, { forwardRef, useCallback, useState } from 'react'
import { Button, IconButton, Dialog, DialogTitle, DialogContent, DialogActions } from '@material-ui/core'
import Markdown from './Markdown'
import PropTypes from 'prop-types'
import HelpIcon from '@material-ui/icons/Help'

export const HelpButton = forwardRef(({heading, content, text, maxWidth, IconProps, children, ...IconButtonProps}, ref) => {
  const [open, setOpen] = useState(false)
  const handleToggleOpen = useCallback(() => {
    setOpen(old => !old)
    IconButtonProps?.onClick?.()
  }, [IconButtonProps])

  return <>
    <IconButton {...IconButtonProps} onClick={handleToggleOpen} ref={ref}>
      {children || <HelpIcon {...IconProps}/>}
    </IconButton>
    <HelpDialog heading={heading} content={content} text={text} open={open} onClose={handleToggleOpen} maxWidth={maxWidth} />
  </>
})

HelpButton.propTypes = {
  heading: PropTypes.string,
  content: PropTypes.string,
  text: PropTypes.string,
  maxWidth: PropTypes.string,
  IconProps: PropTypes.object,
  children: PropTypes.node
}

const HelpDialog = ({heading, content, text, maxWidth, open, onClose}) => {
  return <Dialog
    maxWidth={maxWidth}
    onClose={onClose}
    open={open}
  >
    <DialogTitle>{heading || 'Help'}</DialogTitle>
    <DialogContent>
      <Markdown text={text}>{content}</Markdown>
    </DialogContent>
    <DialogActions>
      <Button onClick={onClose} color="primary">
        Close
      </Button>
    </DialogActions>
  </Dialog>
}

HelpDialog.propTypes = {
  heading: PropTypes.string,
  content: PropTypes.string,
  text: PropTypes.string,
  maxWidth: PropTypes.string,
  open: PropTypes.bool,
  onClose: PropTypes.func
}
