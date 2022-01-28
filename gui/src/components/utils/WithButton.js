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

import React from 'react'
import PropTypes from 'prop-types'
import { IconButton, makeStyles, Tooltip } from '@material-ui/core'
import CopyToClipboard from 'react-copy-to-clipboard'
import ClipboardIcon from '@material-ui/icons/Assignment'

const useStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
    flexDirection: 'row',
    alignItems: 'center',
    '& div': {
      marginRight: theme.spacing(1)
    }
  }
}))

export default function WithButton({children, button, icon, tooltip, onClick, clipboard, ...rest}) {
  const classes = useStyles()

  if (!icon && clipboard) {
    icon = <ClipboardIcon style={{fontSize: 16}} />
    rest.size = rest.size || 'small'
  }

  if (!button) {
    button = <IconButton onClick={onClick} {...rest}>
      {icon}
    </IconButton>
  }

  if (tooltip) {
    button = <Tooltip title={tooltip}>
      {button}
    </Tooltip>
  }

  if (clipboard) {
    button = <CopyToClipboard
      text={clipboard} onCopy={() => null}
    >
      <div onClick={e => e.stopPropagation()}>
        {button}
      </div>
    </CopyToClipboard>
  }

  return <div className={classes.root}>
    <div>{children}</div>
    {button}
  </div>
}

WithButton.propTypes = {
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ]),
  button: PropTypes.node,
  icon: PropTypes.node,
  tooltip: PropTypes.string,
  onClick: PropTypes.func,
  clipboard: PropTypes.string
}
