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
import clsx from 'clsx'
import { makeStyles } from '@material-ui/core'

/**
 * Component for displaying text that should be truncated with an Ellipsis
 * either in front of the text or after the text.
 */
const useStyles = makeStyles(theme => ({
  root: {
    dislay: 'inline-block',
    whiteSpace: 'nowrap',
    overflow: 'hidden',
    textOverflow: 'ellipsis'
  },
  // This is a trick to make the ellipsis appear in front. Note the special sign
  // injected in front that makes the content flow form left-to-right.
  ellipsisFront: {
    direction: 'rtl',
    textAlign: 'left',
    '&::before': {
      content: '"‎"'
    }
  }
}))

const Ellipsis = React.memo(({front, children}) => {
  const styles = useStyles()
  return <div className={clsx(styles.root, front && styles.ellipsisFront)}>{children}</div>
})

Ellipsis.propTypes = {
  children: PropTypes.node,
  front: PropTypes.bool
}

export default Ellipsis
