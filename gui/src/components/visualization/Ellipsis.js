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
import React, { useEffect, useState, useCallback } from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { makeStyles, Tooltip } from '@material-ui/core'
import { useResizeDetector } from 'react-resize-detector'

/**
 * Component for displaying text that should be truncated with an Ellipsis
 * either in front of the text or after the text. Can lso show the full text on
 * hover as a tooltip.
 */
const useStyles = makeStyles(theme => ({
  root: {
    display: 'block',
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

const Ellipsis = React.memo(({front, children, tooltip}) => {
  const styles = useStyles()
  const {width, ref} = useResizeDetector()
  const [showTooltip, setShowTooltip] = useState(false)

  // Used to determine whether to show the tooltip or not
  const compareSize = useCallback(() => {
    if (!ref?.current) return
    const compare =
      ref.current.scrollWidth > ref.current.clientWidth
    setShowTooltip(compare)
  }, [ref])

  // Whenever the component width changes, check if tooltip should be shown
  useEffect(() => { compareSize() }, [width, compareSize])

  // Define state and function to update the value
  return <Tooltip
      title={showTooltip ? tooltip : ''}
      disableHoverListener={!showTooltip}
    >
    <span ref={ref} className={clsx(styles.root, front && styles.ellipsisFront)}>{children}</span>
  </Tooltip>
})

Ellipsis.propTypes = {
  children: PropTypes.node,
  tooltip: PropTypes.string,
  front: PropTypes.bool
}

Ellipsis.defaultProps = {
  tooltip: ''
}

export default Ellipsis
