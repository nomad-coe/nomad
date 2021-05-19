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
import React, { useEffect, useRef, useState } from 'react'
import { makeStyles, fade } from '@material-ui/core/styles'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ExpandLessIcon from '@material-ui/icons/ExpandLess'
import { useResizeDetector } from 'react-resize-detector'
import PropTypes from 'prop-types'
import clsx from 'clsx'

/**
 * Component for displaying a scrollable view.
 */

const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    height: '100%'
  },
  containerOuter: {
    width: '100%',
    height: '100%',
    position: 'relative'
  },
  containerInner: {
    width: '100%',
    height: '100%',
    overflowY: 'scroll',
    '&::-webkit-scrollbar': {
      display: 'none'
    },
    '-ms-overflow-style': 'none',
    scrollbarWidth: 'none',
    paddingTop: theme.spacing(1.5),
    paddingBottom: theme.spacing(1.5),
    boxSizing: 'border-box'
  },
  scrollHint: {
    zIndex: 3,
    position: 'absolute',
    width: '100%',
    height: '2rem',
    backgroundColor: fade(theme.palette.background.default, 0.8),
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center'
  },
  scrollUp: {
    bottom: 0
  },
  scrollDown: {
    top: 0
  }
}))
const Scrollable = React.memo(({
  className,
  classes,
  children,
  'data-testid': testID
}) => {
  const refRoot = useRef()
  const { height, ref } = useResizeDetector()
  const [isOverflowing, setIsOverflowing] = useState(false)
  const [isTop, setIsTop] = useState(true)
  const [isBottom, setIsBottom] = useState(false)
  const styles = useStyles({classes: classes})

  // When the height of children changes, redetermine whether to show
  // scrollhints
  useEffect(() => {
    const element = refRoot?.current
    if (element) {
      const scrollHeight = element.scrollHeight - element.offsetHeight
      const overflowing = scrollHeight > 0
      setIsOverflowing(overflowing)
    }
  }, [height])

  // Attach event listener for checking if scrolled to top or bottom
  useEffect(() => {
    const element = refRoot?.current

    function updateScrollPosition() {
      const scrollHeight = element.scrollHeight - element.offsetHeight
      const bottom = element.scrollTop === scrollHeight
      const top = element.scrollTop === 0
      setIsTop(top || null)
      setIsBottom(bottom || null)
    }

    if (element) {
      element.addEventListener('scroll', updateScrollPosition, false)
      return function cleanup() {
        element.removeEventListener('scroll', updateScrollPosition, false)
      }
    }
  }, [])

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <div className={clsx(styles.containerOuter)}>
      {(isOverflowing && !isTop) && <div className={clsx(styles.scrollHint, styles.scrollDown)}>
        <ExpandLessIcon color="action"/>
      </div>}
      <div ref={refRoot} className={clsx(styles.containerInner)}>
        <div ref={ref}>
          {children}
        </div>
      </div>
      {(isOverflowing && !isBottom) && <div className={clsx(styles.scrollHint, styles.scrollUp)}>
        <ExpandMoreIcon color="action"/>
      </div>}
    </div>
  </div>
})

Scrollable.propTypes = {
  className: PropTypes.string,
  classes: PropTypes.object,
  children: PropTypes.any,
  'data-testid': PropTypes.string
}

export default Scrollable
