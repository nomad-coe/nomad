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
import { makeStyles, alpha } from '@material-ui/core/styles'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ExpandLessIcon from '@material-ui/icons/ExpandLess'
import { useResizeDetector } from 'react-resize-detector'
import { useWindowSize } from '../../hooks'
import PropTypes from 'prop-types'
import clsx from 'clsx'

/**
 * Component for displaying a scrollable view.
 */
const useStyles = makeStyles(theme => {
  const startColor = alpha(theme.palette.grey[100], 0.0)
  const endColor = alpha(theme.palette.grey[300], 0.4)
  return {
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
      boxSizing: 'border-box'
    },
    scrollHint: {
      zIndex: 3,
      position: 'absolute',
      width: '100%',
      height: '1.5rem',
      background: alpha(theme.palette.background.default, 0.8),
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'center',
      pointerEvents: 'none'
    },
    scrollUp: {
      background: `linear-gradient(0deg, ${endColor} 0%, ${startColor} 100%)`,
      bottom: 0
    },
    scrollDown: {
      borderTop: `1px solid ${theme.palette.grey[300]}`,
      background: `linear-gradient(0deg, ${startColor} 0%, ${endColor} 100%)`,
      top: 0
    }
  }
})
const useScrollStyles = makeStyles((theme) => ({
  scrollBar: {
    overflowY: 'scroll',
    '&::-webkit-scrollbar': {
      display: 'none'
    },
    '-ms-overflow-style': 'none',
    scrollbarWidth: 'none'
  }
}))

const Scrollable = React.memo(({
  onBottom,
  scrollBar,
  hints,
  footer,
  className,
  classes,
  children,
  'data-testid': testID
}) => {
  const refRoot = useRef()
  const { height, ref } = useResizeDetector()
  const { windowHeight } = useWindowSize()
  const [isOverflowing, setIsOverflowing] = useState(false)
  const [isTop, setIsTop] = useState(true)
  const [isBottom, setIsBottom] = useState(false)
  const styles = useStyles({classes: classes})
  const scrollStyles = useScrollStyles({classes: {scrollBar: classes?.scrollBar}})

  // When the height of children changes or the window size changes, redetermine
  // whether to show scrollhints
  useEffect(() => {
    const element = refRoot?.current
    if (element) {
      const scrollHeight = element.scrollHeight - element.offsetHeight
      const overflowing = scrollHeight > 0
      setIsOverflowing(overflowing)
    }
  }, [height, windowHeight])

  // Attach event listener for checking if scrolled to top or bottom
  useEffect(() => {
    const element = refRoot?.current

    function updateScrollPosition() {
      const scrollHeight = element.scrollHeight - element.offsetHeight
      const bottom = Math.round(element.scrollTop) >= scrollHeight
      const top = Math.round(element.scrollTop) <= 0
      setIsTop(top || null)
      setIsBottom(bottom || null)
      if (bottom && onBottom) {
        onBottom()
      }
    }

    if (element) {
      element.addEventListener('scroll', updateScrollPosition, false)
      return function cleanup() {
        element.removeEventListener('scroll', updateScrollPosition, false)
      }
    }
  }, [onBottom])

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <div className={clsx(styles.containerOuter)}>
      {(hints && isOverflowing && !isTop) && <div className={clsx(styles.scrollHint, styles.scrollDown)}>
        <ExpandLessIcon color="action"/>
      </div>}
      <div ref={refRoot} className={clsx(styles.containerInner, !scrollBar && scrollStyles.scrollBar)}>
        <div ref={ref}>
          {children}
          {footer}
        </div>
      </div>
      {(hints && isOverflowing && !isBottom) && <div className={clsx(styles.scrollHint, styles.scrollUp)}>
        <ExpandMoreIcon color="action"/>
      </div>}
    </div>
  </div>
})

Scrollable.propTypes = {
  onBottom: PropTypes.func, // Callback that is fired when bottom reached
  scrollBar: PropTypes.bool, // Whether to show native scrollbar
  hints: PropTypes.bool, // Whether to show scroll hints
  fixed: PropTypes.bool, // Whether to show scroll hints
  footer: PropTypes.object, // A fixed component that is shown at the bottom
  className: PropTypes.string,
  classes: PropTypes.object,
  children: PropTypes.any,
  'data-testid': PropTypes.string
}
Scrollable.defaultProps = {
  scrollBar: false,
  hints: true
}

export default Scrollable
