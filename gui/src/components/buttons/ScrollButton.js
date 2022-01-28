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
import React, { useCallback, useContext } from 'react'
import PropTypes from 'prop-types'
import { IconButton } from '@material-ui/core'
import { ScrollContext } from '../nav/Navigation'

/**
 * Button used for manually scrolling the current page.
 */
export default function ScrollButton({scrollAmount, relative, instant, onClick, children, ...other}) {
  const scrollContext = useContext(ScrollContext)

  // Handle click action
  const handleClick = useCallback(() => {
    const container = scrollContext.scrollParentRef.current
    const options = {
      top: scrollAmount,
      behavior: instant ? 'instant' : 'smooth'
    }
    relative
      ? container.scrollBy(options)
      : container.scrollTo(options)
    onClick && onClick()
  }, [scrollContext, scrollAmount, instant, relative, onClick])

  return <IconButton
    onClick={handleClick}
    {...other}
  >
    {children}
  </IconButton>
}

ScrollButton.propTypes = {
  /* Defines the scroll amount. */
  scrollAmount: PropTypes.number.isRequired,
  /* Whether the scrollAmount is relative to current position. If false,
   * scrollAmount is absolute. */
  relative: PropTypes.bool,
  /* Use to disable smoothing. */
  instant: PropTypes.bool,
  /* Callback function for clicks. */
  onClick: PropTypes.func,
  children: PropTypes.node
}
