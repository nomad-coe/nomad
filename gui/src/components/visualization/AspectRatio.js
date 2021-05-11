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
import { makeStyles } from '@material-ui/core/styles'
import PropTypes from 'prop-types'
import clsx from 'clsx'

/**
 * Component that wraps it's children in a container that can be provided a
 * fixed aspect ratio. The width of the root component is used to determine the
 * height.
 */
export default function AspectRatio({className, classes, children, aspectRatio}) {
  const useStyles = makeStyles((theme) => {
    return {
      root: {
      },
      containerOuter: {
        width: '100%',
        height: 0,
        paddingBottom: `${100 / aspectRatio}%`, // CSS hack for fixed aspect ratio
        position: 'relative',
        boxSizing: 'border-box'
      },
      containerInner: {
        position: 'absolute',
        top: 0,
        right: 0,
        bottom: 0,
        left: 0,
        display: 'flex',
        flexDirection: 'column',
        boxSizing: 'border-box'
      }
    }
  })
  const style = useStyles({classes: classes})

  return aspectRatio
    ? <div className={clsx(style.root, className)}>
      <div className={style.containerOuter}>
        <div className={style.containerInner}>
          {children}
        </div>
      </div>
    </div>
    : <div className={clsx(style.root, className)}>
      {children}
    </div>
}

AspectRatio.propTypes = {
  aspectRatio: PropTypes.number, // width/height
  children: PropTypes.any,
  className: PropTypes.string,
  classes: PropTypes.object
}
