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
import { makeStyles, Typography, useTheme } from '@material-ui/core'
import AutoSizer from 'react-virtualized-auto-sizer'

export const limitedWidthWidth = 1024

const useStyles = makeStyles(theme => ({
  fullWidth: {
    marginLeft: theme.spacing(2),
    marginRight: theme.spacing(2)
  },
  limitedWidth: {
    width: '100%',
    maxWidth: limitedWidthWidth + theme.spacing(2),
    paddingLeft: theme.spacing(2),
    paddingRight: theme.spacing(2),
    marginLeft: 'auto',
    marginRight: 'auto'
  },
  limitedHeight: {
    height: '100%'
  },
  fullHeight: {
    marginTop: theme.spacing(2),
    marginBottom: theme.spacing(2)
  }
}))

export default function Page({limitedWidth, limitedHeight, loading, children, ...moreProps}) {
  const classes = useStyles()
  const theme = useTheme()

  const pageWithWidth = <div
    className={clsx({
      [classes.fullWidth]: !limitedWidth,
      [classes.limitedWidth]: limitedWidth,
      [classes.limitedHeight]: limitedHeight
    })}
    {...moreProps}
  >
    {loading && <Typography>loading ...</Typography>}
    {!loading && children}
  </div>

  if (limitedHeight) {
    return <div {...moreProps} className={classes.limitedHeight}>
      <AutoSizer disableWidth>
        {({height}) => <div style={{height: height}}>
          <div style={{height: theme.spacing(2)}}/>
          <div style={{height: height - theme.spacing(4)}}>
            {pageWithWidth}
          </div>
        </div>}
      </AutoSizer>
    </div>
  } else {
    return <div className={classes.fullHeight} {...moreProps}>
      {pageWithWidth}
    </div>
  }
}

Page.propTypes = {
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ]),
  /** Limits the width of the page content to a standard width. */
  limitedWidth: PropTypes.bool,
  /** Limits the height of the page content to fit into the current view port. */
  limitedHeight: PropTypes.bool,
  loading: PropTypes.bool
}
