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
import { makeStyles, Typography } from '@material-ui/core'

export const limitedWidthWidth = 1024

const useStyles = makeStyles(theme => ({
  fullWidth: {
    margin: theme.spacing(2)
  },
  limitedWidth: {
    width: '100%',
    maxWidth: limitedWidthWidth,
    margin: 'auto',
    marginTop: theme.spacing(2)
  }
}))

export default function Page({limitedWidth, loading, children}) {
  const classes = useStyles()

  return <div className={limitedWidth ? classes.limitedWidth : classes.fullWidth}>
    {loading && <Typography>loading ...</Typography>}
    {!loading && children}
  </div>
}

Page.propTypes = {
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ]),
  limitedWidth: PropTypes.bool,
  loading: PropTypes.bool
}
