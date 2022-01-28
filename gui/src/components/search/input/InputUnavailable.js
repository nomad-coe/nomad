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
import { Typography, makeStyles } from '@material-ui/core'

/**
 * Indicates that no input options are available with the current query
 * settings.
 */
const useStyles = makeStyles(theme => ({
  root: {
    fontStyle: 'italic',
    color: theme.palette.text.disabled
  }
}))
export default function InputUnavailable({className}) {
  const styles = useStyles()
  return <Typography className={clsx(className, styles.root)}>
    No options available with current query.
  </Typography>
}

InputUnavailable.propTypes = {
  className: PropTypes.string
}
