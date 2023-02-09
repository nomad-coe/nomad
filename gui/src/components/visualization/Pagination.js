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
import { Box, Button, Tooltip } from '@material-ui/core'
import ArrowDownIcon from '@material-ui/icons/ArrowDropDown'
import ArrowUpIcon from '@material-ui/icons/ArrowDropUp'
import LoadingButton from '../buttons/LoadingButton'

/**
 * Displays simple pagination controls.
 */
const Pagination = React.memo(({
  showMore,
  showLess,
  disableMore,
  disableLess,
  loadingMore,
  onMore,
  onLess,
  variant,
  'data-testid': testID,
  ...rest
}) => {
  const moreIcon = {
    'standard': null,
    'down': <ArrowDownIcon />,
    'up': <ArrowUpIcon />
  }[variant]
  const lessIcon = {
    'standard': null,
    'down': <ArrowUpIcon />,
    'up': <ArrowDownIcon />
  }[variant]
  return <Box display="flex" {...rest} data-testid={testID}>
    {showMore && <Tooltip title={loadingMore ? 'Loading...' : disableMore ? 'No more values available' : ''}>
      <span>
        <LoadingButton
          size="small"
          onClick={onMore}
          loading={loadingMore}
          disabled={disableMore}
        >{moreIcon}Show more
        </LoadingButton>
      </span>
    </Tooltip>}
    {showLess && <Button size="small" disabled={disableLess}
      onClick={onLess}
    >{lessIcon}Show less
    </Button>}
  </Box>
})

Pagination.propTypes = {
  showMore: PropTypes.bool, // Show 'More' button
  showLess: PropTypes.bool, // Show 'Less' button
  disableMore: PropTypes.bool, // Disable more button
  disableLess: PropTypes.bool, // Disable less button
  loadingMore: PropTypes.bool, // Show loading indicator in 'More' button
  onMore: PropTypes.func, // 'More' button handler
  onLess: PropTypes.func, // 'Less' button handler
  /*
   * Visualization variant. 'standard' means no directional hints, 'up' and
   * 'down' will display directional hint icons in the buttons.
   */
  variant: PropTypes.oneOf(['standard', 'up', 'down']),
  'data-testid': PropTypes.string
}

Pagination.defaultProps = {
  variant: 'standard'
}

export default Pagination
