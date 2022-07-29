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
import { useRecoilValue } from 'recoil'
import MUIAddIcon from '@material-ui/icons/Add'
import AddCircleOutlineIcon from '@material-ui/icons/AddCircleOutline'
import CloseIcon from '@material-ui/icons/Close'
import HighlightOffIcon from '@material-ui/icons/HighlightOff'
import AddCircleIcon from '@material-ui/icons/AddCircle'
import CancelIcon from '@material-ui/icons/Cancel'
import { useSearchContext } from '../SearchContext'
import { Action } from '../../Actions'
import { guiState } from '../../GUIMenu'

/**
 * Toggle for hiding/showing a particular statistic.
 */
const StatisticsToggle = React.memo(({quantity, disabled, 'data-testid': testID}) => {
  const { useStatisticState } = useSearchContext()
  const iconSize = useRecoilValue(guiState('iconSize'))
  const icon = useRecoilValue(guiState('icon'))
  const [statistic, setStatistic] = useStatisticState(quantity)
  const hasStatistic = !!statistic

  let RemoveIcon, AddIcon
  if (icon === 'plain') {
    RemoveIcon = CloseIcon
    AddIcon = MUIAddIcon
  } else if (icon === 'outlined') {
    RemoveIcon = HighlightOffIcon
    AddIcon = AddCircleOutlineIcon
  } else if (icon === 'filled') {
    RemoveIcon = CancelIcon
    AddIcon = AddCircleIcon
  }

  return <Action
    tooltip={hasStatistic
      ? 'Remove statistics from the results panel.'
      : 'Display statistics for this filter in the results panel.'}
    disabled={disabled}
    data-testid={testID}
    onClick={() => {
      setStatistic(old => {
        return old ? undefined : {index: Date.now()}
      })
    }}
  >
    {hasStatistic
      ? <RemoveIcon fontSize={iconSize}/>
      : <AddIcon fontSize={iconSize}/>}
  </Action>
})

StatisticsToggle.propTypes = {
  quantity: PropTypes.string.isRequired,
  disabled: PropTypes.bool,
  'data-testid': PropTypes.string
}

export default StatisticsToggle
