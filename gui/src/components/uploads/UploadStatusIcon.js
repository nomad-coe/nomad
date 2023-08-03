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
import { Tooltip } from '@material-ui/core'
import PublishedIcon from '@material-ui/icons/Public'
import PrivateIcon from '@material-ui/icons/AccountCircle'
import SharedIcon from '@material-ui/icons/SupervisedUserCircle'
import NotVisibleIcon from '@material-ui/icons/VisibilityOff'

/**
 * Used to display the status of an Upload. Can work both with upload data and
 * entry data which have different data models.
 */
const UploadStatusIcon = React.memo(({data, user, ...props}) => {
  const coauthors = data?.coauthors || data?.authors?.map(user => user.user_id)
  const reviewers = data?.reviewers || data?.viewers?.map(user => user.user_id)
  const shared = data?.coauthors?.length > 0 || data?.reviewers?.length > 0 || data?.viewers?.length > 1
  const isMainAuthor = user && (data.main_author?.user_id === user.sub || data.main_author === user.sub)
  const isReviewer = user && reviewers?.find(user_id => user_id === user.sub)
  const isCoauthor = user && coauthors?.find(user_id => user_id === user.sub)
  let Icon = shared ? SharedIcon : PrivateIcon

  let tooltip, color, role
  if (!data) {
    tooltip = 'Upload status not available'
    color = 'action'
  } else if (data.published) {
    if (data.with_embargo) {
      color = 'error'
      if (isMainAuthor) {
        tooltip = "Published with embargo by you and only accessible by you, coauthors and reviewers"
      } else if (isCoauthor) {
        tooltip = "Published with embargo and accessible by you as a coauthor"
      } else if (isReviewer) {
        tooltip = "Published with embargo and accessible by you as a reviewer"
      } else {
        Icon = NotVisibleIcon
        if (user) {
          tooltip = "Published with embargo and not accessible by you"
        } else {
          tooltip = "Published with embargo and might become accessible after login"
        }
      }
    } else {
        role = 'published-upload-icon'
        tooltip = "Published and accessible by everyone"
        Icon = PublishedIcon
        color = 'primary'
    }
  } else {
      if (isMainAuthor) {
        tooltip = "Unpublished, only accessible by you, coauthors and reviewers"
      } else if (isCoauthor) {
        tooltip = "Unpublished, accessible by you as a coauthor"
      } else if (isReviewer) {
        tooltip = "Unpublished, accessible by you as a reviewer"
      } else {
        tooltip = "Unpublished"
      }
      color = 'error'
  }

  return <Tooltip title={tooltip}>
      <Icon color={color} role={role} {...props} />
  </Tooltip>
})
UploadStatusIcon.propTypes = {
  data: PropTypes.object, // The upload/entry data
  user: PropTypes.object // The user object from API
}

export default UploadStatusIcon
