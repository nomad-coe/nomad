import GroupIcon from '@material-ui/icons/Group'
import PersonIcon from '@material-ui/icons/Person'
import _ from 'lodash'
import PropTypes from 'prop-types'
import React from 'react'
import { fetchUsersByIds } from '../utils/apiUtils.js'

export const ROLES = Object.freeze({
  MAIN_AUTHOR: 'Main author',
  COAUTHOR: 'Co-author',
  REVIEWER: 'Reviewer'
})

export const TYPES = Object.freeze({
  USER: 'user',
  GROUP: 'group'
})

export function MemberIcon({ type }) {
  return (type === TYPES.USER ? <PersonIcon /> : <GroupIcon />)
}
MemberIcon.propTypes = {
  type: PropTypes.string.isRequired
}

export function userToMember(user, role = ROLES.REVIEWER) {
  return {
    type: TYPES.USER,
    id: user.user_id,
    name: user.name,
    affiliation: user.affiliation,
    role: role
  }
}

export function groupToMember(group, usersDict = {}, role = ROLES.REVIEWER) {
  return {
    type: TYPES.GROUP,
    id: group.group_id,
    name: group.group_name,
    affiliation: usersDict[group.owner]?.name ?? group.owner,
    role: role
  }
}

export async function groupsToMembers(groups, api) {
  if (!groups.length) return []

  let ownersDict = {}
  if (api) {
    const ownerIds = groups.map(group => group.owner)
    const ownerUsers = await fetchUsersByIds(api, ownerIds)
    ownersDict = _.keyBy(ownerUsers, 'user_id')
  }
  return groups.map(group => groupToMember(group, ownersDict))
}

export const getMemberIdsByTypeAndRole = (members, type, role) => {
  return members
    .filter(member => member.type === type && member.role === role)
    .map(member => member.id)
}
