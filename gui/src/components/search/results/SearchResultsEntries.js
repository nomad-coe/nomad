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
import React, { useState, useCallback, useMemo } from 'react'
import PropTypes from 'prop-types'
import { makeStyles } from '@material-ui/core/styles'
import { Link, Typography, Tooltip, IconButton, Button } from '@material-ui/core'
import { useHistory, Link as RouterLink } from 'react-router-dom'
import NavigateNextIcon from '@material-ui/icons/NavigateNext'
import PublicIcon from '@material-ui/icons/Public'
import UploaderIcon from '@material-ui/icons/AccountCircle'
import EditUserMetadataDialog from '../../EditUserMetadataDialog'
import SharedIcon from '@material-ui/icons/SupervisedUserCircle'
import PrivateIcon from '@material-ui/icons/VisibilityOff'
import DownloadButton from '../../DownloadButton'
import { domainData } from '../../domainData'
import { domainComponents } from '../../domainComponents'
import { authorList, nameList } from '../../../utils'
import NewDataTable from '../../NewDataTable'
import { useApi } from '../../apiV1'
import Quantity from '../../Quantity'
import searchQuantities from '../../../searchQuantities'

/**
 * Displays the list of search results for entries.
 */
export function Published(props) {
  const {user} = useApi()
  const {entry} = props
  if (entry.published) {
    if (entry.with_embargo) {
      if (user && entry.uploader.user_id === user.sub) {
        if (entry.owners.length === 1) {
          return <Tooltip title="published with embargo by you and only accessible by you">
            <UploaderIcon color="error" />
          </Tooltip>
        } else {
          return <Tooltip title="published with embargo by you and only accessible to you and users you shared the data with">
            <SharedIcon color="error" />
          </Tooltip>
        }
      } else if (user && entry.owners.find(user => user.user_id === user.sub)) {
        return <Tooltip title="published with embargo and shared with you">
          <SharedIcon color="error" />
        </Tooltip>
      } else {
        if (user) {
          return <Tooltip title="published with embargo and not accessible by you">
            <PrivateIcon color="error" />
          </Tooltip>
        } else {
          return <Tooltip title="published with embargo and might become accessible after login">
            <PrivateIcon color="error" />
          </Tooltip>
        }
      }
    } else {
      return <Tooltip title="published and accessible by everyone">
        <PublicIcon color="primary" />
      </Tooltip>
    }
  } else {
    return <Tooltip title="you have not published this entry yet">
      <UploaderIcon color="error"/>
    </Tooltip>
  }
}

const columns = {
  formula: {
    label: 'Formula',
    render: row => row.results.material.chemical_formula_hill,
    supportsSort: true,
    description: searchQuantities['results.material.chemical_formula_hill'].description
  },
  method: {
    label: 'Method name',
    render: row => row.results.method.method_name,
    supportsSort: true,
    description: searchQuantities['results.method.method_name'].description
  },
  program_name: {
    label: 'Program name',
    render: row => row.results.method.simulation?.program_name,
    supportsSort: true,
    description: searchQuantities['results.method.simulation.program_name'].description
  },
  basis_set: {
    label: 'Basis set',
    render: row => row.results.method.simulation?.dft?.basis_set_name,
    supportsSort: true,
    description: searchQuantities['results.method.simulation.dft.basis_set_name'].description
  },
  functional_family: {
    label: 'XC functional type',
    render: row => row.results.method.simulation?.dft?.xc_functional_type,
    supportsSort: true,
    description: searchQuantities['results.method.simulation.dft.xc_functional_type'].description
  },
  structural_type: {
    label: 'Structural type',
    render: row => row.results.material.structural_type,
    supportsSort: true,
    description: searchQuantities['results.material.structural_type'].description
  },
  crystal_system: {
    label: 'Crystal system',
    render: row => row.results.material.symmetry?.crystal_system,
    supportsSort: true,
    ellipsisFront: true,
    description: searchQuantities['results.material.symmetry.crystal_system'].description
  },
  space_group_symbol: {
    label: 'Space group symbol',
    render: row => row.results.material.symmetry?.space_group_symbol,
    supportsSort: true,
    ellipsisFront: true,
    description: searchQuantities['results.material.symmetry.space_group_symbol'].description
  },
  space_group_number: {
    label: 'Space group number',
    render: row => row.results.material.symmetry?.space_group_number,
    supportsSort: true,
    ellipsisFront: true,
    description: searchQuantities['results.material.symmetry.space_group_number'].description
  },
  mainfile: {
    label: 'Mainfile',
    render: row => row.mainfile,
    supportsSort: true,
    ellipsisFront: true,
    description: searchQuantities['mainfile'].description
  },
  upload_time: {
    label: 'Upload time',
    render: row => new Date(row.upload_time).toLocaleString(),
    supportsSort: true,
    description: searchQuantities['upload_time'].description
  },
  authors: {
    label: 'Authors',
    render: row => authorList(row),
    supportsSort: true,
    description: 'The authors of this entry. This includes the uploader and its co-authors.'
  },
  co_authors: {
    label: 'co-Authors',
    render: row => nameList(row.authors),
    supportsSort: false,
    description: 'The people that this entry was co authored with'
  },
  shared_with: {
    label: 'Shared with',
    render: row => nameList(row.authors),
    supportsSort: false,
    description: 'The people that this entry was shared with'
  },
  uploader: {
    label: 'Uploader',
    render: row => row.uploader.name,
    supportsSort: true,
    description: 'The uploader of this entry.'
  },
  comment: {
    label: 'Comment',
    render: row => row.comment,
    supportsSort: false,
    description: 'User provided comment on this entry'
  },
  references: {
    label: 'References',
    render: row => {
      const refs = row.references || []
      if (refs.length > 0) {
        return (
          <div style={{display: 'inline'}}>
            {refs.map((ref, i) => <span key={ref}>
              <Link href={ref}>{ref}</Link>{(i + 1) < refs.length ? ', ' : <React.Fragment/>}
            </span>)}
          </div>
        )
      } else {
        return <i>no references</i>
      }
    },
    supportsSort: true
  },
  datasets: {
    label: 'Datasets',
    render: entry => {
      const datasets = entry.datasets || []
      if (datasets.length > 0) {
        return datasets.map(dataset => dataset.name).join(', ')
      } else {
        return <i>no datasets</i>
      }
    },
    supportsSort: false,
    description: 'The dataset names that this entry belongs to.'
  }
}

const selectedColumns = ['formula', 'method', 'program_name', 'structural_type', 'crystal_system']

const useStyles = makeStyles(theme => ({
  root: {
    height: '100%'
  },
  entryDetails: {
    paddingTop: theme.spacing(2),
    paddingLeft: theme.spacing(2),
    paddingRight: theme.spacing(2)
  },
  entryDetailsContents: {
    display: 'flex',
    maxWidth: 1024,
    margin: 'auto'
  },
  entryDetailsRow: {
    paddingRight: theme.spacing(3)
  },
  entryDetailsActions: {
    display: 'flex',
    flexBasis: 'auto',
    flexGrow: 0,
    flexShrink: 0,
    justifyContent: 'flex-end',
    marginBottom: theme.spacing(1),
    marginTop: theme.spacing(2)
  }
}))
const SearchResultsEntries = React.memo(({
  data,
  query,
  editable,
  className,
  ...rest
}) => {
  const {user} = useApi()
  const [selected, setSelected] = useState([])
  const styles = useStyles()
  const total = data?.pagination && data.pagination.total
  const history = useHistory()
  const entryPagePathPrefix = undefined

  // The access column is hidden if user is not logged in
  const visibleColumns = useMemo(() => {
    if (user) {
      return {...columns,
        published: {
          label: 'Access',
          align: 'center',
          render: (entry) => <Published entry={entry} />
        }
      }
    }
    return columns
  }, [user])

  // Decides whether actions should be shown for an entry
  const showEntryActions = useCallback((row) => {
    if (row.with_embargo && !(user && row.owners.find(owner => owner.user_id === user.sub))) {
      return false
    }
    return true
  }, [user])

  const handleViewEntryPage = useCallback((event, row) => {
    event.stopPropagation()
    const prefix = entryPagePathPrefix || ''
    const url = `${prefix}/entry/id/${row.upload_id}/${row.calc_id}`
    history.push(url)
  }, [entryPagePathPrefix, history])

  const renderEntryActions = useCallback((row, selected) => {
    if (showEntryActions(row)) {
      return <Tooltip title="Show raw files and archive">
        <IconButton style={selected ? {color: 'white'} : null} onClick={event => handleViewEntryPage(event, row)}>
          <NavigateNextIcon />
        </IconButton>
      </Tooltip>
    } else {
      return ''
    }
  }, [handleViewEntryPage, showEntryActions])

  const renderEntryDetails = useCallback((row) => {
    const domain = (row.domain && domainData[row.domain]) || domainData.dft
    const domainComponent = (row.domain && domainComponents[row.domain]) || domainComponents.dft

    return (<div className={styles.entryDetails}>
      <div className={styles.entryDetailsContents}>
        <div className={styles.entryDetailsRow}>
          <domainComponent.EntryDetails data={row} />
        </div>

        <div className={styles.entryDetailsRow} style={{flexGrow: 1}}>
          <Quantity className={styles.entryDetailsRow} column>
            <Quantity quantity='comment' placeholder='no comment' data={row} />
            <Quantity quantity='references' placeholder='no references' data={row}>
              {row.references && <div style={{display: 'inline-grid'}}>
                {(row.references || []).map(ref => <Typography key={ref} noWrap>
                  <Link href={ref}>{ref}</Link>
                </Typography>)}
              </div>}
            </Quantity>
            <Quantity quantity='authors' data={row}>
              <Typography>
                {authorList(row)}
              </Typography>
            </Quantity>
            <Quantity quantity='datasets' placeholder='no datasets' data={row}>
              <div>
                {(row.datasets || []).map(ds => (
                  <Typography key={ds.dataset_id}>
                    <Link component={RouterLink} to={`/dataset/id/${ds.dataset_id}`}>{ds.name}</Link>
                    {ds.doi ? <span>&nbsp; (<Link href={`https://dx.doi.org/${ds.doi}`}>{ds.doi}</Link>)</span> : <React.Fragment/>}
                  </Typography>))}
              </div>
            </Quantity>
          </Quantity>
        </div>

        <div className={styles.entryDetailsRow} style={{maxWidth: '33%', paddingRight: 0}}>
          <Quantity column >
            {/* <Quantity quantity="pid" label='PID' placeholder="not yet assigned" noWrap data={row} withClipboard /> */}
            <Quantity quantity="calc_id" label={`${domain ? domain.entryLabel : 'entry'} id`} noWrap withClipboard data={row} />
            <Quantity quantity="raw_id" label={`raw id`} noWrap withClipboard data={row} />
            <Quantity quantity="external_id" label={`external id`} noWrap withClipboard data={row} />
            <Quantity quantity='mainfile' noWrap ellipsisFront data={row} withClipboard />
            <Quantity quantity="upload_id" label='upload id' data={row} noWrap withClipboard />
          </Quantity>
        </div>
      </div>

      <div className={styles.entryDetailsActions}>
        {showEntryActions(row) &&
          <Button color="primary" onClick={event => handleViewEntryPage(event, row)}>
            Show raw files and archive
          </Button>
        }
      </div>
    </div>)
  }, [handleViewEntryPage, showEntryActions, styles])

  const totalNumber = total || 0
  const example = selected && selected.length > 0 ? data?.data.find(d => d.calc_id === selected[0]) : data?.data[0]
  const selectQuery = (selected && selected.length > 0) ? {calc_id: selected, owner: query['owner']} : query
  const createActions = useCallback((props, moreActions) => <>
    {example && editable ? <EditUserMetadataDialog
      example={example} total={selected === null ? totalNumber : selected.length}
      onEditComplete={() => this.props.onEdit()}
      {...props}
    /> : ''}
    <DownloadButton
      tooltip="Download files"
      {...props}/>
    {moreActions}
  </>, [editable, example, selected, totalNumber])
  const selectActions = createActions({query: selectQuery, buttonProps: {color: 'secondary'}})

  return <NewDataTable
    entityLabels={['entry', 'entries']}
    selectActions={selectActions}
    id={row => row.entry_id}
    total={total}
    columns={visibleColumns}
    selectedColumns={selectedColumns}
    selectedColumnsKey="entries"
    entryDetails={renderEntryDetails}
    entryActions={renderEntryActions}
    data={data?.data || []}
    rows={data?.data.length || 0}
    selected={selected}
    onSelectionChanged={setSelected}
    {...rest}
  />
})
SearchResultsEntries.propTypes = {
  data: PropTypes.object,
  query: PropTypes.object,
  editable: PropTypes.bool,
  className: PropTypes.string
}

export default SearchResultsEntries
