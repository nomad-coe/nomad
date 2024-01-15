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
import React, {useCallback, useEffect, useState, useMemo, useContext, useRef} from 'react'
import PropTypes from 'prop-types'
import Markdown from '../Markdown'
import {
  Paper, IconButton, Tooltip, Typography,
  Box, Divider, makeStyles
} from '@material-ui/core'
import ClipboardIcon from '@material-ui/icons/Assignment'
import { HelpButton } from '../Help'
import { CopyToClipboard } from 'react-copy-to-clipboard'
import { guiBase, servicesUploadLimit } from '../../config'
import NewUploadButton from './NewUploadButton'
import ExampleUploadButton from './ExampleUploadButton'
import { useApi, withLoginRequired } from '../api'
import Page from '../Page'
import { useErrors } from '../errors'
import DetailsIcon from '@material-ui/icons/ArrowForward'
import { UploadButton } from '../nav/Routes'
import {
  addColumnDefaults, combinePagination, Datatable, DatatablePagePagination,
  DatatableTable, DatatableToolbar, DatatableToolbarActions
} from '../datatable/Datatable'
import TooltipButton from '../utils/TooltipButton'
import ReloadIcon from '@material-ui/icons/Replay'
import Quantity from '../Quantity'
import { SourceApiCall, SourceApiDialogButton } from '../buttons/SourceDialogButton'
import { useHistory } from 'react-router-dom'
import DeleteUploadsButton from './DeleteUploadsButton'
import UploadStatusIcon from './UploadStatusIcon'
import { formatTimestamp } from '../../utils'
import { SupportedCodes, UploadDocumentation } from './UploadOverview'
import { defaultFilterData } from '../search/FilterRegistry'

export const help = `
NOMAD allows you to upload data. After upload, NOMAD will process your data: it will
identify parseable files (*mainfiles*) and parse these.
The result will be a list of entries (one per identified mainfile).
Each entry is associated with metadata. This is data that NOMAD acquired from your files and that
describe your entry (e.g. chemical formula, used code, system type and symmetry, etc.).
Furthermore, you can provide your own metadata (comments, references, co-authors, etc.).
When an upload is created, it is created in the *staging area*, and it will only be visible
to you (the main author). You can add other users as *co-authors*. They will be listed
as co-authors of the upload and all its entries, and they can view and edit the upload.
You can also add *reviewers*, which are users who are allowed the view the upload. No
other users can see your data until you *publish* the upload.

#### Prepare and upload files

While the upload is in the staging area, you and the co-authors can add and remove files
to it. You can upload plain files, or \`*.zip\` or \`*.tar.gz\` archives. We encourage you
to add all input and output files, as well as any other auxiliary files that you might
have created, and to organize your files in a folder structure. Ideally, you should have a
separate directory for each entry, as NOMAD considers all files located in a directory where
there is a mainfile to be auxiliary files to this entry. If users want to download an entry,
they can download all files in the directory where the entry's mainfile is located. The
directory structure can be nested.

Drop your archive file(s) on the dropbox. You can also click the dropbox to select a file from
your hard drive. Alternatively, you can upload files via the given shell command.
Replace \`<local_file>\` with your archive file. After executing the command,
return here and press the reload button below).

### Limits

There is a limit of 10 unpublished uploads per user. Please accumulate all data into as
few uploads as possible. But, there is a also an upper limit of 32 GB per upload.
Please upload multiple archives, if you have more than 32 GB of data to upload.
Only uploads with at least one recognized entry can be published.

#### The uploads page

Here you will find all your unpublished and published uploads. You can see the
progress on the processing, you can review your uploads, and publish or delete them.

Click on an upload to see more details about its contents. Click on processed entries
to see their metadata, archive data, and a processing log. In the details view, you also
find buttons for editing user metadata, deleting uploads, and publishing uploads. Entries
cannot be published individually, only the upload as a whole, with all its entries, can be published.

#### Publishing and embargo

When you publish, you can set an *embargo* or publish your data as *Open Access* right away.
The embargo allows you to share data with selected users, create a DOI for your data, etc,
before the full upload data is made public. The embargo might last up to 36 months, and when
it expires the data becomes public automatically. While embargoed, some metadata (and datasets
created from this data) are publicly visible and findable, but only you and users you share
the upload with (i.e. users added as co-authors or reviewers) can view and download the raw
data and archive.

#### Processing errors

We distinguish between uploads that fail processing completely and uploads that contain
entries that could not be processed. The former might be caused by issues during the
upload, bad file formats, etc. The latter (far more common) case means that not all of the
identified mainfiles could be parsed. The processing logs of the failed entries might
provide some insight.

You cannot publish uploads that failed processing completely. In most
cases they wouldn't have any entries to publish anyway. In the case of failed processing of
some entries however, the data can still be published. You will be able to share it and create
DOIs for it, etc. The only shortcomings will be missing metadata (labeled *not processed*
or *unavailable*) and missing archive data. We continuously improve our parsers and
the now missing information might become available in the future automatically.

#### Co-Authors, References, Comments, Datasets, DOIs

You can edit additional *user metadata*. Some of these fields are set on the upload level,
others on the entry level. You can select and edit many entries at once. Edit buttons for
user metadata are available in many views on this web-page. For example, you can edit user
metadata when you click on an upload to open its details, and press the edit button there.
The metadata fields cannot be changed after the upload has been published (except for dataset members).
The documentation on the [user data page](${guiBase}/userdata) contains more information.
`
export const uploadsPageContext = React.createContext()
export function useUploadsPageContext() {
  return useContext(uploadsPageContext)
}

const useUploadCommandStyles = makeStyles(theme => ({
  root: {
    width: '100%'
  },
  commandContainer: {
    display: 'flex',
    flexDirection: 'row',
    alignItems: 'center'
  },
  commandMarkup: {
    flexGrow: 1,
    marginRight: theme.spacing(1),
    overflow: 'hidden'
  }
}))

function UploadCommands({uploadCommands}) {
  const classes = useUploadCommandStyles()

  return <div className={classes.root}>
    <div className={classes.commandContainer} role='upload-commands'>
      <div className={classes.commandMarkup}>
        <Markdown>{`
          \`\`\`
            ${uploadCommands.upload_command}
          \`\`\`
        `}</Markdown>
      </div>
      <CopyToClipboard text={uploadCommands.upload_command} onCopy={() => null}>
        <Tooltip title="Copy command to clipboard">
          <IconButton>
            <ClipboardIcon />
          </IconButton>
        </Tooltip>
      </CopyToClipboard>
      <Tooltip title="Alternative shell commands">
        <HelpButton maxWidth="md" heading="Alternative shell commands" content={`
          As an experienced shell and *curl* user, you can modify the commands to
          your liking.

          The given command can be modified. To see progress on large files, use
          \`\`\`
            ${uploadCommands.upload_progress_command}
          \`\`\`
          To \`tar\` and upload multiple folders in one command, use
          \`\`\`
          ${uploadCommands.upload_tar_command}
          \`\`\`

          ### Form data vs. streaming
          NOMAD accepts stream data (\`-T <local_file>\`) (like in the
          examples above) or multi-part form data (\`-X PUT -f file=@<local_file>\`):
          \`\`\`
          ${uploadCommands.upload_command_form}
          \`\`\`
          We generally recommend to use streaming, because form data can produce very
          large HTTP request on large files. Form data has the advantage of carrying
          more information (e.g. the file name) to our servers (see below).

          #### Upload names
          With multi-part form data (\`-X PUT -f file=@<local_file>\`), your upload will
          be named after the file by default. With stream data (\`-T <local_file>\`)
          there will be no default name. To set a custom name, you can use the URL
          parameter \`name\`:
          \`\`\`
          ${uploadCommands.upload_command_with_name}
          \`\`\`
          Make sure to user proper [URL encoding](https://www.w3schools.com/tags/ref_urlencode.asp)
          and shell encoding, if your name contains spaces or other special characters.
        `}/>
      </Tooltip>
    </div>
  </div>
}

const UploadActions = React.memo(function UploadActions({data}) {
  return <Tooltip title="Open this upload">
    <UploadButton component={IconButton} uploadId={data.upload_id}>
      <DetailsIcon />
    </UploadButton>
  </Tooltip>
})
UploadActions.propTypes = {
  data: PropTypes.object.isRequired
}

UploadCommands.propTypes = {
  uploadCommands: PropTypes.object.isRequired
}

export function UploadsPage() {
  const {api, user} = useApi()
  const errors = useErrors()
  const [apiData, setApiData] = useState(null)
  const data = useMemo(() => apiData?.response, [apiData])
  const [unpublished, setUnpublished] = useState(null)
  const [uploadCommands, setUploadCommands] = useState(null)
  const [selected, setSelected] = useState(new Set())
  const fetchTimer = useRef(-1)
  const [pagination, setPagination] = useState({
    page_size: 10,
    page: 1,
    order_by: 'upload_create_time',
    order: 'desc'
  })
  const history = useHistory()

  const fetchData = useCallback(() => {
    const {page_size, page, order_by, order} = pagination
    const url = `/uploads?page_size=${page_size}&page=${page}&order_by=${(order_by === 'published' ? 'publish_time' : order_by)}&order=${order}`
    api.get(url, null, {returnRequest: true})
      .then(uploads => {
        setApiData(uploads)
        const deleteStatus = uploads.response.data.filter(upload => upload?.current_process === 'delete_upload' && (upload?.process_status === 'PENDING' || upload?.process_status === 'RUNNING'))
        if (deleteStatus.length === 0) {
          if (fetchTimer.current !== -1) {
            clearInterval(fetchTimer.current)
          }
          fetchTimer.current = -1
        } else {
          if (fetchTimer.current === -1) {
            fetchTimer.current = setInterval(fetchData, 2000)
          }
        }
      })
      .catch(errors.raiseError)
    api.get(`/uploads?is_published=false&roles=main_author&page_size=10000&order_by=${(order_by === 'published' ? 'publish_time' : order_by)}&order=${order}`)
      .then(response => {
        setUnpublished(response)
      })
      .catch(errors.raiseError)
  }, [pagination, setApiData, setUnpublished, errors, api])

  const handleReload = () => {
    fetchData()
  }

  // This history dependency makes sure that the uploads are reloaded after page navigating
  // back to the UploadsPage. E.g. if an upload was deleted on the UploadPage, we are
  // history.pushed back to the UploadsPage and want to reload the list of uploads.
  useEffect(() => {
    fetchData()
    return () => {
      if (fetchTimer.current !== -1) clearInterval(fetchTimer.current)
      fetchTimer.current = -1
    }
  }, [fetchData, history])

  const isDisabled = unpublished
    ? (unpublished.pagination ? unpublished.pagination.total >= servicesUploadLimit : false)
    : false

  useEffect(() => {
    api.get('/uploads/command-examples')
      .then(setUploadCommands)
      .catch(errors.raiseError)
  }, [api, errors, setUploadCommands])

  const selectedUploads = useMemo(() => selected === 'all'
    ? unpublished.data
    : data?.data?.filter(upload => selected.has(upload.upload_id)
  ), [data, selected, unpublished])

  // A map to store whether each of the upload object's main author is the current user
  const currentUserUploadsMap = useMemo(() => {
    const map = new Map()
    if (data && user) {
      data.data.forEach((d) => {
        map[d.upload_id] = d.main_author === user.sub
      })
    }
    return map
  }, [data, user])

  // boolean to disable the Delete button if any of the upload's main author is not the current user
  const disableDeleteButton = useMemo(() => {
    if (selected === 'all' && data) {
      return data.data.some(({upload_id}) => !currentUserUploadsMap[upload_id])
    } else {
      return Array.from(selected).some((id) => !currentUserUploadsMap[id])
    }
  }, [data, selected, currentUserUploadsMap])

  const deleteUploads = useCallback(() => {
    const promises = selectedUploads.map(upload => api.delete(`uploads/${upload.upload_id}`))
    return Promise.allSettled(promises)
  }, [selectedUploads, api])

  const handleDeleteSubmitted = useCallback(() => {
    deleteUploads().then(responses => {
      responses.forEach((response) => {
        if (response.status === 'fulfilled') {
          // No error occurred
        } else if (response.status === 'rejected') {
          errors.raiseError(new Error(response.reason.apiMessage))
        } else {
          errors.raiseError(new Error('Unexpected error!'))
        }
      })
      setSelected(new Set([]))
      fetchData()
    })
  }, [deleteUploads, errors, fetchData])

  // Setup columns
  const columns = useMemo(() => {
    const columns = [
      {
        key: 'upload_create_time',
        render: upload => formatTimestamp(upload.upload_create_time)
      },
      {key: 'upload_name'},
      {
        key: 'upload_id',
        render: upload => <Quantity quantity={'upload_id'} noLabel noWrap withClipboard data={upload}/>
      },
      {key: 'last_status_message', label: 'Status'},
      {key: 'entries', render: upload => upload.entries, align: 'center'},
      {key: 'published', render: upload => <UploadStatusIcon data={upload} user={user} />, align: 'center'}
    ]
    addColumnDefaults(columns, {align: 'left'}, defaultFilterData)
    return columns
  }, [user])

  return (
    <Page loading={!(data && uploadCommands)}>
      <Box marginBottom={2}>
        <Typography>
          You can either create an upload and upload files through the browser:
        </Typography>
      </Box>
      <Box alignItems='center' style={{display: 'flex'}}>
        <NewUploadButton color="primary" disabled={isDisabled}/>
        <Typography color='initial' style={{padding: '10px'}} >or</Typography>
        <ExampleUploadButton color="inherit" disabled={isDisabled} />
        <Box display="inline-block" marginLeft={2}>
          {isDisabled && <Typography color="error" role='error-maximum-number-of-unpublished'>
            You have reached maximum number of unpublished uploads!
          </Typography>}
        </Box>
      </Box>
      <Box marginTop={4}>
        <Typography>
          Or, you can create an upload by sending a file-archive via shell command:
        </Typography>
      </Box>
      <Box>
        {uploadCommands && <UploadCommands uploadCommands={uploadCommands}/>}
      </Box>
      <Box>
        <Typography>
          Visit our documentation on <UploadDocumentation>uploading files</UploadDocumentation> or show <SupportedCodes>supported codes</SupportedCodes>.
        </Typography>
      </Box>
      {(data?.pagination?.total || 0) > 0 && <React.Fragment>
        <Box marginTop={2} marginBottom={2}>
          <Divider/>
        </Box>
        <Paper>
          <Datatable
            columns={columns} shownColumns={columns.map(column => column.key)}
            sortingColumns={['upload_create_time', 'upload_name', 'last_status_message', 'published']}
            data={data.data || []}
            pagination={combinePagination(pagination, data.pagination)}
            onPaginationChanged={setPagination}
            selected={selected}
            getId={option => option.upload_id}
            onSelectedChanged={setSelected}
          >
            <DatatableToolbar title="Your existing uploads">
              <DatatableToolbarActions>
                <TooltipButton
                  title="Reload the uploads"
                  component={IconButton}
                  onClick={handleReload}
                >
                  <ReloadIcon/>
                </TooltipButton>
                <SourceApiDialogButton maxWidth="lg" fullWidth>
                  <SourceApiCall {...apiData} />
                </SourceApiDialogButton>
              </DatatableToolbarActions>
              <DatatableToolbarActions selection>
                <DeleteUploadsButton
                  isIcon
                  disabled={disableDeleteButton}
                  uploads={selectedUploads}
                  onConfirm={handleDeleteSubmitted}
                />
              </DatatableToolbarActions>
            </DatatableToolbar>
            <DatatableTable actions={UploadActions}>
              <DatatablePagePagination />
            </DatatableTable>
          </Datatable>
        </Paper>
      </React.Fragment>}
    </Page>
  )
}

export default withLoginRequired(UploadsPage)
