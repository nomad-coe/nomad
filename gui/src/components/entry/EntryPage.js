import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Tab, Tabs } from '@material-ui/core'
import ArchiveEntryView from './ArchiveEntryView'
import ArchiveLogView from './ArchiveLogView'
import RepoEntryView from './RepoEntryView'
import { withApi } from '../api'
import { compose } from 'recompose'
import KeepState from '../KeepState'
import { guiBase } from '../../config'
import { matchPath } from 'react-router-dom'

export const help = `
The *raw files* tab, will show you all files that belong to the entry and offers a download
on individual, or all files. The files can be selected and downloaded. You can also
view the contents of some files directly here on this page.

The *archive* tab, shows you the parsed data as a tree
data structure. This view is connected to NOMAD's [meta-info](${guiBase}/metainfo), which acts a schema for
all parsed data.

The *log* tab, will show you a log of the entry's processing.
`

class EntryPage extends React.Component {
  static styles = theme => ({
    root: {
    },
    content: {
      padding: `0 ${theme.spacing(3)}px`,
      maxWidth: 1024,
      margin: 'auto'
    }
  })

  static propTypes = {
    classes: PropTypes.object.isRequired,
    location: PropTypes.object,
    match: PropTypes.object.isRequired,
    history: PropTypes.object.isRequired
  }

  render() {
    const { classes, location, match, history } = this.props
    const { path } = match

    const entryMatch = matchPath(location.pathname, {
      path: `${path}/:uploadId?/:calcId?/:tab?`
    })
    const { tab, calcId, uploadId } = entryMatch.params

    if (calcId && uploadId) {
      const calcProps = { calcId: calcId, uploadId: uploadId }
      return (
        <div className={classes.root}>
          <Tabs
            className={classes.tabs}
            value={tab || 'raw'}
            onChange={(_, value) => history.push(`${match.url}/${uploadId}/${calcId}/${value}`)}
            indicatorColor="primary"
            textColor="primary"
            variant="fullWidth"
          >
            <Tab label="Raw data" value="raw" />
            <Tab label="Archive" value="archive"/>
            <Tab label="Logs" value="logs"/>
          </Tabs>

          <div className={classes.content}>
            <KeepState visible={tab === 'raw' || tab === undefined} render={props => <RepoEntryView {...props} />} {...calcProps} />
            <KeepState visible={tab === 'archive'} render={props => <ArchiveEntryView {...props} />} {...calcProps} />
            <KeepState visible={tab === 'logs'} render={props => <ArchiveLogView {...props} />} {...calcProps} />
          </div>
        </div>
      )
    } else {
      return ''
    }
  }
}

export default compose(withApi(false, true), withStyles(EntryPage.styles))(EntryPage)
