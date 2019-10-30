import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Divider, Card, CardContent, Grid, CardHeader, Fab, Typography, Link } from '@material-ui/core'
import { withApi } from '../api'
import { compose } from 'recompose'
import Download from './Download'
import DownloadIcon from '@material-ui/icons/CloudDownload'
import ApiDialogButton from '../ApiDialogButton'
import Quantity from '../Quantity'
import { withDomain } from '../domains'
import { Link as RouterLink } from 'react-router-dom'

class RepoEntryView extends React.Component {
  static styles = theme => ({
    root: {},
    error: {
      marginTop: theme.spacing.unit * 2
    },
    title: {
      marginBottom: theme.spacing.unit * 3
    },
    content: {
      marginTop: theme.spacing.unit * 3
    },
    downloadFab: {
      zIndex: 1,
      right: 32,
      bottom: 32,
      position: 'fixed !important'
    },
    cardContent: {
      paddingTop: 0
    },
    entryCards: {
      marginTop: theme.spacing.unit * 3
    }
  })

  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired,
    uploadId: PropTypes.string.isRequired,
    calcId: PropTypes.string.isRequired,
    domain: PropTypes.object.isRequired
  }

  state = {
    calcData: null,
    doesNotExist: false
  }

  componentDidMount() {
    this.update()
  }

  componentDidUpdate(prevProps) {
    if (prevProps.api !== this.props.api) {
      this.update()
    }
  }

  update() {
    const {uploadId, calcId} = this.props
    this.props.api.repo(uploadId, calcId).then(data => {
      this.setState({calcData: data})
    }).catch(error => {
      this.setState({calcData: null})
      if (error.name === 'DoesNotExist') {
        this.setState({doesNotExist: true})
      } else {
        this.props.raiseError(error)
      }
    })
  }

  render() {
    const { classes, domain, ...calcProps } = this.props
    const calcData = this.state.calcData || calcProps
    const loading = !this.state.calcData
    const { uploadId, calcId } = calcProps
    const quantityProps = {data: calcData, loading: loading}

    const mainfile = calcData.mainfile
    const calcPath = mainfile ? mainfile.substring(0, mainfile.lastIndexOf('/')) : null

    const authors = loading ? null : calcData.authors

    if (this.state.doesNotExist) {
      return <Typography className={classes.error}>
          This entry does not exist.
      </Typography>
    }

    return (
      <div className={classes.root}>
        <div className={classes.content}>
          <Grid container spacing={24}>
            <Grid item xs={7}>
              <Card>
                <CardHeader
                  title="Metadata"
                  action={<ApiDialogButton title="Repository JSON" data={calcData} />}
                />
                <CardContent classes={{root: classes.cardContent}}>
                  <domain.EntryOverview data={calcData} loading={loading} />
                </CardContent>
                <Divider />
                <CardContent>
                  <Quantity column>
                    <Quantity quantity='comment' placeholder='no comment' {...quantityProps} />
                    <Quantity quantity='references' placeholder='no references' {...quantityProps}>
                      <div>
                        {(calcData.references || []).map(ref => <Typography key={ref} noWrap>
                          <a href={ref}>{ref}</a>
                        </Typography>)}
                      </div>
                    </Quantity>
                    <Quantity quantity='authors' {...quantityProps}>
                      <Typography>
                        {(authors || []).map(author => author.name).join('; ')}
                      </Typography>
                    </Quantity>
                    <Quantity quantity='datasets' placeholder='no datasets' {...quantityProps}>
                      <div>
                        {(calcData.datasets || []).map(ds => (
                          <Typography key={ds.id}>
                            <Link component={RouterLink} to={`/dataset/id/${ds.id}`}>{ds.name}</Link>
                            {ds.doi ? <span>&nbsp; (<Link href={ds.doi}>{ds.doi}</Link>)</span> : ''}
                          </Typography>))}
                      </div>
                    </Quantity>
                  </Quantity>
                </CardContent>
              </Card>
            </Grid>

            <Grid item xs={5}>
              <Card>
                <CardHeader title="Ids / processing" />
                <CardContent classes={{root: classes.cardContent}}>
                  <Quantity column style={{maxWidth: 350}}>
                    <Quantity quantity="pid" label='PID' loading={loading} placeholder="not yet assigned" noWrap {...quantityProps} withClipboard />
                    <Quantity quantity="upload_id" label='upload id' {...quantityProps} noWrap withClipboard />
                    <Quantity quantity="upload_time" label='upload time' noWrap {...quantityProps} >
                      <Typography noWrap>
                        {new Date(calcData.upload_time * 1000).toLocaleString()}
                      </Typography>
                    </Quantity>
                    <Quantity quantity="calc_id" label={`${domain.entryLabel} id`} noWrap withClipboard {...quantityProps} />
                    <Quantity quantity='mainfile' loading={loading} noWrap {...quantityProps} withClipboard />
                    <Quantity quantity="calc_hash" label={`${domain.entryLabel} hash`} loading={loading} noWrap {...quantityProps} />
                    <Quantity quantity="raw_id" label='raw id' loading={loading} noWrap {...quantityProps} withClipboard />
                    <Quantity quantity="external_id" label='external id' loading={loading} noWrap {...quantityProps} withClipboard />
                    <Quantity quantity="last_processing" label='last processing' loading={loading} placeholder="not processed" noWrap {...quantityProps}>
                      <Typography noWrap>
                        {new Date(calcData.last_processing * 1000).toLocaleString()}
                      </Typography>
                    </Quantity>
                    <Quantity quantity="last_processing" label='processing version' loading={loading} noWrap placeholder="not processed" {...quantityProps}>
                      <Typography noWrap>
                        {calcData.nomad_version}/{calcData.nomad_commit}
                      </Typography>
                    </Quantity>
                  </Quantity>
                </CardContent>
              </Card>
            </Grid>
          </Grid>

          <domain.EntryCards data={calcData} classes={{root: classes.entryCards}} />

          <Download
            disabled={!mainfile} tooltip="download all raw files for calculation"
            classes={{root: classes.downloadFab}}
            component={Fab} className={classes.downloadFab} color="primary" size="medium"
            url={`raw/${uploadId}/${calcPath}/*`} fileName={`${calcId}.zip`}
          >
            <DownloadIcon />
          </Download>

        </div>
      </div>
    )
  }
}

export default compose(withApi(false, true), withDomain, withStyles(RepoEntryView.styles))(RepoEntryView)
