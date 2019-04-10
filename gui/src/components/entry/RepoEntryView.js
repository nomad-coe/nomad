import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Divider, Card, CardContent, Grid, CardHeader, Fab } from '@material-ui/core'
import { withApi } from '../api'
import { compose } from 'recompose'
import RawFiles from './RawFiles'
import Download from './Download'
import DownloadIcon from '@material-ui/icons/CloudDownload'
import ApiDialogButton from '../ApiDialogButton'
import Quantity from '../Quantity'
import { withDomain } from '../domains'

class RepoEntryView extends React.Component {
  static styles = theme => ({
    root: {},
    title: {
      marginBottom: theme.spacing.unit * 3
    },
    content: {
      marginTop: theme.spacing.unit * 3
    },
    quantityContainer: {
      display: 'flex'
    },
    quantityColumn: {
      display: 'flex',
      flexDirection: 'column'
    },
    quantityRow: {
      display: 'flex',
      flexDirection: 'row',
      flexWrap: 'wrap',
      marginBottom: theme.spacing.unit
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
    calcData: null
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
      this.props.raiseError(error)
    })
  }

  render() {
    const { classes, domain, ...calcProps } = this.props
    const calcData = this.state.calcData || calcProps
    const loading = !this.state.calcData
    const { uploadId, calcId } = calcProps

    const mainfile = calcData.mainfile
    const calcPath = mainfile ? mainfile.substring(0, mainfile.lastIndexOf('/')) : null

    const authors = loading ? null : calcData.authors

    return (
      <div className={classes.root}>
        <div className={classes.content}>

          <div className={classes.title}>
            <Quantity label="chemical formula" typography="h3" loading={loading}>
              {calcData.formula}
            </Quantity>
          </div>

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
                <CardContent classes={{root: classes.cardContent}}>
                  <div className={classes.quantityColumn}>
                    <div className={classes.quantityColumn}>
                      <Quantity label='comment' loading={loading} placeholder='no comment'>
                        {calcData.comment}
                      </Quantity>
                      <Quantity label='references' loading={loading} placeholder='no references'>
                        {calcData.references ? calcData.references.map(ref => (<a key={ref.id} href={ref.value}>{ref.value}</a>)) : null}
                      </Quantity>
                      <Quantity label='authors' loading={loading}>
                        {authors ? authors.map(author => author.name) : null}
                      </Quantity>
                      <Quantity label='datasets' loading={loading} placeholder='no datasets'>
                        {calcData.datasets ? calcData.datasets.map(ds => `${ds.name}${ds.doi ? ` (${ds.doi})` : ''}`).join(', ') : null}
                      </Quantity>
                    </div>
                  </div>
                </CardContent>
              </Card>
            </Grid>

            <Grid item xs={5}>
              <Card>
                <CardHeader title="Ids / processing" />
                <CardContent classes={{root: classes.cardContent}}>
                  <div className={classes.quantityColumn} style={{maxWidth: 350}}>
                    <Quantity label='PID' loading={loading} noWrap>
                      {calcData.pid ? <b>{calcData.pid}</b> : <i>not yet assigned</i>}
                    </Quantity>
                    <Quantity label='upload id' noWrap>
                      {calcData.upload_id}
                    </Quantity>
                    <Quantity label='upload time' noWrap>
                      {new Date(calcData.upload_time * 1000).toLocaleString()}
                    </Quantity>
                    <Quantity label='calculation id' noWrap>
                      {calcData.calc_id}
                    </Quantity>
                    <Quantity label='mainfile' loading={loading} noWrap>
                      {mainfile}
                    </Quantity>
                    <Quantity label='calculation hash' loading={loading} noWrap>
                      {calcData.calc_hash}
                    </Quantity>
                    <Quantity label='last processing' loading={loading} noWrap>
                      {calcData.last_processing ? new Date(calcData.last_processing * 1000).toLocaleString() : <i>not processed</i>}
                    </Quantity>
                    <Quantity label='processing version' loading={loading} noWrap>
                      {calcData.last_processing ? `${calcData.nomad_version}/${calcData.nomad_commit}` : <i>not processed</i>}
                    </Quantity>
                  </div>
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
