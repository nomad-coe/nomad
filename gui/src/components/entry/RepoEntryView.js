import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Divider, Card, CardContent, Grid, CardHeader, Typography, Link } from '@material-ui/core'
import { withApi } from '../api'
import { compose } from 'recompose'
import ApiDialogButton from '../ApiDialogButton'
import Quantity from '../Quantity'
import { withDomain } from '../domains'
import { Link as RouterLink } from 'react-router-dom'
import { DOI } from '../search/DatasetList'

class RepoEntryView extends React.Component {
  static styles = theme => ({
    root: {
      marginTop: theme.spacing.unit * 2
    },
    error: {
      marginTop: theme.spacing.unit * 2
    },
    cardContent: {
      paddingTop: 0
    },
    entryCards: {
      marginTop: theme.spacing.unit * 2
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

  static defaultState = {
    calcData: null,
    doesNotExist: false
  }

  state = {...RepoEntryView.defaultState}

  componentDidMount() {
    this.update()
  }

  componentDidUpdate(prevProps) {
    if (prevProps.api !== this.props.api ||
        prevProps.uploadId !== this.props.uploadId ||
        prevProps.calcId !== this.props.calcId) {
      this.setState({...RepoEntryView.defaultState})
      this.update()
    }
  }

  update() {
    const {uploadId, calcId} = this.props
    this.props.api.repo(uploadId, calcId).then(data => {
      this.setState({calcData: data, doesNotExist: false})
    }).catch(error => {
      this.setState({calcData: null, doesNotExist: false})
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

    const authors = loading ? null : calcData.authors

    if (this.state.doesNotExist) {
      return <Typography className={classes.error}>
          This entry does not exist.
      </Typography>
    }

    return (
      <div className={classes.root}>
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
                    <div style={{display: 'inline-grid'}}>
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
                          {ds.doi ? <span>&nbsp; (<DOI doi={ds.doi}/>)</span> : ''}
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
                  <Quantity quantity="calc_id" label={`${domain.entryLabel} id`} noWrap withClipboard {...quantityProps} />
                  <Quantity quantity="pid" label='PID' loading={loading} placeholder="not yet assigned" noWrap {...quantityProps} withClipboard />
                  <Quantity quantity="raw_id" label='raw id' loading={loading} noWrap {...quantityProps} withClipboard />
                  <Quantity quantity="external_id" label='external id' loading={loading} noWrap {...quantityProps} withClipboard />
                  <Quantity quantity="mainfile" loading={loading} noWrap ellipsisFront {...quantityProps} withClipboard />
                  <Quantity quantity="calc_hash" label={`${domain.entryLabel} hash`} loading={loading} noWrap {...quantityProps} />
                  <Quantity quantity="upload_id" label='upload id' {...quantityProps} noWrap withClipboard />
                  <Quantity quantity="upload_time" label='upload time' noWrap {...quantityProps} >
                    <Typography noWrap>
                      {new Date(calcData.upload_time * 1000).toLocaleString()}
                    </Typography>
                  </Quantity>
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

        <domain.EntryCards data={calcData} calcId={calcId} uploadId={uploadId} classes={{root: classes.entryCards}} />
      </div>
    )
  }
}

export default compose(withApi(false, true), withDomain, withStyles(RepoEntryView.styles))(RepoEntryView)
