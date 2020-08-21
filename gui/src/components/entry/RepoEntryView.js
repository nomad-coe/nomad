import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Divider, Card, CardContent, Grid, CardHeader, Typography, Link } from '@material-ui/core'
import { withApi } from '../api'
import { compose } from 'recompose'
import ApiDialogButton from '../ApiDialogButton'
// import Structure from '../visualization/Structure'
import Quantity from '../Quantity'
import { Link as RouterLink } from 'react-router-dom'
import { DOI } from '../search/DatasetList'
import { domains } from '../domains'
import { EntryPageContent } from './EntryPage'

class RepoEntryView extends React.Component {
  static styles = theme => ({
    root: {
      marginTop: theme.spacing(2)
    },
    error: {
      marginTop: theme.spacing(2)
    },
    cardContent: {
      paddingTop: 0
    },
    entryCards: {
      marginTop: theme.spacing(2)
    },
    structureViewer: {
      height: '25rem',
      padding: '0'
    }
  })

  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired,
    uploadId: PropTypes.string.isRequired,
    calcId: PropTypes.string.isRequired
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
    const { classes, ...calcProps } = this.props
    const calcData = this.state.calcData || calcProps
    const loading = !this.state.calcData
    const { uploadId, calcId } = calcProps
    const quantityProps = {data: calcData, loading: loading}

    const authors = loading ? null : calcData.authors
    const domain = calcData.domain && domains[calcData.domain]

    let entryHeader = 'Entry metadata'
    if (domain) {
      entryHeader = domain.entryTitle(calcData)
    }

    if (this.state.doesNotExist) {
      return <Typography className={classes.error}>
          This entry does not exist.
      </Typography>
    }

    return (
      <EntryPageContent className={classes.root} fixed>
        <Grid container spacing={2}>
          <Grid item xs={7}>
            <Card>
              <CardHeader
                title={entryHeader}
                action={<ApiDialogButton title="Repository JSON" data={calcData} />}
              />
              <CardContent classes={{root: classes.cardContent}}>
                {domain && <domain.EntryOverview data={calcData} loading={loading} />}
              </CardContent>
              <Divider />
              <CardContent>
                <Quantity column>
                  <Quantity quantity='comment' placeholder='no comment' {...quantityProps} />
                  <Quantity quantity='references' placeholder='no references' {...quantityProps}>
                    {calcData.references &&
                      <div style={{display: 'inline-grid'}}>
                        {calcData.references.map(ref => <Typography key={ref} noWrap>
                          <a href={ref}>{ref}</a>
                        </Typography>)}
                      </div>}
                  </Quantity>
                  <Quantity quantity='authors' {...quantityProps}>
                    <Typography>
                      {(authors || []).map(author => author.name).join('; ')}
                    </Typography>
                  </Quantity>
                  <Quantity quantity='datasets' placeholder='no datasets' {...quantityProps}>
                    {calcData.datasets &&
                      <div>
                        {calcData.datasets.map(ds => (
                          <Typography key={ds.dataset_id}>
                            <Link component={RouterLink} to={`/dataset/id/${ds.dataset_id}`}>{ds.name}</Link>
                            {ds.doi ? <span>&nbsp; (<DOI doi={ds.doi}/>)</span> : ''}
                          </Typography>))}
                      </div>}
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
                  <Quantity quantity="calc_id" label={`${domain ? domain.entryLabel : 'entry'} id`} noWrap withClipboard {...quantityProps} />
                  <Quantity quantity={entry => entry.encyclopedia.material.material_id} label='material id' loading={loading} noWrap {...quantityProps} withClipboard />
                  <Quantity quantity="raw_id" label='raw id' loading={loading} noWrap {...quantityProps} withClipboard />
                  <Quantity quantity="external_id" label='external id' loading={loading} noWrap {...quantityProps} withClipboard />
                  <Quantity quantity="mainfile" loading={loading} noWrap ellipsisFront {...quantityProps} withClipboard />
                  <Quantity quantity="upload_id" label='upload id' {...quantityProps} noWrap withClipboard />
                  <Quantity quantity="upload_time" label='upload time' noWrap {...quantityProps} >
                    <Typography noWrap>
                      {new Date(calcData.upload_time).toLocaleString()}
                    </Typography>
                  </Quantity>
                  <Quantity quantity="last_processing" label='last processing' loading={loading} placeholder="not processed" noWrap {...quantityProps}>
                    <Typography noWrap>
                      {new Date(calcData.last_processing).toLocaleString()}
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

        {domain && <domain.EntryCards data={calcData} calcId={calcId} uploadId={uploadId} classes={{root: classes.entryCards}} />}
      </EntryPageContent>
    )
  }
}

export default compose(withApi(false, true), withStyles(RepoEntryView.styles))(RepoEntryView)
