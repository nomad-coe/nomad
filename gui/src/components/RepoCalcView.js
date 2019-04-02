import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Typography, Divider, Card, CardContent, Grid, CardHeader, Fab } from '@material-ui/core'
import { withApi } from './api'
import { compose } from 'recompose'
import RawFiles from './RawFiles'
import Download from './Download'
import DownloadIcon from '@material-ui/icons/CloudDownload'

function CalcQuantity(props) {
  const {children, label, typography, loading, placeholder, noWrap} = props
  const content = (!children || children.length === 0) ? null : children
  return (
    <div style={{margin: '8px 24px 0px 0'}}>
      <Typography variant="caption">{label}</Typography>
      <Typography noWrap={noWrap} variant={typography || 'body1'}>{content || <i>{loading ? 'loading...' : placeholder || 'unavailable'}</i>}</Typography>
    </div>
  )
}

CalcQuantity.propTypes = {
  classes: PropTypes.object,
  children: PropTypes.node,
  label: PropTypes.string,
  typography: PropTypes.string,
  loading: PropTypes.bool,
  placeholder: PropTypes.string,
  noWrap: PropTypes.bool
}

class RepoCalcView extends React.Component {
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
    }
  })

  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired,
    uploadId: PropTypes.string.isRequired,
    calcId: PropTypes.string.isRequired
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
    const { classes, ...calcProps } = this.props
    const calcData = this.state.calcData || calcProps
    const loading = !this.state.calcData
    const { uploadId, calcId } = calcProps

    const filePaths = calcData.files || []
    const mainfile = calcData.mainfile
    const calcPath = mainfile ? mainfile.substring(0, mainfile.lastIndexOf('/')) : null

    const authors = loading ? null : calcData.authors

    return (
      <div className={classes.root}>
        <div className={classes.content}>

          <div className={classes.title}>
            <CalcQuantity label="chemical formula" typography="h3" loading={loading}>
              {calcData.formula}
            </CalcQuantity>
          </div>

          <Grid container spacing={24}>
            <Grid item xs={7}>
              <Card>
                <CardHeader title="Metadata" />
                <CardContent classes={{root: classes.cardContent}}>
                  <div className={classes.quantityColumn}>
                    <div className={classes.quantityRow}>
                      <CalcQuantity label='dft code' loading={loading}>
                        {calcData.code_name}
                      </CalcQuantity>
                      <CalcQuantity label='dft code version' loading={loading}>
                        {calcData.code_version}
                      </CalcQuantity>
                    </div>
                    <div className={classes.quantityRow}>
                      <CalcQuantity label='basis set' loading={loading}>
                        {calcData.basis_set}
                      </CalcQuantity>
                      <CalcQuantity label='xc functional' loading={loading}>
                        {calcData.xc_functional}
                      </CalcQuantity>
                    </div>
                    <div className={classes.quantityRow}>
                      <CalcQuantity label='system type' loading={loading}>
                        {calcData.system}
                      </CalcQuantity>
                      <CalcQuantity label='crystal system' loading={loading}>
                        {calcData.crystal_system}
                      </CalcQuantity>
                      <CalcQuantity label='spacegroup' loading={loading}>
                        {calcData.spacegroup_symbol} ({calcData.spacegroup})
                      </CalcQuantity>
                    </div>
                  </div>
                </CardContent>
                <Divider />
                <CardContent classes={{root: classes.cardContent}}>
                  <div className={classes.quantityColumn}>
                    <div className={classes.quantityColumn}>
                      <CalcQuantity label='comment' loading={loading} placeholder='no comment'>
                        {calcData.comment}
                      </CalcQuantity>
                      <CalcQuantity label='references' loading={loading} placeholder='no references'>
                        {calcData.references ? calcData.references.map(ref => (<a key={ref.id} href={ref.value}>{ref.value}</a>)) : null}
                      </CalcQuantity>
                      <CalcQuantity label='authors' loading={loading}>
                        {authors ? authors.map(author => author.name) : null}
                      </CalcQuantity>
                      <CalcQuantity label='datasets' loading={loading} placeholder='no datasets'>
                        {calcData.datasets ? calcData.datasets.map(ds => ds.name).join(', ') : null}
                      </CalcQuantity>
                    </div>
                  </div>
                </CardContent>
              </Card>
            </Grid>

            <Grid item xs={5}>
              <Card>
                <CardHeader title="Identification" />
                <CardContent classes={{root: classes.cardContent}}>
                  <div className={classes.quantityColumn} style={{maxWidth: 350}}>
                    <CalcQuantity label='PID' loading={loading} noWrap>
                      <b>{calcData.pid}</b>
                    </CalcQuantity>
                    <CalcQuantity label='upload id' noWrap>
                      {calcData.upload_id}
                    </CalcQuantity>
                    <CalcQuantity label='calculation id' noWrap>
                      {calcData.calc_id}
                    </CalcQuantity>
                    <CalcQuantity label='mainfile' loading={loading} noWrap>
                      {mainfile}
                    </CalcQuantity>
                    <CalcQuantity label='calculation hash' loading={loading} noWrap>
                      {calcData.calc_hash}
                    </CalcQuantity>
                  </div>
                </CardContent>
              </Card>
            </Grid>

            <Grid item xs={12}>
              <Card>
                <CardHeader title="Raw files" />
                <CardContent classes={{root: classes.cardContent}}>
                  <RawFiles {...calcProps} files={filePaths} />
                </CardContent>
              </Card>
            </Grid>
          </Grid>

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

export default compose(withApi(false, true), withStyles(RepoCalcView.styles))(RepoCalcView)
