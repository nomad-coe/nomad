import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Typography, Divider, LinearProgress, Fab } from '@material-ui/core'
import { withApi } from './api'
import { compose } from 'recompose'
import RawFiles from './RawFiles'
import Download from './Download'
import DownloadIcon from '@material-ui/icons/CloudDownload'
import { withErrors } from './errors'

function CalcQuantity(props) {
  const {children, label, typography} = props
  return (
    <div style={{margin: '0px 24px 8px 0'}}>
      <Typography variant="caption">{label}</Typography>
      <Typography variant={typography || 'body1'}>{children || 'loading...'}</Typography>
    </div>
  )
}

CalcQuantity.propTypes = {
  classes: PropTypes.object,
  children: PropTypes.node,
  label: PropTypes.string,
  typography: PropTypes.string
}

class RepoCalcView extends React.Component {
  static styles = theme => ({
    root: {},
    content: {
      marginTop: theme.spacing.unit * 3
    },
    quantityRow: {
      display: 'flex',
      flexDirection: 'row',
      flexWrap: 'wrap',
      marginBottom: theme.spacing.unit
    },
    downloadFab: {
      position: 'absolute',
      zIndex: 1,
      top: theme.spacing.unit,
      right: theme.spacing.unit * 3
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
    const { uploadId, calcId } = calcProps
    const calcData = this.state.calcData || {}

    const filePaths = calcData.files || []
    const mainfile = calcData.mainfile
    const calcPath = mainfile ? mainfile.substring(0, mainfile.lastIndexOf('/')) : null

    return (
      <div className={classes.root}>
        {!this.state.calcData ? <LinearProgress /> : ''}
        <div className={classes.content}>
          <Download
            disabled={!mainfile} tooltip="download all raw files for calculation"
            classes={{root: classes.downloadFab}}
            component={Fab} className={classes.downloadFab} color="primary" size="medium"
            url={`raw/${uploadId}/${calcPath}/*`} fileName={`${calcId}.zip`}
          >
            <DownloadIcon />
          </Download>
          <div className={classes.quantityRow}>
            <CalcQuantity label="chemical formula" typography="h4">
              {calcData.formula}
            </CalcQuantity>
          </div>
          <div className={classes.quantityRow}>
            <CalcQuantity label='dft code'>
              {calcData.code_name}
            </CalcQuantity>
            <CalcQuantity label='dft code version'>
              {calcData.code_version}
            </CalcQuantity>
          </div>
          <div className={classes.quantityRow}>
            <CalcQuantity label='basis set'>
              {calcData.basis_set}
            </CalcQuantity>
            <CalcQuantity label='xc functional'>
              {calcData.xc_functional}
            </CalcQuantity>
          </div>
          <div className={classes.quantityRow}>
            <CalcQuantity label='system type'>
              {calcData.system}
            </CalcQuantity>
            <CalcQuantity label='crystal system'>
              {calcData.crystal_system}
            </CalcQuantity>
            <CalcQuantity label='spacegroup'>
              {calcData.spacegroup}
            </CalcQuantity>
          </div>
          <div className={classes.quantityRow}>
            <CalcQuantity label='upload id'>
              {calcData.upload_id}
            </CalcQuantity>
            <CalcQuantity label='calculation id'>
              {calcData.calc_id}
            </CalcQuantity>
            <CalcQuantity label='mainfile'>
              {mainfile}
            </CalcQuantity>
            <CalcQuantity label='calculation hash'>
              {calcData.calc_hash}
            </CalcQuantity>
          </div>
          <Divider />
          <RawFiles {...calcProps} files={filePaths} />
        </div>
      </div>
    )
  }
}

export default compose(withApi(false), withErrors, withStyles(RepoCalcView.styles))(RepoCalcView)
