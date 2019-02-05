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

  data(quantity) {
    const path = quantity.split('.')
    let data = this.state.calcData
    for (let i = 0; i < path.length; i++) {
      if (data) {
        data = data[path[i]]
      }
    }
    return data
  }

  renderQuantity(quantity, label, defaultValue) {
    const value = this.data(quantity) || defaultValue || ''

    return (
      <div key={quantity}>
        <Typography variant="caption">{label}</Typography>
        <Typography variant="body1">{value}</Typography>
      </div>
    )
  }

  render() {
    const { classes, ...calcProps } = this.props
    const { uploadId, calcId } = calcProps

    const filePaths = this.data('section_repository_info.repository_filepaths') || []
    const mainfile = this.data('section_calculation_info.main_file')
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
              {this.data('section_repository_info.section_repository_parserdata.repository_chemical_formula')}
            </CalcQuantity>
          </div>
          <div className={classes.quantityRow}>
            <CalcQuantity label='dft code'>
              {this.data('section_repository_info.section_repository_parserdata.repository_program_name')}
            </CalcQuantity>
            <CalcQuantity label='dft code version'>
              {this.data('section_repository_info.section_repository_parserdata.repository_code_version')}
            </CalcQuantity>
          </div>
          <div className={classes.quantityRow}>
            <CalcQuantity label='basis set'>
              {this.data('section_repository_info.section_repository_parserdata.repository_basis_set_type')}
            </CalcQuantity>
            <CalcQuantity label='xc functional'>
              {this.data('section_repository_info.section_repository_parserdata.repository_xc_treatment')}
            </CalcQuantity>
          </div>
          <div className={classes.quantityRow}>
            <CalcQuantity label='system type'>
              {this.data('section_repository_info.section_repository_parserdata.repository_system_type')}
            </CalcQuantity>
            <CalcQuantity label='crystal system'>
              {this.data('section_repository_info.section_repository_parserdata.repository_crystal_system')}
            </CalcQuantity>
            <CalcQuantity label='spacegroup'>
              {this.data('section_repository_info.section_repository_parserdata.repository_spacegroup_nr')}
            </CalcQuantity>
          </div>
          <div className={classes.quantityRow}>
            <CalcQuantity label='upload id'>
              {this.data('section_calculation_info.upload_id')}
            </CalcQuantity>
            <CalcQuantity label='calculation id'>
              {this.data('section_calculation_info.calc_id')}
            </CalcQuantity>
            <CalcQuantity label='mainfile'>
              {mainfile}
            </CalcQuantity>
            <CalcQuantity label='calculation hash'>
              {this.data('section_calculation_info.calc_hash')}
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
