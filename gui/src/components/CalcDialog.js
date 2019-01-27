import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Dialog, DialogContent, DialogActions, Button, DialogTitle, Tab, Tabs,
  Typography, FormGroup, FormControlLabel, Checkbox, Divider, FormLabel, IconButton,
  LinearProgress } from '@material-ui/core'
import SwipeableViews from 'react-swipeable-views'
import ArchiveCalcView from './ArchiveCalcView'
import DownloadIcon from '@material-ui/icons/CloudDownload'
import ArchiveLogView from './ArchiveLogView'
import { withApi } from './api'
import { compose } from 'recompose'

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

class CalcDialog extends React.Component {
  static styles = theme => ({
    dialog: {

    },
    dialogTitle: {
      padding: 0
    },
    dialogContent: {
      paddingBottom: 0
    },
    tabContent: {
      paddingTop: theme.spacing.unit * 3,
      overflowY: 'auto',
      height: '70vh',
      zIndex: 1
    },
    formLabel: {
      padding: theme.spacing.unit * 2
    },
    quantityRow: {
      display: 'flex',
      flexDirection: 'row',
      flexWrap: 'wrap',
      marginBottom: theme.spacing.unit
    }
  })

  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired,
    uploadId: PropTypes.string.isRequired,
    calcId: PropTypes.string.isRequired,
    onClose: PropTypes.func.isRequired
  }

  state = {
    open: false,
    calcData: null,
    viewIndex: 0
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
    const { classes, onClose, ...calcProps } = this.props
    const { viewIndex } = this.state

    const filePaths = this.data('section_repository_info.repository_filepaths') || []
    const files = filePaths.map(filePath => filePath.substring(filePath.lastIndexOf('/') + 1))

    return (
      <Dialog className={classes.dialog} open={true} onClose={onClose} fullWidth={true} maxWidth={'md'} >
        <DialogTitle disableTypography classes={{root: classes.dialogTitle}}>
          {(!this.state.calcData) ? <LinearProgress /> : ''}
          <Tabs
            className={classes.tabs}
            value={viewIndex}
            onChange={(event, state) => this.setState({viewIndex: state})}
            indicatorColor="primary"
            textColor="primary"
            variant="fullWidth"
          >
            <Tab label="Raw data" />
            <Tab label="Archive" />
            <Tab label="Logs" />
          </Tabs>
        </DialogTitle>
        <DialogContent classes={{root: classes.dialogContent}}>
          <SwipeableViews
            // axis={theme.direction === 'rtl' ? 'x-reverse' : 'x'}
            index={viewIndex}
            onChangeIndex={() => null}
          >
            <div className={classes.tabContent}>
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
                  {this.data('section_calculation_info.main_file')}
                </CalcQuantity>
                <CalcQuantity label='calculation hash'>
                  {this.data('section_calculation_info.calc_hash')}
                </CalcQuantity>
              </div>
              <Divider />
              <FormGroup row>
                <FormControlLabel label="select all" control={<Checkbox checked={false} value="select_all" />} style={{flexGrow: 1}}/>
                <FormLabel className={classes.formLabel}>0/10 files selected</FormLabel>
                <IconButton><DownloadIcon /></IconButton>
              </FormGroup>
              <Divider />
              <FormGroup row>
                {files.map((file, index) => (
                  <FormControlLabel key={index} label={file}
                    control={<Checkbox checked={false} onChange={() => true} value={file} />}
                  />
                ))}
              </FormGroup>
            </div>
            <div className={classes.tabContent}>
              <ArchiveCalcView {...calcProps} />
            </div>
            <div className={classes.tabContent}>
              <ArchiveLogView {...calcProps} />
            </div>
          </SwipeableViews>
        </DialogContent>
        <DialogActions>
          <Button onClick={onClose} color="primary" autoFocus>
              Close
          </Button>
        </DialogActions>
      </Dialog>
    )
  }
}

export default compose(withApi(false), withStyles(CalcDialog.styles))(CalcDialog)
