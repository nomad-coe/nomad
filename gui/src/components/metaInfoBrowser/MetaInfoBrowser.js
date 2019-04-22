
import React, { Component } from 'react'
import { withRouter } from 'react-router-dom'
import Viewer from './Viewer'
import PropTypes from 'prop-types'
import { withApi } from '../api'
import { Help } from '../help'
import MetainfoSearch from './MetainfoSearch'
import { FormControl, withStyles, Select, Input, MenuItem, ListItemText, InputLabel, FormGroup } from '@material-ui/core'
import { compose } from 'recompose'

const ITEM_HEIGHT = 48
const ITEM_PADDING_TOP = 8
const MenuProps = {
  PaperProps: {
    style: {
      maxHeight: ITEM_HEIGHT * 4.5 + ITEM_PADDING_TOP,
      width: 300
    }
  }
}

class MetaInfoBrowser extends Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    metainfo: PropTypes.string,
    api: PropTypes.object.isRequired,
    loading: PropTypes.number,
    raiseError: PropTypes.func.isRequired,
    history: PropTypes.object.isRequired
  }

  static styles = theme => ({
    root: {},
    forms: {
      padding: `${theme.spacing.unit * 3}px ${theme.spacing.unit * 3}px 0 ${theme.spacing.unit * 3}px`
    },
    packageSelect: {
      width: 300
    },
    search: {
      marginLeft: theme.spacing.unit * 3
    }
  })

  state = {
    metainfos: null,
    allMetainfos: null,
    selectedPackage: null,
    loadedPackage: null
  }

  constructor(props) {
    super(props)
    this.handleSelectedPackageChanged = this.handleSelectedPackageChanged.bind(this)
    this.handleSearch = this.handleSearch.bind(this)
  }

  update(pkg) {
    this.props.api.getMetaInfo(pkg).then(metainfos => {
      const metainfoName = this.props.metainfo || 'section_run'
      const definition = metainfos.get(metainfoName)
      if (!definition) {
        this.props.history.push('/metainfo/section_run')
      } else {
        this.setState({loadedPackage: pkg, metainfos: metainfos})
      }
    }).catch(error => {
      this.props.raiseError(error)
    })
  }

  init() {
    this.props.api.getMetaInfo('all.nomadmetainfo.json').then(metainfos => {
      const metainfoName = this.props.metainfo || 'section_run'
      const definition = metainfos.get(metainfoName)
      this.setState({allMetainfos: metainfos, selectedPackage: definition.package.name})
      this.update(definition.package.name)
    }).catch(error => {
      this.props.raiseError(error)
    })
  }

  componentDidUpdate(prevProps) {
    if (this.props.metainfo !== prevProps.metainfo) {
      this.init()
    }
  }

  componentDidMount() {
    this.init()
  }

  handleSelectedPackageChanged(event) {
    this.setState({selectedPackage: event.target.value})
    this.update(event.target.value)
  }

  handleSearch(term) {
    if (this.state.metainfos.get(term)) {
      this.props.history.push(`/metainfo/${term}`)
    }
  }

  render() {
    const { classes, loading } = this.props
    const { metainfos, selectedPackage, allMetainfos, loadedPackage } = this.state

    if (!metainfos || !allMetainfos) {
      return <div />
    }

    const metainfoName = this.props.metainfo || 'section_run'
    const metainfo = metainfos.resolve(metainfos.createProxy(metainfoName))

    return <div>
      <div className={classes.forms}>
        <Help cookie="uploadList">{`
          The nomad *metainfo* is the schema that nomad uses to define the quantities that
          comprise all nomad data. The nomad metainfo is consist of definitions for:

          - **sections**: A section are nested groups of quantities that allow a hierarchical data structure
          - **values**: Actual quantities that contain data
          - **references**: References that allow to connect related sections.

          All definitions have a name that you can search for. Furthermore, all definitions
          are organized in packages. There is a *common* package with definitions that are
          used by all codes and there are packages for each code with code specific definitions.
        `}</Help>
        <form style={{ display: 'flex' }}>
          <FormControl disabled={loading > 0}>
            <InputLabel htmlFor="select-multiple-checkbox">Package</InputLabel>
            <Select
              className={classes.packageSelect}
              value={selectedPackage}
              onChange={this.handleSelectedPackageChanged}
              input={<Input id="select-multiple-checkbox" />}
              MenuProps={MenuProps}
            >
              {allMetainfos.contents
                .map(pkg => pkg.name)
                .map(name => {
                  return <MenuItem key={name} value={name}>
                    <ListItemText primary={name.substring(0, name.length - 19)} />
                  </MenuItem>
                })
              }
            </Select>
          </FormControl>
          <MetainfoSearch classes={{root: classes.search}}
            suggestions={Object.values(metainfos.names)}
            onChange={this.handleSearch}
          />
        </form>
      </div>
      <Viewer key={loadedPackage} rootElement={metainfo} packages={metainfos.contents} />
    </div>
  }
}

export default compose(withRouter, withApi(false), withStyles(MetaInfoBrowser.styles))(MetaInfoBrowser)
