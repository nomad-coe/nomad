
import React, { useContext, useState, useEffect, useRef } from 'react'
import PropTypes from 'prop-types'
import { matchPath, useLocation, useHistory, useRouteMatch } from 'react-router-dom'
import Viewer from './Viewer'
import { apiContext } from '../api'
import MetainfoSearch from './MetainfoSearch'
import { FormControl, Select, Input, MenuItem, ListItemText, InputLabel, makeStyles } from '@material-ui/core'
import { schema } from '../MetaInfoRepository'
import { errorContext } from '../errors'

export const help = `
The NOMAD *metainfo* defines all quantities used to represent archive data in
NOMAD. You could say it is the archive *schema*. You can browse this schema and
all its definitions here.

The NOMAD metainfo contains three different *kinds* of definitions:

- **sections**: A section is a nested groups of quantities that allow a hierarchical data structure
- **values**: Actual quantities that contain data
- **references**: References that allow to connect related sections.

All definitions have a name that you can search for. Furthermore, all definitions
are organized in packages. There is a *common* pkg with definitions that are
used by all codes and there are packages for each code with code specific definitions.
You can select the pkg to browse below.

Depending on the selected pkg, there are quite a large number of definitions.
You can use the *definition* field to search based on definition names.

All definitions are represented as *cards* below. Click on the various card items
to expand sub-sections, open values or references, hide and show compartments, or
collapse cards again. The highlighted *main* card cannot be collapsed. The
shapes in the background represent section containment (grey) and
reference (blue) relations.

If you bookmark this page, you can save the definition represented by the highlighted
*main* card.

To learn more about the meta-info, visit the [meta-info homepage](https://metainfo.nomad-coe.eu/nomadmetainfo_public/archive.html).
`
const MenuProps = {
  PaperProps: {
    style: {
      width: 300, maxHeight: '90vh'
    }
  }
}

const useStyles = makeStyles(theme => ({
  root: {},
  forms: {
    padding: `${theme.spacing(3)}px ${theme.spacing(3)}px 0 ${theme.spacing(3)}px`
  },
  packageSelect: {
    width: 300, height: 24
  },
  search: {
    width: 450,
    marginRight: theme.spacing(2)
  }
}))

export default function MetaInfoBrowser({visible}) {
  const classes = useStyles()

  const routingRef = useRef({})
  const routing = {
    location: useLocation(),
    history: useHistory(),
    routeMatch: useRouteMatch()
  }
  if (visible) {
    routingRef.current = routing
  }
  const {location, history, routeMatch} = routingRef.current

  const match = matchPath(location.pathname, {
    path: `${routeMatch.path}/:pkg?/:metainfo?`
  })
  const pkg = match.params.pkg || 'general'
  const metainfoName = match.params.metainfo || 'section_run'

  const {api, loading} = useContext(apiContext)
  const {raiseError} = useContext(errorContext)

  const [metainfos, setMetainfos] = useState(null)
  const [packages, setPackages] = useState(null)

  useEffect(() => {
    api.getInfo().then(info => {
      setPackages(info.metainfo_packages)
    }).catch(raiseError)
  }, [api, raiseError])

  useEffect(() => {
    api.getMetaInfo(pkg).then(metainfos => {
      const definition = metainfos.get(metainfoName)
      if (!definition) {
        history.push(`/metainfo/${pkg}/section_run`)
      } else {
        setMetainfos(metainfos)
      }
    }).catch(raiseError)
  }, [pkg, metainfoName, api, history, raiseError])

  const handleSelectedPackageChanged = event => {
    history.push(`/metainfo/${event.target.value}/section_run`)
  }

  const handleSearch = metainfoName => {
    if (metainfos.get(metainfoName)) {
      history.push(`/metainfo/${pkg}/${metainfoName}`)
    }
  }

  if (!metainfos || !packages) {
    return <div />
  }

  const metainfo = metainfos.resolve(metainfos.createProxy(metainfoName))
  console.log(metainfoName)

  return <div style={{display: visible ? 'block' : 'none'}}>
    <div className={classes.forms}>
      <form style={{ display: 'flex' }}>
        <MetainfoSearch
          classes={{container: classes.search}}
          suggestions={Object.values(metainfos.names).filter(metainfo => !schema.isPackage(metainfo))}
          onChange={handleSearch}
        />
        <FormControl disabled={loading > 0}>
          <InputLabel htmlFor="select-multiple-checkbox">Package</InputLabel>
          <Select
            classes={{root: classes.packageSelect}}
            value={pkg}
            onChange={handleSelectedPackageChanged}
            input={<Input id="select-multiple-checkbox" />}
            MenuProps={MenuProps}
          >
            {packages
              .map(name => {
                return <MenuItem key={name} value={name}>
                  <ListItemText primary={name} style={{margin: 0}} />
                </MenuItem>
              })
            }
          </Select>
        </FormControl>
      </form>
    </div>
    <Viewer key={`${metainfo.package.name}/${metainfo.name}`} rootElement={metainfo} packages={metainfos.contents} />
  </div>
}
MetaInfoBrowser.propTypes = {
  visible: PropTypes.bool
}
