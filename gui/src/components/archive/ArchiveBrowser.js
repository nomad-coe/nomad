
import React, { useContext, useRef, useLayoutEffect, useMemo, useState } from 'react'
import PropTypes from 'prop-types'
import { RecoilRoot, atom, useRecoilState } from 'recoil'
import { makeStyles, Card, CardContent, Box, Typography, FormGroup, FormControlLabel, Checkbox } from '@material-ui/core'
import grey from '@material-ui/core/colors/grey'
import ArrowRightIcon from '@material-ui/icons/ArrowRight'
import ArrowDownIcon from '@material-ui/icons/ArrowDropDown'
import archiveAdaptorFactory from './archiveAdaptors'
import classNames from 'classnames'
import { useLocation, useRouteMatch, Link } from 'react-router-dom'

export const configState = atom({
  key: 'config',
  default: {
    'showMeta': true,
    'showCodeSpecific': false,
    'showAllDefined': true
  }
})

const useBrowserStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
    flexFlow: 'column',
    margin: `-${theme.spacing(2)}px`,
    marginBottom: `-${theme.spacing(3)}px`
  },
  lanesContainer: {
    flex: '1 1 auto',
    height: '100%',
    overflowX: 'auto',
    overflowY: 'hidden',
    scrollBehavior: 'smooth'
  },
  lanes: {
    display: 'table',
    overflow: 'scroll',
    height: '100%',
    overflowY: 'hidden',
    width: 'fit-content'
  }
}))

export default function ArchiveBrowser({data}) {
  const classes = useBrowserStyles()
  const rootRef = useRef()
  const outerRef = useRef()
  const innerRef = useRef()

  useLayoutEffect(() => {
    const height = window.innerHeight - outerRef.current.getBoundingClientRect().top - 24
    rootRef.current.style.height = `${height}px`
    const scrollAmmount = innerRef.current.clientWidth - outerRef.current.clientWidth
    outerRef.current.scrollLeft = Math.max(scrollAmmount, 0)
  })

  const { pathname } = useLocation()
  const { url } = useRouteMatch()
  const archivePath = pathname.substring(url.length).split('/')

  const root = useMemo(() => ({
    key: 'root',
    path: url.endsWith('/') ? url.substring(0, url.length - 1) : url,
    adaptor: archiveAdaptorFactory(data),
    next: null
  }), [data, url])

  const lanes = [root]
  archivePath.filter(segment => segment).forEach((segment, i) => {
    const prev = lanes[i]
    const lane = {
      key: segment,
      path: `${prev.path}/${segment}`,
      adaptor: prev.adaptor.itemAdaptor(segment)
    }
    lanes.push(lane)
    prev.next = lane
  })

  return <RecoilRoot>
    <ArchiveBrowserConfig />
    <Card>
      <CardContent>
        <div className={classes.root} ref={rootRef} >
          <div className={classes.lanesContainer} ref={outerRef} >
            <div className={classes.lanes} ref={innerRef} >
              {lanes.map((lane, index) => (
                <Lane key={index} lane={lane} />
              ))}
            </div>
          </div>
        </div>
      </CardContent>
    </Card>
  </RecoilRoot>
}
ArchiveBrowser.propTypes = ({
  data: PropTypes.object.isRequired
})

function ArchiveBrowserConfig() {
  const [config, setConfig] = useRecoilState(configState)
  const handleConfigChange = event => {
    const changes = {[event.target.name]: event.target.checked}
    if (changes.showCodeSpecific) {
      changes.showAllDefined = !changes.showCodeSpecific
    } else if (changes.showAllDefined) {
      changes.showCodeSpecific = !changes.showAllDefined
    }
    setConfig({...config, ...changes})
  }

  return <FormGroup row>
    <FormControlLabel
      control={
        <Checkbox
          checked={config.showCodeSpecific}
          onChange={handleConfigChange}
          name="showCodeSpecific"
        />
      }
      label="include code specific"
    />
    <FormControlLabel
      control={
        <Checkbox
          checked={config.showAllDefined}
          onChange={handleConfigChange}
          name="showAllDefined"
        />
      }
      label="show all defined metadata"
    />
    <FormControlLabel
      control={
        <Checkbox
          checked={config.showMeta}
          onChange={handleConfigChange}
          name="showMeta" />
      }
      label="show metainfo definitions"
    />
  </FormGroup>
}

const laneContext = React.createContext()
const useLaneStyles = makeStyles(theme => ({
  root: {
    minWidth: 200,
    maxWidth: 512,
    borderRight: `solid 1px ${grey[500]}`,
    display: 'table-cell'
  },
  container: {
    display: 'block',
    height: '100%',
    overflowY: 'scroll'
  }
}))
function Lane({lane}) {
  const classes = useLaneStyles()
  const { adaptor } = lane

  return <div className={classes.root}>
    <div className={classes.container}>
      <laneContext.Provider value={lane}>
        {adaptor.render()}
      </laneContext.Provider>
    </div>
  </div>
}
Lane.propTypes = ({
  lane: PropTypes.object.isRequired
})

const useItemStyles = makeStyles(theme => ({
  root: {
    color: theme.palette.text.primary,
    textDecoration: 'none',
    margin: `0 -${theme.spacing(1)}px`,
    padding: `0 0 0 ${theme.spacing(1)}px`,
    whiteSpace: 'nowrap',
    display: 'flex'
  },
  rootSelected: {
    backgroundColor: theme.palette.primary.main,
    color: theme.palette.primary.contrastText,
    whiteSpace: 'nowrap'
  },
  rootUnSelected: {
    '&:hover': {
      backgroundColor: grey[300]
    }
  },
  disabled: {
    color: 'grey'
  },
  childContainer: {
    flexGrow: 1,
    overflow: 'hidden',
    textOverflow: 'ellipsis',
    whiteSpace: 'nowrap'
  }
}))

export function Item({children, itemKey, disabled}) {
  const classes = useItemStyles()
  const lane = useContext(laneContext)
  const selected = lane.next && lane.next.key
  if (disabled) {
    return <div className={classNames(classes.childContainer, classes.disabled)}>{children}</div>
  }
  return <Link
    className={classNames(
      classes.root,
      selected === itemKey ? classes.rootSelected : classes.rootUnSelected
    )}
    to={`${lane.path}/${itemKey}`}
  >
    <span className={classes.childContainer}>{children}</span>
    <ArrowRightIcon/>
  </Link>
}
Item.propTypes = ({
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ]).isRequired,
  itemKey: PropTypes.string.isRequired,
  disabled: PropTypes.bool
})

const useListStyles = makeStyles(theme => ({
  title: {
    color: theme.palette.text.primary,
    textDecoration: 'none',
    margin: `0 -${theme.spacing(1)}px`,
    whiteSpace: 'nowrap',
    display: 'flex',
    fontWeight: 'bold'
  },
  selected: {
    backgroundColor: theme.palette.primary.main,
    color: theme.palette.primary.contrastText,
    whiteSpace: 'nowrap'
  },
  unSelected: {
    '&:hover': {
      backgroundColor: grey[300]
    }
  }
}))
export function List({title, itemKey}) {
  const classes = useListStyles()
  const [open, setOpen] = useState(false)
  const lane = useContext(laneContext)
  const selected = lane.next && lane.next.key
  const values = lane.adaptor.e[itemKey]
  return <div>
    <Typography onClick={() => setOpen(!open)} className={classNames(
      classes.title,
      (!open && selected && selected.startsWith(itemKey + ':')) ? classes.selected : classes.unSelected
    )}>
      {open ? <ArrowDownIcon/> : <ArrowRightIcon/>}
      <span>{title || 'list'}</span>
    </Typography>
    {open &&
      <div>
        {values.map((_, index) => (
          <Item key={index} itemKey={`${itemKey}:${index}`}>
            <Box component="span" marginLeft={2}>
              <Typography component="span">{index}</Typography>
            </Box>
          </Item>
        ))}
      </div>
    }
  </div>
}
List.propTypes = ({
  itemKey: PropTypes.string.isRequired,
  title: PropTypes.string
})

export function Content({children}) {
  return <Box padding={1}>
    {children}
  </Box>
}
Content.propTypes = ({
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ]).isRequired
})

export function Compartment({title, children, color}) {
  if (!React.Children.count(children)) {
    return ''
  }
  return <React.Fragment>
    <Box paddingTop={1}>
      {title && <Typography color={color} variant="overline">{title}</Typography>}
    </Box>
    {children}
  </React.Fragment>
}
Compartment.propTypes = ({
  title: PropTypes.string,
  color: PropTypes.string,
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ]).isRequired
})
