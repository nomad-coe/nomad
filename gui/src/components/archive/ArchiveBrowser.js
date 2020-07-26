
import React, { useState, useContext, useRef, useLayoutEffect } from 'react'
import PropTypes from 'prop-types'
import { RecoilRoot, atom, useRecoilState } from 'recoil'
import { makeStyles, Card, CardContent, Box, Typography, FormGroup, FormControlLabel, Checkbox } from '@material-ui/core'
import grey from '@material-ui/core/colors/grey'
import ArrowRightIcon from '@material-ui/icons/ArrowRight'
import archiveAdaptorFactory from './archiveAdaptors'
import classNames from 'classnames'

export const viewConfigState = atom({
  key: 'viewConfig',
  default: {
    'showDescriptions': false,
    'showDefinitions': false
  }
})

export const filterConfigState = atom({
  key: 'filterConfig',
  default: {
    'showCodeSpecific': false
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

  const [lanes, setLanes] = useState([{key: 'root', adaptor: archiveAdaptorFactory(data)}])
  return <RecoilRoot>
    <ArchiveBrowserConfig />
    <Card>
      <CardContent>
        <div className={classes.root} ref={rootRef} >
          <div className={classes.lanesContainer} ref={outerRef} >
            <div className={classes.lanes} ref={innerRef} >
              {lanes.map((lane, index) => (
                <Lane
                  key={`${lane.key}:${index}`} adaptor={lane.adaptor}
                  onSetNext={next => setLanes(
                    [...lanes.slice(0, index + 1), {key: next, adaptor: lane.adaptor.itemAdaptor(next)}]
                  )}
                />
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
  const [viewConfig, setViewConfig] = useRecoilState(viewConfigState)
  const handleViewConfigChange = event => setViewConfig({
    ...viewConfig,
    [event.target.name]: event.target.checked
  })

  const [filterConfig, setFilterConfig] = useRecoilState(filterConfigState)
  const handleFilterConfigChange = event => setFilterConfig({
    ...filterConfig,
    [event.target.name]: event.target.checked
  })

  return <FormGroup row>
    <FormControlLabel
      control={
        <Checkbox
          checked={filterConfig.showCodeSpecific}
          onChange={handleFilterConfigChange}
          name="showCodeSpecific"
        />
      }
      label="include code specific"
    />
    <Box margin={2} />
    <FormControlLabel
      control={
        <Checkbox
          checked={viewConfig.showDefinitions}
          onChange={handleViewConfigChange}
          name="showDefinitions" />
      }
      label="show definitions"
    />
    <FormControlLabel
      disabled={!viewConfig.showDefinitions}
      control={
        <Checkbox
          checked={viewConfig.showDescriptions}
          onChange={handleViewConfigChange}
          name="showDescriptions" />
      }
      label="show descriptions"
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
function Lane({adaptor, onSetNext}) {
  const classes = useLaneStyles()
  const [selected, setSelected] = useState()

  const onSelect = key => {
    setSelected(key)
    onSetNext(key)
  }

  return <div className={classes.root}>
    <div className={classes.container}>
      <laneContext.Provider value={[selected, onSelect]}>
        {adaptor.render()}
      </laneContext.Provider>
    </div>
  </div>
}
Lane.propTypes = ({
  adaptor: PropTypes.object.isRequired,
  onSetNext: PropTypes.func.isRequired
})

const useItemStyles = makeStyles(theme => ({
  root: {
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
  childContainer: {
    flexGrow: 1,
    overflow: 'hidden',
    textOverflow: 'ellipsis',
    whiteSpace: 'nowrap'
  }
}))

export function Item({children, itemKey}) {
  const classes = useItemStyles()
  const [selected, setSelected] = useContext(laneContext)
  return <span
    className={classNames(
      classes.root,
      selected === itemKey ? classes.rootSelected : classes.rootUnSelected
    )}
    onClick={() => setSelected(itemKey)}
  >
    <span className={classes.childContainer}>{children}</span>
    <ArrowRightIcon/>
  </span>
}
Item.propTypes = ({
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ]).isRequired,
  itemKey: PropTypes.string.isRequired
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

export function Compartment({title, children}) {
  if (!React.Children.count(children)) {
    return ''
  }
  return <React.Fragment>
    <Box paddingTop={1}>
      {title && <Typography variant="overline">{title}</Typography>}
    </Box>
    {children}
  </React.Fragment>
}
Compartment.propTypes = ({
  title: PropTypes.string.isRequired,
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ]).isRequired
})
