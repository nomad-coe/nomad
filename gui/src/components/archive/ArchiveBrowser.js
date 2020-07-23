
import React, { useState, useContext, useRef, useLayoutEffect } from 'react'
import { makeStyles, Typography, Box } from '@material-ui/core'
import grey from '@material-ui/core/colors/grey'
import ArrowRightIcon from '@material-ui/icons/ArrowRight'
import { SectionAdaptor, entryArchiveDef } from './archive'

const useBrowserStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
    flexFlow: 'column'
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

  const contextData = {
    archive: data
  }

  const [lanes, setLanes] = useState([{key: 'root', adaptor: new SectionAdaptor(data, entryArchiveDef, contextData)}])
  return (
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
  )
}

const laneContext = React.createContext()
const useLaneStyles = makeStyles(theme => ({
  root: {
    minWidth: 200,
    maxWidth: 512,
    borderRight: `solid 1px ${grey[500]}`,
    display: 'table-cell',
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

const useItemStyles = makeStyles(theme => ({
  root: {
    margin: `0 -${theme.spacing(1)}px`,
    padding: `0 0 0 ${theme.spacing(1)}px`,
    '&:hover': {
      backgroundColor: grey[300]
    },
    whiteSpace: 'nowrap'
  },
  selected: {
    margin: `0 -${theme.spacing(1)}px`,
    padding: `0 0 0 ${theme.spacing(1)}px`,
    backgroundColor: theme.palette.primary.main,
    color: theme.palette.primary.contrastText,
    whiteSpace: 'nowrap'
  }
}))
export function Item({children, itemKey}) {
  const classes = useItemStyles()
  const [selected, setSelected] = useContext(laneContext)
  return <Box
    display="flex" classes={{root: selected === itemKey ? classes.selected : classes.root}}
    onClick={() => setSelected(itemKey)}
  >
    <Box flexGrow={1}>
      {children}
    </Box>
    <ArrowRightIcon/>
  </Box>
}

export function Content({children}) {
  return <Box padding={1}>
    {children}
  </Box>
}
