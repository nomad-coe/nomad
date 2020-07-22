
import React, { useEffect, useState, useContext, useCallback, useRef, useMemo, useLayoutEffect } from 'react'
import { makeStyles, Typography, Box } from '@material-ui/core'
import grey from '@material-ui/core/colors/grey';
import ArrowRightIcon from '@material-ui/icons/ArrowRight';
import { json } from 'd3';

const exampleData = {
  "section_one": {
    "value1": "1.1",
    "value2": "1.2"
  },
  "section_two": [
    {
      "value1": "2.1",
      "value2": "2.2"
    },
    {
      "value1": "3.1",
      "value2": "3.2"
    }
  ],
  "value1": "1.1",
  "value2": "1.2"
}

function jsonAdaptorFactory(child) {
  if (Array.isArray(child)) {
    if (child.length === 1) {
      return new ObjectAdaptor(child[0])
    }
    return new ArrayAdaptor(child)
  } else if (typeof child === 'string') {
    return new ValueAdaptor(child)
  } else if (typeof child === 'object') {
    return new ObjectAdaptor(child)
  } else {
    return new ValueAdaptor(child)
  }
}

class Adaptor {
  constructor(e) {
    this.e = e
  }

  itemAdaptor(key) {
    return null
  }
}

class ObjectAdaptor extends Adaptor {
  itemAdaptor(key) {
    return jsonAdaptorFactory(this.e[key])
  }
  render() {
    return <React.Fragment>
      {Object.keys(this.e).map(key => (
        <Item key={key} itemKey={key}>
          <Typography>
            {key}
          </Typography>
        </Item>
      ))}
    </React.Fragment>
  }
}

class ValueAdaptor extends Adaptor {
  render() {
    return <Content>
      <Typography>{String(this.e)}</Typography>
    </Content>
  }
}

class ArrayAdaptor extends Adaptor {
  itemAdaptor(index) {
    return jsonAdaptorFactory(this.e[index])
  }
  render() {
    return <React.Fragment>
      {this.e.map((_, index) => (
        <Item key={index} itemKey={index}>
          <Typography>
            {index + 1}
          </Typography>
        </Item>
      ))}
    </React.Fragment>
  }
}

const useBrowserStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
    flexFlow: 'column'
  },
  lanesContainer: {
    flex: '1 1 auto',
    height: '100%',
    overflowX: 'auto',
    scrollBehavior: 'smooth'
  },
  lanes: {
    display: 'table',
    overflow: 'scroll',
    height: '100%',
    overflowY: 'hidden',
    width: 'fit-content',
    margin: `${theme.spacing(1)}px 0`
  }
}))

export default function ArchiveBrowser({data}) {
  data = data || exampleData

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

  const [lanes, setLanes] = useState([{key: 'root', adaptor: jsonAdaptorFactory(data)}])
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
    padding: `0 0 0 ${theme.spacing(1)}px`,
    '&:hover': {
      backgroundColor: grey[300]
    }
  },
  selected: {
    padding: `0 0 0 ${theme.spacing(1)}px`,
    backgroundColor: theme.palette.primary.main,
    color: theme.palette.primary.contrastText
  }
}))
function Item({children, itemKey}) {
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

const useContentStyles = makeStyles(theme => ({
  root: {
    padding: `0 ${theme.spacing(1)}px`
  }
}))
function Content({children}) {
  const classes = useContentStyles()
  return <div className={classes.root}>
    {children}
  </div>
}