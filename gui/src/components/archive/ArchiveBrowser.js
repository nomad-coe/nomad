
import React, { useEffect, useState, useContext, useCallback, useRef, useMemo, useLayoutEffect } from 'react'
import { makeStyles, Typography, Box } from '@material-ui/core'
import grey from '@material-ui/core/colors/grey';
import ArrowRightIcon from '@material-ui/icons/ArrowRight';

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

class Adaptor {
  constructor(e) {
    this.e = e
  }

  hasItems() {
    return false
  }

  item(key) {
    return null
  }

  itemAdaptor(key) {
    return null
  }

  items() {
    return []
  }

  renderItem(key) {
    return null
  }

  renderContent() {
    return null
  }
}

class ObjectAdaptor extends Adaptor {
  hasItems() {
    return  this.e.length > 0
  }

  item(key) {
    return this.e[key]
  }

  itemAdaptor(key) {
    const child = this.e[key]
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

  items() {
    return Object.keys(this.e)
  }

  renderItem(key) {
    return <Typography>{key}</Typography>
  }

  renderContent() {
    return null
  }
}

class ValueAdaptor extends Adaptor {
  renderContent() {
    return <Typography>{String(this.e)}</Typography>
  }
}

class ArrayAdaptor extends Adaptor {
  renderContent() {
    return <Typography>ARRAY</Typography>
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

  const [lanes, setLanes] = useState([{key: 'root', adaptor: new ObjectAdaptor(data)}])
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

  return <div className={classes.root}>
    <div className={classes.container}>
      {adaptor.items().map(itemKey => (
        <Item
          key={itemKey}
          selected={selected === itemKey}
          onClick={() => {
            setSelected(itemKey)
            onSetNext(itemKey)
          }}
        >
          {adaptor.renderItem(itemKey)}
        </Item>
      ))}&nbsp;
    </div>
  </div>
}

const useItemStyles = makeStyles(theme => ({
  root: {
    padding: `0 ${theme.spacing(1)}px`,
    '&:hover': {
      backgroundColor: grey[300]
    }
  },
  selected: {
    padding: `0 ${theme.spacing(1)}px`,
    backgroundColor: theme.palette.primary.main,
    color: theme.palette.primary.contrastText
  }
}))
function Item({children, selected, hasItems, ...rest}) {
  const classes = useItemStyles()
  return <Box
    display="flex" classes={{root: selected ? classes.selected : classes.root}}
    {...rest}
  >
    <Box flexGrow={1}>
      {children}
    </Box>
    <ArrowRightIcon/>
  </Box>
}