/*
 * Copyright The NOMAD Authors.
 *
 * This file is part of NOMAD. See https://nomad-lab.eu for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
import React, { useRef, useEffect } from 'react'
import PropTypes from 'prop-types'
import { IconButton, Card, CardActions, CardHeader, makeStyles } from '@material-ui/core'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import clsx from 'clsx'

const useStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
    flexDirection: 'column',
    transition: theme.transitions.create('height', {
      duration: theme.transitions.duration.shortest
    })
  },
  cardHeader: {
    flex: '0 0 1.5rem',
    paddingBottom: theme.spacing(1.5),
    height: '3rem'
  },
  cardContent: {
    padding: theme.spacing(2),
    paddingBottom: 0,
    paddingTop: 0,
    position: 'relative',
    display: 'flex',
    flexDirection: 'column',
    flex: '0 1 auto',
    overflow: 'hidden'
  },
  cardFixedContent: {
    padding: theme.spacing(2),
    paddingTop: 0,
    paddingBottom: 0,
    flex: '0 0 auto'
  },
  cardActions: {
    flex: '0 0 1.5rem'
  },
  vspace: {
    flex: '1 1 0'
  },
  expand: {
    transform: 'rotate(0deg)',
    marginLeft: 'auto',
    transition: theme.transitions.create('transform', {
      duration: theme.transitions.duration.shortest
    })
  },
  limiter: {
    position: 'absolute',
    bottom: 0,
    left: 0,
    right: 0,
    height: '1rem',
    background: 'linear-gradient(180deg, rgba(255,255,255,0) 0%, rgba(255,255,255,1) 100%)'
  },
  expandOpen: {
    transform: 'rotate(180deg)'
  }
}))

/**
 * Shows an informative overview about the selected entry.
 */
export default function CollapsibleCard({height, title, action, content, fixedContent}) {
  const classes = useStyles()
  const [expanded, setExpanded] = React.useState(false)
  const [expandable, setExpandable] = React.useState(false)
  const [expandedHeight, setExpandedHeight] = React.useState(false)
  const root = useRef(null)
  const collapsible = useRef(null)
  const handleExpandClick = () => {
    setExpanded(!expanded)
  }

  // Figure out the expanded size and whether the card is expandable once on
  // startup
  useEffect(() => {
    setExpandedHeight(root.current.offsetHeight + collapsible.current.scrollHeight - collapsible.current.offsetHeight)
    const hasOverflowingChildren = collapsible.current.offsetHeight < collapsible.current.scrollHeight ||
                                    collapsible.current.offsetWidth < collapsible.current.scrollWidth
    setExpandable(hasOverflowingChildren)
  }, [])

  return <Card ref={root} className={classes.root} style={{height: expanded ? expandedHeight : height}}>
    <CardHeader
      title={title}
      className={classes.cardHeader}
      action={action}
    />
    <div ref={collapsible} className={classes.cardContent}>
      <div style={{boxSizing: 'border-box'}}>
        {content}
        { expandable && !expanded
          ? <div className={classes.limiter}></div>
          : null
        }
      </div>
    </div>
    <div className={classes.vspace}></div>
    {fixedContent
      ? <div className={classes.cardFixedContent}>
        {fixedContent}
      </div>
      : null
    }
    <div className={classes.vspace}></div>
    <CardActions
      disableSpacing
      className={classes.cardActions}
    >
      {expandable
        ? <IconButton
          size='small'
          className={clsx(classes.expand, {
            [classes.expandOpen]: expanded
          })}
          disabled={!expandable}
          onClick={handleExpandClick}
          aria-expanded={expanded}
          aria-label="show more"
        >
          <ExpandMoreIcon />
        </IconButton>
        : null
      }
    </CardActions>
  </Card>
}

CollapsibleCard.propTypes = {
  height: PropTypes.string,
  title: PropTypes.string,
  action: PropTypes.object,
  content: PropTypes.object,
  fixedContent: PropTypes.object
}
CollapsibleCard.defaultProps = {
  height: '32rem'
}
