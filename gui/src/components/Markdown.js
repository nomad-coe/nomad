import React from 'react'
import PropTypes from 'prop-types'
import marked from 'marked'
import { withStyles } from '@material-ui/core'
import extend from '@babel/runtime/helpers/extends'

/**
 * A simple markdown component.
 *
 * It uses marked with styled 'borrowed' from the materials ui docs system.
 */

var styles = theme => ({
  root: {
    fontFamily: theme.typography.fontFamily,
    fontSize: 16,
    color: theme.palette.text.primary,
    '& .anchor-link': {
      marginTop: -96,
      // Offset for the anchor.
      position: 'absolute'
    },
    '& pre, & pre[class*="language-"]': {
      margin: '24px 0',
      padding: '12px 18px',
      backgroundColor: theme.palette.background.paper,
      borderRadius: theme.shape.borderRadius,
      overflow: 'auto',
      WebkitOverflowScrolling: 'touch' // iOS momentum scrolling.

    },
    '& code': {
      display: 'inline-block',
      lineHeight: 1.6,
      fontFamily: 'Consolas, "Liberation Mono", Menlo, Courier, monospace',
      padding: '3px 6px',
      color: theme.palette.text.primary,
      backgroundColor: theme.palette.background.paper,
      fontSize: 14
    },
    '& p code, & ul code, & pre code': {
      fontSize: 14,
      lineHeight: 1.6
    },
    '& h1': (0, extend)({}, theme.typography.display2, {
      color: theme.palette.text.secondary,
      margin: '32px 0 16px'
    }),
    '& .description': (0, extend)({}, theme.typography.headline, {
      margin: '0 0 40px'
    }),
    '& h2': (0, extend)({}, theme.typography.display1, {
      color: theme.palette.text.secondary,
      margin: '32px 0 24px'
    }),
    '& h3': (0, extend)({}, theme.typography.headline, {
      color: theme.palette.text.secondary,
      margin: '32px 0 24px'
    }),
    '& h4': (0, extend)({}, theme.typography.title, {
      color: theme.palette.text.secondary,
      margin: '24px 0 16px'
    }),
    '& p, & ul, & ol': {
      lineHeight: 1.6
    },
    '& h1, & h2, & h3, & h4': {
      '& code': {
        fontSize: 'inherit',
        lineHeight: 'inherit',
        // Remove scroll on small screens.
        wordBreak: 'break-word'
      },
      '& .anchor-link-style': {
        opacity: 0,
        // To prevent the link to get the focus.
        display: 'none'
      },
      '&:hover .anchor-link-style': {
        display: 'inline-block',
        opacity: 1,
        padding: '0 8px',
        color: theme.palette.text.hint,
        '&:hover': {
          color: theme.palette.text.secondary
        },
        '& svg': {
          width: '0.55em',
          height: '0.55em',
          fill: 'currentColor'
        }
      }
    },
    '& table': {
      width: '100%',
      display: 'block',
      overflowX: 'auto',
      WebkitOverflowScrolling: 'touch',
      // iOS momentum scrolling.
      borderCollapse: 'collapse',
      borderSpacing: 0,
      overflow: 'hidden',
      '& .prop-name': {
        fontSize: 13,
        fontFamily: 'Consolas, "Liberation Mono", Menlo, monospace'
      },
      '& .required': {
        color: theme.palette.type === 'light' ? '#006500' : '#9bc89b'
      },
      '& .prop-type': {
        fontSize: 13,
        fontFamily: 'Consolas, "Liberation Mono", Menlo, monospace',
        color: theme.palette.type === 'light' ? '#932981' : '#dbb0d0'
      },
      '& .prop-default': {
        fontSize: 13,
        fontFamily: 'Consolas, "Liberation Mono", Menlo, monospace',
        borderBottom: '1px dotted '.concat(theme.palette.text.hint)
      }
    },
    '& thead': {
      fontSize: 14,
      fontWeight: theme.typography.fontWeightMedium,
      color: theme.palette.text.secondary
    },
    '& tbody': {
      fontSize: 14,
      lineHeight: 1.5,
      color: theme.palette.text.primary
    },
    '& td': {
      borderBottom: '1px solid '.concat(theme.palette.divider),
      padding: '8px 16px 8px 8px',
      textAlign: 'left'
    },
    '& td:last-child': {
      paddingRight: 24
    },
    '& td compact': {
      paddingRight: 24
    },
    '& td code': {
      fontSize: 13,
      lineHeight: 1.6
    },
    '& th': {
      whiteSpace: 'pre',
      borderBottom: '1px solid '.concat(theme.palette.divider),
      fontWeight: theme.typography.fontWeightMedium,
      padding: '0 16px 0 8px',
      textAlign: 'left'
    },
    '& th:last-child': {
      paddingRight: 24
    },
    '& tr': {
      height: 48
    },
    '& thead tr': {
      height: 64
    },
    '& strong': {
      fontWeight: theme.typography.fontWeightMedium
    },
    '& blockquote': {
      borderLeft: '5px solid '.concat(theme.palette.text.hint),
      backgroundColor: theme.palette.background.paper,
      padding: '4px 24px',
      margin: '24px 0'
    },
    '& a, & a code': {
      // Style taken from the Link component
      color: theme.palette.secondary.main,
      textDecoration: 'none',
      '&:hover': {
        textDecoration: 'underline'
      }
    },
    '& img': {
      maxWidth: '100%'
    }
  }
})

function Markdown(props) {
  const { classes, text, children } = props

  let content = text
  if (children) {
    content = children.replace(/^\s+/gm, '')
  }

  return (
    <div>
      <div
        className={classes.root}
        dangerouslySetInnerHTML={{__html: marked(content)}}
      />
    </div>
  )
}

Markdown.propTypes = {
  classes: PropTypes.object.isRequired,
  text: PropTypes.string,
  children: PropTypes.string
}

export default withStyles(styles)(Markdown)
