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
import React from 'react'
import PropTypes from 'prop-types'
// import remark from 'remark'
import { withStyles, Link } from '@material-ui/core'
import { Link as RouterLink } from 'react-router-dom'
import extend from '@babel/runtime/helpers/extends'
import ReactMarkdown from 'react-markdown'
import MathJax from 'react-mathjax'
import RemarkMathPlugin from 'remark-math'
import { useGlobalMetainfo } from './archive/metainfo'
import { appBase } from '../config'

/**
 * A simple markdown component.
 *
 * It uses marked with styled 'borrowed' from the materials ui docs system.
 */

const styles = theme => {
  return {
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
        backgroundColor: theme.palette.primary.veryLight,
        borderRadius: theme.shape.borderRadius,
        overflow: 'auto',
        WebkitOverflowScrolling: 'touch' // iOS momentum scrolling.

      },
      '& code': {
        display: 'inline-block',
        lineHeight: 1,
        fontFamily: 'Consolas, "Liberation Mono", Menlo, Courier, monospace',
        padding: '5px 4px 2px 4px',
        color: theme.palette.text.primary.dark,
        backgroundColor: theme.palette.primary.veryLight,
        borderRadius: theme.shape.borderRadius,
        fontSize: 14
      },
      '& p code, & ul code': {
        fontSize: 14
      },
      '& pre code': {
        fontSize: 14,
        lineHeight: 1.6
      },
      '& p:first-child': {
        marginTop: 0
      },
      '& p:last-child': {
        marginBottom: 0
      },
      '& h1': (0, extend)({}, theme.typography.h3, {
        color: theme.palette.text.primary,
        margin: '32px 0 16px'
      }),
      '& .description': (0, extend)({}, theme.typography.h5, {
        margin: '0 0 40px'
      }),
      '& h2': (0, extend)({}, theme.typography.h4, {
        color: theme.palette.text.primary,
        margin: '32px 0 24px'
      }),
      '& h3': (0, extend)({}, theme.typography.h5, {
        color: theme.palette.text.primary,
        margin: '32px 0 24px'
      }),
      '& h4': (0, extend)({}, theme.typography.h6, {
        color: theme.palette.text.primary,
        margin: '24px 0 16px'
      }),
      '& p, & ul, & ol': {
        lineHeight: theme.typography.lineHeight,
        marginBottom: theme.spacing(1)
      },
      '& ul': {
        paddingLeft: 0,
        '& li': {
          listStyleType: 'none',
          fontSize: 'inherit',
          paddingLeft: theme.spacing(4),
          '&:before': {
            content: '\'■\'',
            fontSize: 'x-large',
            marginLeft: -theme.spacing(1) * 4,
            paddingRight: theme.spacing(4) - 14
          }
        }
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
        color: theme.palette.primary.main,
        textDecoration: 'none',
        '&:hover': {
          textDecoration: 'underline'
        }
      },
      '& img': {
        maxWidth: '100%'
      }
    }
  }
}

function MarkdownLink({href, ...props}) {
  if (href.match(/^https?:\/\//) || href.startsWith(appBase)) {
    return <Link href={href} {...props} />
  } else {
    return <Link component={RouterLink} to={href} {...props} />
  }
}
MarkdownLink.propTypes = {
  href: PropTypes.string.isRequired
}

function Markdown(props) {
  const { classes, text, children, ...moreProps } = props
  const globalMetainfo = useGlobalMetainfo()

  let content = text
  if (children) {
    const state = []
    // We modify the children for the following reasons:
    // - remove unnecessary whitespaces that might break the markdown layout
    // - escape _ (usually in metainfo names) with sensitivity to _ used in latex math
    // - put metainfo names into code quotes
    // - turn metainfo names into links into the metainfo
    let word = ''
    const letters = children.replace(/^ +/gm, '').split('')
    content = letters.map((c, i) => {
      if (c === '`' || c === '$') {
        if (state.peek === c) {
          state.pop()
        } else {
          state.push(c)
        }
      }

      if (c === '[') {
        state.push(c)
      }
      if (c === ']' && state.peek === '[') {
        state.pop()
      }

      if (c.match(/[a-zA-Z0-9_]/)) {
        word += c
        if (letters.length - 1 === i) {
          return word
        }
      } else {
        if (state.length === 0 && (word.match(/_/g) || word.match(/[a-z]+[A-Z]/g)) && word.match(/^[a-zA-Z0-9_]+$/g) && c !== ']') {
          const path = globalMetainfo?.path(word)
          if (path) {
            word = `[\`${word}\`](/metainfo/${path})`
          } else if (word.includes('_')) {
            word = `\`${word}\``
          }
        }
        const result = word + c
        word = ''
        return result
      }
      return ''
    }).join('')
  }

  const math = ({value}) => <MathJax.Node formula={value} />
  const inlineMath = ({value}) => <MathJax.Node inline formula={value} />
  const newProps = {
    ...moreProps,
    children: content,
    remarkPlugins: [
      RemarkMathPlugin
    ],
    components: {
      ...moreProps.components,
      math: math,
      inlineMath: inlineMath,
      a: props => <MarkdownLink {...props} />
    }
  }
  const md = (
    <MathJax.Provider input="tex">
      <ReactMarkdown {...newProps} />
    </MathJax.Provider>
  )

  return (
    <div className={classes.root}>{md}</div>
  )
}

Markdown.propTypes = {
  classes: PropTypes.object.isRequired,
  text: PropTypes.string,
  children: PropTypes.string
}

export default withStyles(styles)(Markdown)
