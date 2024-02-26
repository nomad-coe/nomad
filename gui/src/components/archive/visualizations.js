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
import React, { useRef, useLayoutEffect, useState } from 'react'
import PropTypes from 'prop-types'
import { FixedSizeGrid as Grid } from 'react-window'
import { Typography, makeStyles, Button, Grid as MuiGrid, Box } from '@material-ui/core'
import AutoSizer from 'react-virtualized-auto-sizer'

export function Number({value, exp, variant, unit, ...props}) {
  variant = variant || 'body2'
  exp = exp || 2
  const fixed = 5
  let html = '-'
  if (value !== undefined && value !== null) {
    if (typeof value === 'number') {
      const str = value.toExponential(exp)
      const [f, e] = str.split('e')
      if (e <= fixed && e >= -fixed) {
        html = value.toFixed(fixed).replace(/[.,]0*$/, '')
      } else {
        html = <span>{f.toString().replace(/0*$/, '')}&middot;10<sup>{e}</sup></span>
      }
    } else {
      html = value.toString()
    }
  }
  return <Typography {...props} variant={variant} data-testid="scientific-number-with-unit">{html}{unit && <span>&nbsp;{unit}</span>}</Typography>
}
Number.propTypes = ({
  value: PropTypes.any,
  variant: PropTypes.string,
  exp: PropTypes.number,
  unit: PropTypes.string
})

function MatrixPagination({length, page, onChange}) {
  return <MuiGrid
    spacing={2}
    container
    direction="row"
    justifyContent="center"
    alignItems="center"
  >
    <MuiGrid item>
      <Button disabled={page === 0} size="small" onClick={() => onChange(page - 1)}>
        prev
      </Button>
    </MuiGrid>
    <MuiGrid item>
      <Typography>{page}</Typography>
    </MuiGrid>
    <MuiGrid item>
      <Button size="small" disabled={page + 1 === length} onClick={() => onChange(page + 1)}>
      next
      </Button>
    </MuiGrid>
  </MuiGrid>
}
MatrixPagination.propTypes = ({
  length: PropTypes.number.isRequired,
  page: PropTypes.number.isRequired,
  onChange: PropTypes.func.isRequired
})

const useMatrixStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
    flexDirection: 'row',
    alignItems: 'stretch',
    minWidth: 300,
    justifyContent: 'center'
  },
  leftBracket: {
    width: theme.spacing(1),
    border: 'solid 2px black',
    borderRight: 'none',
    marginRight: -theme.spacing(1)
  },
  matrix: {
    width: '100%'
  },
  rightBracket: {
    width: theme.spacing(1),
    border: 'solid 2px black',
    borderLeft: 'none',
    marginLeft: -theme.spacing(1)
  }
}))
export function Matrix({values, shape, invert, type}) {
  const rootRef = useRef()
  const matrixRef = useRef()
  const [pages, setPages] = useState(new Array(Math.max(0, shape.length - 2)).fill(0))
  const pageLengths = []
  const classes = useMatrixStyles()

  let ii = 0
  for (let i = shape.length; i > 2; i--) {
    pageLengths.push(values.length)
    values = values[pages[ii++]]
  }

  const columnWidth = useRef(92)
  const rowHeight = 24
  const rowCount = invert ? values.length : shape.length > 1 ? values[0].length : 1
  const columnCount = invert ? shape.length > 1 ? values[0].length : 1 : values.length
  const height = Math.min(300, (rowCount - 1) * rowHeight + 24)

  useLayoutEffect(() => {
    if (type === 'str') {
      matrixRef.current.style.width = '100%'
    } else {
      matrixRef.current.style.width = `${rootRef.current.clientWidth - 4}px`
    }
    columnWidth.current = Math.max(92, (rootRef.current.clientWidth - 4) / columnCount)
  })

  let value = shape.length > 1 ? ({rowIndex, columnIndex}) => values[columnIndex][rowIndex] : ({columnIndex}) => values[columnIndex]
  if (invert) {
    value = shape.length > 1 ? ({rowIndex, columnIndex}) => values[rowIndex][columnIndex] : ({rowIndex}) => values[rowIndex]
  }

  return <React.Fragment>
    <div ref={rootRef} className={classes.root}>
      <div className={classes.leftBracket} style={{height: height}}>&nbsp;</div>
      <div ref={matrixRef} className={classes.matrix}>
        <AutoSizer>
          {({width}) => (
            <Grid
              columnCount={columnCount}
              columnWidth={columnWidth.current}
              height={height}
              rowCount={rowCount}
              rowHeight={rowHeight}
              width={width}
            >
              {({style, ...props}) => <Number style={{whiteSpace: 'nowrap', ...style}} value={value(props)} />}
            </Grid>
          )}
        </AutoSizer>
      </div>
      <div className={classes.rightBracket} style={{height: height}}>&nbsp;</div>
    </div>
    {pages.map((page, index) => <Box margin={1} key={index}>
      <MatrixPagination
        length={pageLengths[index]} page={page}
        onChange={(page) => setPages([...pages.slice(0, index), page, ...pages.slice(index + 1)])}
      />
    </Box>)}
  </React.Fragment>
}
Matrix.propTypes = ({
  values: PropTypes.array.isRequired,
  shape: PropTypes.arrayOf(PropTypes.any).isRequired,
  invert: PropTypes.bool,
  type: PropTypes.string
})
