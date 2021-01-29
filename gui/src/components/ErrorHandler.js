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
import clsx from 'clsx'
import { makeStyles } from '@material-ui/core/styles'
import PropTypes from 'prop-types'
import {
  Box,
  Button,
  Card,
  CardContent,
  Typography
} from '@material-ui/core'
import {
  Error
} from '@material-ui/icons'

export class ErrorHandler extends React.Component {
  state = {
    hasError: false
  }

  static getDerivedStateFromError(error) {
    return {
      hasError: true,
      error: error
    }
  }

  componentDidCatch(error, errorInfo) {
    console.log(error, errorInfo)
  }

  render() {
    if (this.state.hasError) {
      let msg = typeof this.props.message === 'string' ? this.props.message : this.props.message(this.state.error)
      return <ErrorCard
        message={msg}
        className={this.props.className}
        classes={this.props.classes}
      ></ErrorCard>
    }
    return this.props.children
  }
}
ErrorHandler.propTypes = ({
  children: PropTypes.object,
  message: PropTypes.oneOf([PropTypes.string, PropTypes.func]), // Provide either a fixed error message or a callback that will receive the error details.
  classes: PropTypes.object,
  className: PropTypes.string
})

export function ErrorCard({message, className, classes, actions}) {
  const useStyles = makeStyles((theme) => {
    return {
      root: {
        color: theme.palette.error.main
      },
      content: {
        paddingBottom: '16px'
      },
      'content:last-child': {
        paddingBottom: '16px !important'
      },
      title: {
        marginBottom: 0
      },
      pos: {
        marginBottom: 12
      },
      row: {
        display: 'flex'
      },
      actions: {
        display: 'flex',
        justifyContent: 'flex-end'
      },
      column: {
        display: 'flex',
        flexDirection: 'column'
      },
      errorIcon: {
        marginRight: theme.spacing(1)
      }
    }
  })

  const style = useStyles(classes)

  return <Card className={clsx(style.root, className)}>
    <CardContent className={[style.content, style['content:last-child']].join(' ')}>
      <Box className={style.row}>
        <Error className={style.errorIcon}/>
        <Box className={style.column}>
          <Typography className={style.title} color="error" gutterBottom>
            {message}
          </Typography>
          {actions
            ? <Box className={style.actions}>
              {actions.map((action) => <Button key={action.label} onClick={action.onClick}>
                {action.label}
              </Button>
              )}
            </Box>
            : ''
          }
        </Box>
      </Box>
    </CardContent>
  </Card>
}

ErrorCard.propTypes = ({
  message: PropTypes.string,
  classes: PropTypes.object,
  className: PropTypes.string,
  actions: PropTypes.array
})

export const withErrorHandler = (WrappedComponent, message) => props => (
  <ErrorHandler message={message}>
    <WrappedComponent {...props}></WrappedComponent>
  </ErrorHandler>)
