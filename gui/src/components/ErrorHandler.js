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
      let msg = this.props.errorHandler ? this.props.errorHandler(this.state.error) : this.props.message
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
  message: PropTypes.string, // Fixed error message. Provide either this or errorHandler
  errorHandler: PropTypes.func, // Function that is called once an error is caught. It recveives the error object as argument and should return an error message as string.
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
  console.log(actions)

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
