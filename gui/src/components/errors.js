import React from "react"
import { SnackbarContent, IconButton, Snackbar, withStyles } from "@material-ui/core";
import ErrorIcon from '@material-ui/icons/Error'
import CloseIcon from '@material-ui/icons/Close'

export const ErrorContext = React.createContext({
  errors: [],
  raiseError: (err) => { console.error('Wrong usage of ErrorContext.')},
});

class ErrorSnacksUnstyled extends React.Component {
  static styles = theme => ({
    root: {},
    errorSnack: {
      backgroundColor: theme.palette.error.dark,
    },
    icon: {
      fontSize: 20,
      opacity: 0.9,
      marginRight: theme.spacing.unit,
    },
    message: {
      display: 'flex',
      alignItems: 'center',
    },
    errors: {
      marginLeft: theme.spacing.unit
    }
  })
  state = {
    errors: [],
    raiseError: (error) => {
      this.onError(error)
    }
  }

  onError(error) {
    this.setState({errors: [error, ...this.state.errors]})
  }

  onClose() {
    this.setState({errors: []})
  }

  render() {
    const { children, classes } = this.props
    return (
      <ErrorContext.Provider value={this.state}>
        {children}
        <Snackbar
          anchorOrigin={{vertical: 'bottom', horizontal: 'left'}}
          open={this.state.errors.length > 0}
          // autoHideDuration={6000}
          onClose={this.handleClose}
        >
          <SnackbarContent
            className={classes.errorSnack}
            message={
              <span id="client-snackbar" className={classes.message}>
                <ErrorIcon className={classes.icon} />
                <div className={classes.errors}>
                  {this.state.errors.map((error, index) => (
                    <p key={index}>{'' + error}</p>
                  ))}
                </div>
              </span>
            }
            action={[
              <IconButton color="inherit" onClick={this.onClose.bind(this)}>
                <CloseIcon />
              </IconButton>,
            ]}
          />
        </Snackbar>
      </ErrorContext.Provider>
    )
  }
}

export const ErrorSnacks = withStyles(ErrorSnacksUnstyled.styles)(ErrorSnacksUnstyled)

export function withErrors(Component) {
  return function WithErrorComponent(props) {
    return (
      <ErrorContext.Consumer>
        {errorContext => <Component {...props} raiseError={errorContext.raiseError} />}
      </ErrorContext.Consumer>
    )
  }
}