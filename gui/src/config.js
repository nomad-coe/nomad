import repo from '@material-ui/core/colors/deepPurple'
import archive from '@material-ui/core/colors/teal'
import enc from '@material-ui/core/colors/amber'
import analytics from '@material-ui/core/colors/lightGreen'
import secondary from '@material-ui/core/colors/blueGrey'
import { createMuiTheme } from '@material-ui/core'

window.nomad_env = window.nomad_env || {}
export const apiBase = process.env.REACT_APP_API_BASE || window.nomad_env.apiBase
export const appBase = process.env.REACT_APP_APP_BASE || window.nomad_env.appBase
export const appStaticBase = process.env.REACT_APP_APP_STATIC_BASE || window.nomad_env.appStaticBase
export const debug = process.env.REACT_APP_DEBUG ? process.env.REACT_APP_DEBUG === 'true' : window.nomad_env.debug

export const genTheme = createMuiTheme({
  palette: {
    primary: secondary,
    secondary: secondary
  }
})

export const repoTheme = createMuiTheme({
  palette: {
    primary: repo,
    secondary: repo
  }
})

export const archiveTheme = createMuiTheme({
  palette: {
    primary: archive,
    secondary: repo
  }
})

export const encTheme = createMuiTheme({
  palette: {
    primary: enc,
    secondary: repo
  }
})

export const analyticsTheme = createMuiTheme({
  palette: {
    primary: analytics,
    secondary: repo
  }
})
