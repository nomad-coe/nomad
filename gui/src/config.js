import { createMuiTheme } from '@material-ui/core'

window.nomadEnv = window.nomadEnv || {}
export const version = window.nomadEnv.version
export const appBase = window.nomadEnv.appBase.replace(/\/$/, '')
// export const apiBase = 'http://nomad-lab.eu/prod/rae/api'
export const apiBase = `${appBase}/api`
export const optimadeBase = `${appBase}/optimade`
export const guiBase = process.env.PUBLIC_URL
export const matomoUrl = window.nomadEnv.matomoUrl
export const matomoSiteId = window.nomadEnv.matomoSiteId
export const keycloakBase = window.nomadEnv.keycloakBase
export const keycloakRealm = window.nomadEnv.keycloakRealm
export const keycloakClientId = window.nomadEnv.keycloakClientId
export const debug = window.nomadEnv.debug || false
export const matomoEnabled = window.nomadEnv.matomoEnabled || false
export const encyclopediaEnabled = window.nomadEnv.encyclopediaEnabled || false
export const aitoolkitEnabled = window.nomadEnv.aitoolkitEnabled || false
export const oasis = window.nomadEnv.oasis || false
export const email = 'support@nomad-lab.eu'
export const maxLogsToShow = 50

export const consent = `
By using this web-site and by uploading and downloading data, you agree to the
[terms of use](https://nomad-lab.eu/terms).

Uploaded data is licensed under the Creative Commons Attribution license
([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)). You can publish
uploaded data with an *embargo*. Data with an *embargo* is only visible to
you and users you share your data with. The *embargo period* lasts up to 36 months.
After the *embargo* your published data will be public. **Note that public data
is visible to others and files become downloadable by everyone.**

This web-site uses *cookies*. We use cookies to track you login status for all NOMAD services
and optionally to store information about your use of NOMAD. None of this information is
shared with other parties. By using this web-site you agree to the described use of *cookies*.
[Learn more](https://www.cookiesandyou.com/).
`
export const nomadPrimaryColor = {
  main: '#008DC3',
  light: '#03B9FF',
  dark: '#005271',
  veryLight: '#10BAFB'
}

export const nomadSecondaryColor = {
  main: '#00CED1',
  light: '#54DCDC',
  veryLight: '#B5F0F0',
  dark: '#007C7C'
}

export const nomadFontFamily = [
  'Titillium Web',
  'sans-serif'
].join(',')

export const nomadTheme = createMuiTheme({
  typography: {
    useNextVariants: true,
    fontFamily: nomadFontFamily,
    fontWeightRegular: 400,
    fontWeightMedium: 600
  },
  palette: {
    primary: nomadPrimaryColor,
    secondary: nomadSecondaryColor
  }
})

export const formatQuantity = (x) => {
  const parts = x.toString().split('.')
  parts[0] = parts[0].replace(/\B(?=(\d{3})+(?!\d))/g, ',')
  return parts.join('.')
}
