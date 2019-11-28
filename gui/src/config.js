import { createMuiTheme } from '@material-ui/core'

window.nomadEnv = window.nomadEnv || {}
export const appBase = window.nomadEnv.appBase.replace(/\/$/, '')
export const apiBase = `${appBase}/api`
export const optimadeBase = `${appBase}/optimade`
export const guiBase = process.env.PUBLIC_URL
export const kibanaBase = window.nomadEnv.kibanaBase
export const matomoUrl = window.nomadEnv.matomoUrl
export const matomoSiteId = window.nomadEnv.matomoSiteId
export const keycloakBase = window.nomadEnv.keycloakBase
export const keycloakRealm = window.nomadEnv.keycloakRealm
export const keycloakClientId = window.nomadEnv.keycloakClientId
export const debug = window.nomadEnv.debug || false
export const sendTrackingData = window.nomadEnv.sendTrackingData

export const consent = `
By using this web-site and by uploading and downloading data, you agree to the
[terms of use](https://www.nomad-coe.eu/the-project/nomad-repository/nomad-repository-terms).

Uploaded data is licensed under the Creative Commons Attribution license
([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)). You can publish
uploaded data with an *embargo*. Data with an *embargo* is only visible to
you and users you share your data with. The *embargo period* lasts up to 36 month.
After the *embargo* your published data will be public. **Note that public data
is visible to others and files become downloadable by everyone.**

This web-site uses *cookies*. By using this web-site you agree to our use
of *cookies*. [Learn more](https://www.cookiesandyou.com/).
`
export const nomadPrimaryColor = {
  main: '#294277',
  veryLight: '#cfdeff'
}

export const nomadTheme = createMuiTheme({
  typography: {
    useNextVariants: true
  },
  palette: {
    primary: nomadPrimaryColor,
    secondary: nomadPrimaryColor
  }
})

export const formatQuantity = (x) => {
  const parts = x.toString().split('.')
  parts[0] = parts[0].replace(/\B(?=(\d{3})+(?!\d))/g, ',')
  return parts.join('.')
}
