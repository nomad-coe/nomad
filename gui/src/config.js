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
import { createTheme } from '@material-ui/core'

window.nomadEnv = window.nomadEnv || {}
export const version = window.nomadEnv.version
export const appBase = window.nomadEnv.appBase.replace(/\/$/, '')
// export const apiBase = 'http://nomad-lab.eu/prod/rae/api'
export const apiBase = `${appBase}/api`
export const northBase = window.nomadEnv.northBase
export const guiBase = process.env.PUBLIC_URL
export const servicesUploadLimit = window.nomadEnv.servicesUploadLimit
export const keycloakBase = window.nomadEnv.keycloakBase
export const keycloakRealm = window.nomadEnv.keycloakRealm
export const keycloakClientId = window.nomadEnv.keycloakClientId
export const debug = window.nomadEnv.debug || false
export const encyclopediaBase = window.nomadEnv.encyclopediaBase
export const aitoolkitEnabled = window.nomadEnv.aitoolkitEnabled || false
export const oasis = window.nomadEnv.oasis || false
export const globalLoginRequired = window.nomadEnv.globalLoginRequired || false
export const email = 'support@nomad-lab.eu'
export const maxLogsToShow = 50

export const nomadPrimaryColor = {
  main: '#008DC3',
  light: '#03B9FF',
  dark: '#005271',
  veryLight: '#99e2ff'
}

export const nomadSecondaryColor = {
  main: '#00CED1',
  light: '#54DCDC',
  dark: '#007C7C',
  veryLight: '#B5F0F0'
}

export const nomadFontFamily = [
  'Titillium Web',
  'sans-serif'
].join(',')

export const nomadTheme = createTheme({
  typography: {
    useNextVariants: true,
    fontFamily: nomadFontFamily,
    fontWeightRegular: 400,
    fontWeightMedium: 600
  },
  palette: {
    primary: nomadPrimaryColor,
    secondary: nomadSecondaryColor
  },
  overrides: {
    // This is used to inject global css styles through the CssBaseline
    // component, see: https://v4.mui.com/customization/globals/#global-css
    MuiCssBaseline: {
      '@global': {
        '.react-grid-item.react-grid-placeholder': {
            backgroundColor: nomadSecondaryColor.main
        },
        // This rule is overwritten. The normalization was added quite late into
        // the app and the 'box-sizing: inherit' -rule would mess up the layout
        // quite a bit.
        '*, *::before, *::after': {
          boxSizing: 'revert'
        }
      }
    },
    MuiTooltip: {
      tooltip: {
        fontWeight: 'normal',
        fontSize: '0.75rem'
      }
    },
    MuiTableRow: {
      root: {
        '&:last-child td': {
          borderBottom: 0
        }
      }
    }
  }
})

export const formatQuantity = (x) => {
  const parts = x.toString().split('.')
  parts[0] = parts[0].replace(/\B(?=(\d{3})+(?!\d))/g, ',')
  return parts.join('.')
}

export function normalizeDisplayValue(value) {
  if (value === 'not processed' || value === 'unavailbale') {
    return '-'
  }
  return value
}

/**
 * The available aggregation sizes for an ES aggregation result.
 */
export const aggregationSizes = [10, 20, 30, 40, 50, 100, 200]
/**
 * The range of electronic energies in eV.
 */
export const electronicRange = [-5, 10]
/**
 * The default date format.
 */
export const dateFormat = 'dd/MM/yyyy'
/**
 * Warning shown when energy values could not be normalized.
 */
export const msgNormalizationWarning = `
Energy reference could not be found: energies have an unknown shift with
respect to the highest occupied energy.
`
