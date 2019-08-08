import 'react-app-polyfill/ie11'
import 'react-app-polyfill/stable'
import React from 'react'
import ReactDOM from 'react-dom'
import './index.css'
import App from './components/App'
import registerServiceWorker from './registerServiceWorker'
import { Router } from 'react-router-dom'
import history from './history'
import PiwikReactRouter from 'piwik-react-router'
import { debug } from './config'

const matomo = PiwikReactRouter({
  url: 'https://labdev-nomad.esc.rzg.mpg.de/fairdi/matomo/',
  siteId: 1
})

console.log(debug)

ReactDOM.render(
  <Router history={debug ? history : matomo.connectToHistory(history)}>
    <App />
  </Router>, document.getElementById('root'))
registerServiceWorker()
