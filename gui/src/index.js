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
import { sendTrackingData, matomoUrl, matomoSiteId } from './config'

const matomo = sendTrackingData ? PiwikReactRouter({
  url: matomoUrl,
  siteId: matomoSiteId
}) : null

ReactDOM.render(
  <Router history={sendTrackingData ? matomo.connectToHistory(history) : history}>
    <App />
  </Router>, document.getElementById('root'))
registerServiceWorker()
