const express = require('express');
const session = require('express-session');
const passport = require('passport');
const OAuth2Strategy = require('passport-oauth2');
const axios = require('axios')

const port = 8888
const app = express();
const baserouter = express.Router();

const clientHubApiUrl = 'http://localhost:9000/fairdi/nomad/latest/north/hub/api'; // process.env.JUPYTERHUB_API_URL
const serverHubApiUrl = 'http://host.docker.internal:9000/fairdi/nomad/latest/north/hub/api';
const baseurl = process.env.JUPYTERHUB_SERVICE_PREFIX || process.env.SUBFOLDER || '/';
const secret = process.env.JUPYTERHUB_API_TOKEN
const user = process.env.JUPYTERHUB_USER

const passportOptions = {
  authorizationURL: `${clientHubApiUrl}/oauth2/authorize`,
  tokenURL: `${serverHubApiUrl}/oauth2/token`,
  clientID: process.env.JUPYTERHUB_CLIENT_ID,
  clientSecret: secret
}

passport.use(new OAuth2Strategy(
  passportOptions,
  function(accessToken, refreshToken, params, profile, done) {
    axios.get(`${serverHubApiUrl}/user`, {
      headers: { 'Authorization': `Bearer ${params['access_token']}`}
    }).then(response => {
      if (!response?.data?.name) {
        done('Cannot info for loggedin user to authorize access.', null);
      } else if (response?.data?.name !== user) {
        done('Logged in user does not match the container\'s user', null);
      } else {
        done(null, response.data);
      }
    }).catch(error => done(error, null))
  }
));

passport.serializeUser(function(user, done) {
  done(null, user);
});

passport.deserializeUser(function(user, done) {
  done(null, user);
});

baserouter.use(session({ secret: secret, cookie: { maxAge: 60000, path: baseurl }}))
baserouter.use(passport.initialize());
baserouter.use(passport.session());

baserouter.get('/', (req, res, next) => {
  if (!req.user) {
    res.redirect(`${baseurl}/login`);
    return
  }
  res.send('Welcome Home');
});

baserouter.get('/login', passport.authenticate('oauth2'));

baserouter.get('/oauth_callback',
  passport.authenticate('oauth2'),
  function(req, res) {
    res.redirect(baseurl);
  });

app.use((req, res, next) => {
  if (!req.path.startsWith(baseurl)) {
    res.redirect(301, baseurl)
  } else {
    next()
  }
})

app.use(baseurl, baserouter);

app.listen(port, function () {
  console.log('Example app listening on port ' + port + '!');
});
