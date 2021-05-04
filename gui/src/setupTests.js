global.nomadEnv = {
  'keycloakBase': 'https://nomad-lab.eu/fairdi/keycloak/auth/',
  'keycloakRealm': 'fairdi_nomad_test',
  'keycloakClientId': 'nomad_gui_dev',
  'appBase': 'http://nomad-lab.eu/prod/rae/beta',
  'debug': false,
  'matomoEnabled': false,
  'matomoUrl': 'https://nomad-lab.eu/fairdi/stat',
  'matomoSiteId': '2',
  'version': {
    'label': '0.10.3',
    'isBeta': false,
    'isTest': true,
    'usesBetaData': true,
    'officialUrl': 'https://nomad-lab.eu/prod/rae/gui'
  },
  'encyclopediaEnabled': true,
  'aitoolkitEnabled': true,
  'oasis': false
}
// Increased the default jest timeout for individual tests
// eslint-disable-next-line no-undef
jest.setTimeout(10000)
