window.nomadEnv = {
  'keycloakBase': 'https://nomad-lab.eu/fairdi/keycloak/auth/',
  // Use the production API
  // 'keycloakRealm': 'fairdi_nomad_prod',
  // 'keycloakClientId': 'nomad_public',
  // 'appBase': 'https://nomad-lab.eu/prod/v1',
  // Use the local API
  'keycloakRealm': 'fairdi_nomad_test',
  'keycloakClientId': 'nomad_gui_dev',
  'appBase': 'http://localhost:8000/fairdi/nomad/latest',
  'encyclopediaBase': 'https://nomad-lab.eu/prod/rae/encyclopedia/#',
  'northBase': 'http://localhost:9000/fairdi/nomad/latest/north',
  'debug': false,
  'version': {
    'label': '1.1.3',
    'isBeta': false,
    'isTest': true,
    'usesBetaData': true,
    'officialUrl': 'https://nomad-lab.eu/prod/rae/gui'
  },
  'aitoolkitEnabled': false,
  'oasis': true,
  'servicesUploadLimit': 10
}
