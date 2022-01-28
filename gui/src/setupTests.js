/* eslint-disable */
/**
 * Code to configure or set up the testing framework before each test file in
 * the suite is executed. Contains e.g. global setup/teardown functionality for
 * tests.
 */
global.nomadEnv = {
  'keycloakBase': 'https://nomad-lab.eu/fairdi/keycloak/auth/',
  'keycloakRealm': 'fairdi_nomad_test',
  'keycloakClientId': 'nomad_gui_dev',
  'appBase': 'http://nomad-lab.eu/prod/rae/beta',
  'debug': false,
  'version': {
    'label': '1.1.0',
    'isBeta': false,
    'isTest': true,
    'usesBetaData': true,
    'officialUrl': 'https://nomad-lab.eu/prod/rae/gui'
  },
  'encyclopediaBase': 'https://nomad-lab.eu/prod/rae/encyclopedia/#',
  'aitoolkitEnabled': true,
  'oasis': false
}
// Increased the default jest timeout for individual tests
// eslint-disable-next-line no-undef
jest.setTimeout(10000)

const { ResizeObserver } = window

beforeAll(() => {
  // ResizeObserver mock init
  delete window.ResizeObserver
  window.ResizeObserver = jest.fn().mockImplementation(() => ({
    observe: jest.fn(),
    unobserve: jest.fn(),
    disconnect: jest.fn()
  }))
})

afterAll(() => {
  // ResizeObserver mock cleanup
  window.ResizeObserver = ResizeObserver
  jest.restoreAllMocks()
})

