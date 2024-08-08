const CracoWorkboxPlugin = require('craco-workbox')

module.exports = {
  plugins: [{
    plugin: CracoWorkboxPlugin
  }],
  webpack: {
    configure: (webpackConfig, { env, paths }) => {
      // Add a rule to handle .mjs files
      webpackConfig.module.rules.push({
        test: /\.mjs$/,
        include: /node_modules/,
        type: 'javascript/auto'
      })

      return webpackConfig
    }
  }
}
