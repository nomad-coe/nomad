module.exports = {
  InjectManifest: options => {
    // We override the default workbox cache max size.
    options.maximumFileSizeToCacheInBytes = 20 * 1024 * 1024
    return options
  }
}
