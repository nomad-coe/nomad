
export default class Adaptor {
  constructor(e) {
    this.e = e

    if (new.target === Adaptor) {
      throw new TypeError('Cannot construct Abstract instances directly')
    }
  }

  itemAdaptor(key) {
    return null
  }

  render() {
    return ''
  }
}
