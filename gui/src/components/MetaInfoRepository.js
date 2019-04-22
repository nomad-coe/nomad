class BadNomadMIError extends Error {

}

class Schema {
  category = 'category'
  section = 'section'
  property = 'property'
  value = 'value'
  reference = 'reference'
  pkg = 'package'

  isCategory = (element) => element.mType === this.category
  isFeature = (element) => element.mType === this.section || this.isProperty(element)
  isSection = (element) => element.mType === this.section
  isReference = (element) => element.mType === this.reference
  isValue = (element) => element.mType === this.value
  isProperty = (element) => this.isValue(element) || this.isReference(element)
  isDefinition = (element) => this.isCategory(element) || this.isFeature(element)
  isPakage = (element) => element.mType === this.pkg
  isElement = (element) => this.isPakage(element) || this.isDefinition(element)

  type = (type) => {
    if (type === 'C') {
      return 'string'
    } else if (type === 'f') {
      return 'float'
    } else if (type === 'i') {
      return 'integer'
    } else {
      return type
    }
  }

  isParent = (child, parent) => {
    let current = child.parent
    while (current) {
      if (parent === current) {
        return true
      }
      current = current.parent
    }
    return false
  }

  allContents = (element, func) => {
    if (this.isPakage(element)) {
      (element.definitions || []).forEach(element => this.allContents(element, func))
    }
    func(element)
  }
}

export const schema = new Schema()

export default class MetaInfoRepository {
  constructor(metaInfoPackages) {
    this.contents = []
    this.names = {}
    this.warnings = []
    this.errors = []
    Object.keys(metaInfoPackages).forEach(packageName => this.addNomadMiJson(packageName, metaInfoPackages[packageName]))

    this.resolveAll()
    this.resolveSuperNames()
    this.addReverseRefs()
    this.addCategoryFeaturesToSections()

    this.warnings.forEach(warning => console.warning(warning))
    if (this.errors.length) {
      console.log(this.errors[0])
      console.log('there might be more errors')
    }
    // this.errors.forEach(error => console.error(error))
  }

  createProxy(reference) {
    return {
      mIsProxy: true,
      mReference: reference
    }
  }

  get(metaInfoName) {
    return this.names[metaInfoName]
  }

  resolve(proxy) {
    if (proxy.mIsProxy) {
      const resolved = this.names[proxy.mReference]
      return resolved || proxy
    } else {
      return proxy
    }
  }

  resolveSuperNames() {
    this.allContents(element => {
      if (element._superNames) {
        element._superNames.forEach(parentOrSuper => {
          if (schema.isSection(parentOrSuper)) {
            if (schema.isFeature(element)) {
              if (element.parent) {
                // this.errors.push(
                //   new BadNomadMIError(`More than one parent in feature ${element.name}`))
              } else {
                element.parent = parentOrSuper
              }
            } else {
              if (element.parent) {
                this.errors.push(
                  new BadNomadMIError(`More than one section for category ${element.name}`))
              } else {
                element.parent = parentOrSuper
              }
            }
          } else if (schema.isCategory(parentOrSuper)) {
            if (schema.isCategory(element)) {
              element.super = element.super || []
              element.super.push(parentOrSuper)
            } else if (schema.isFeature(element)) {
              element.categories = element.categories || []
              element.categories.push(parentOrSuper)
            } else {
              this.errors.push(new BadNomadMIError(`Non feature ${element.name} references category ${parentOrSuper.name}.`))
            }
          } else if (schema.isProperty(parentOrSuper)) {
            this.warnings.push(`SuperName in ${element.name} references property ${parentOrSuper.name}. That is not allowed ?!`)
          } else {
            this.errors.push(new BadNomadMIError(`Referenced parent or super ${parentOrSuper.name} is not a section or category (in ${element.name})`))
          }
        })
      }
    })
  }

  resolveAll() {
    this.allContents(element => {
      if (element._superNames) {
        element._superNames = element._superNames.map((ref) => {
          const resolved = this.resolve(ref)
          if (resolved.mIsProxy) {
            element.problems.push(`Could not resolve parent section ${ref.mReference} in ${element.name}.`)
          }
          return resolved
        })
      }
      if (element.referencedSection) {
        const ref = element.referencedSection
        const resolved = this.resolve(ref)
        if (resolved.mIsProxy) {
          element.problems.push(`Could not resolve referenced section ${ref.mReference} in ${element.name}.`)
        }
        element.referencedSection = resolved
      }
    })
  }

  addReverseRefs() {
    this.allContents(definition => {
      if (definition.parent) {
        definition.parent.features = definition.parent.features || []
        definition.parent.features.push(definition)
      }
      if (schema.isPakage(definition)) {
        definition.definitions.forEach(feature => {
          feature.package = definition
        })
      }
    })
  }

  addCategoryFeaturesToSections() {
    this.allContents(definition => {
      if (schema.isProperty(definition) && definition.categories) {
        definition.categories.forEach(category => {
          const section = category.parent
          if (section) {
            section.features = section.features || []
            if (section.features.indexOf(definition) === -1) {
              section.features.push(definition)
              definition.parent = section
            }
          }
        })
      }
    })
  }

  allContents(func) {
    this.contents.forEach(element => schema.allContents(element, func))
  }

  addName(namedElement) {
    const {name} = namedElement
    if (this.names[name]) {
      this.errors.push(new BadNomadMIError(`Element with name ${namedElement.name} does already exist.`))
    } else {
      this.names[name] = namedElement
    }
  }

  addNomadMiJson(name, json) {
    const transformMetaInfo = (metaInfo) => {
      const isMeta = metaInfo.kindStr === 'type_meta'
      const isSection = metaInfo.kindStr === 'type_section'
      const isCategory = metaInfo.kindStr === 'type_abstract_document_content'
      const isProperty = metaInfo.dtypeStr && metaInfo.dtypeStr !== 'r'
      const isReference = metaInfo.dtypeStr === 'r' || metaInfo.referencedSections
      const isValue = isProperty && !isReference

      const superNames = metaInfo.superNames || []
      const definition = {
        name: metaInfo.name,
        description: metaInfo.description,
        miJson: metaInfo,
        type: schema.type(metaInfo.dtypeStr),
        problems: [],
        _superNames: superNames.map(ref => this.createProxy(ref))
      }

      if (isSection) {
        definition.mType = schema.section
      } else if (isCategory) {
        definition.mType = schema.category
      } else if (isValue) {
        definition.mType = schema.value
      } else if (isReference) {
        definition.mType = schema.reference
        if (!metaInfo.referencedSections || metaInfo.referencedSections.length < 1) {
          definition.problems.push(`Reference ${definition.name} does not reference anything.`)
        } else {
          definition.referencedSection = this.createProxy(metaInfo.referencedSections[0])
        }
      } else if (isMeta) {
        // ignore meta-meta definitions
      } else {
        this.errors.push(new BadNomadMIError(`Cannot determine mType ${metaInfo.kindStr} of feature ${name}:${metaInfo.name}`))
      }

      this.addName(definition)
      return definition
    }

    const metaInfos = json.metaInfos || []
    const pkg = {
      mType: schema.pkg,
      name: name,
      description: json.description,
      definitions: metaInfos.map(transformMetaInfo)
    }

    this.addName(pkg)
    this.contents.push(pkg)
  }
}
