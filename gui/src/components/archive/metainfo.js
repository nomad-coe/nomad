import metainfo from '../../metainfo'

export const sectionDefs = {}
export const packageDefs = {}

metainfo.packages.forEach(pkg => {
  packageDefs[pkg.name] = pkg
  pkg._sections = {}
  if (pkg.category_definitions) {
    pkg.category_definitions.forEach(categoryDef => {
      categoryDef._qualifiedName = `${pkg.name}:${categoryDef.name}`
    })
  }
  pkg.section_definitions.forEach(sectionDef => {
    pkg._sections[sectionDef.name] = sectionDef
    sectionDefs[sectionDef.name] = sectionDef
    sectionDef.quantities = sectionDef.quantities || []
    sectionDef.sub_sections = sectionDef.sub_sections || []
    sectionDef._incomingRefs = sectionDef._incomingRefs || []
    sectionDef._parentSections = sectionDef._parentSections || []
    sectionDef._qualifiedName = `${pkg.name}:${sectionDef.name}`

    const addPropertiesFromSections = sections => sections
      .map(ref => resolveRef(ref)).forEach(extendingSectionDef => {
        if (extendingSectionDef.quantities) {
          sectionDef.quantities.push(...extendingSectionDef.quantities)
        }
        if (extendingSectionDef.sub_sections) {
          sectionDef.sub_sections.push(...extendingSectionDef.sub_sections)
        }
      })
    sectionDef.extending_sections = sectionDef.extending_sections || []
    addPropertiesFromSections(sectionDef.extending_sections)
    if (!sectionDef.extends_base_section) {
      addPropertiesFromSections(sectionDef.base_sections)
      sectionDef.base_sections = sectionDef.base_sections || []
    }

    sectionDef._properties = {}
    const addProperty = property => {
      sectionDef._properties[property.name] = property
      if (!sectionDef.extends_base_section) {
        property._section = sectionDef
        property._qualifiedName = `${sectionDef._qualifiedName}:${property.name}`
      }
    }
    sectionDef.quantities.forEach(quantitiy => {
      addProperty(quantitiy)
      quantitiy.shape = quantitiy.shape || []
      if (isReference(quantitiy)) {
        const referencedSection = resolveRef(quantitiy.type.type_data)
        referencedSection._incomingRefs = referencedSection._incomingRefs || []
        referencedSection._incomingRefs.push(quantitiy)
      }
    })
    sectionDef.sub_sections.forEach(subSection => {
      addProperty(subSection)
      const subSectionsSectionDef = resolveRef(subSection.sub_section)
      subSectionsSectionDef._parentSections = subSectionsSectionDef._parentSections || []
      subSectionsSectionDef._parentSections.push(sectionDef)
    })
  })
})

export function resolveRef(ref, data) {
  data = data || metainfo
  const segments = ref.split('/').filter(segment => segment !== '')
  const reducer = (current, segment) => {
    return isNaN(segment) ? current[segment] : current[parseInt(segment)]
  }
  return segments.reduce(reducer, data)
}

export function metainfoDef(name) {
  return packageDefs['nomad.metainfo.metainfo']._sections[name]
}

export function isReference(property) {
  return property.type && property.type.type_kind === 'reference'
}

/**
 * Constructs a graph from and with the definition. The graph will contain the given nodes,
 * all its outgoing and incomming references, the parents up to root (for sections and categories)
 *
 * @param {Object} def
 */
export function vicinityGraph(def) {
  const nodesMap = {}
  const nodes = []
  const edges = []

  function addEdge(from, to, def) {
    const edge = {
      def: def,
      source: from.id,
      target: to.id,
      value: 1
    }
    edges.push(edge)
  }

  function addNode(def, more) {
    const {x, y, recursive} = more || {}
    if (nodesMap[def._qualifiedName]) {
      console.log('####### B', def._qualifiedName)
      return nodesMap[def._qualifiedName]
    }

    const node = {
      id: def._qualifiedName,
      def: def,
      x: x || 0,
      y: y || 0
    }
    nodes.push(node)
    nodesMap[def._qualifiedName] = node

    if (recursive) {
      if (def.m_def === 'Section') {
          def._parentSections.forEach(parentSection => {
            const parent = addNode(parentSection, {recursive: true})
            addEdge(node, parent, {})
          })
          def.quantities
            .filter(quantity => quantity.type.type_kind === 'reference')
            .forEach(reference => {
              const referencedSectionDef = resolveRef(reference.type.type_data)
              const referenced = addNode(referencedSectionDef)
              addEdge(node, referenced, reference)
            })
      } else if (def.m_def === 'Quantity') {
        const section = addNode(def._section, {recursive: true})
        addEdge(node, section, {})
      }

      if (def.categories) {
        def.categories.forEach(categoryDef => {
          const category = addNode(resolveRef(categoryDef), {recursive: true})
          addEdge(node, category, {})
        })
      }
    }

    return node
  }

  addNode(def, {recursive: true})

  return {
    nodes: nodes,
    links: edges
  }
}

export default metainfo
