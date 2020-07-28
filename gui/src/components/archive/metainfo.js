import metainfo from '../../metainfo'

export const sectionDefs = {}
export const packageDefs = {}

metainfo.packages.forEach(pkg => {
  packageDefs[pkg.name] = pkg
  pkg._sections = {}
  pkg.section_definitions.forEach(sectionDef => {
    pkg._sections[sectionDef.name] = sectionDef
    sectionDefs[sectionDef.name] = sectionDef
    sectionDef.quantities = sectionDef.quantities || []
    sectionDef.sub_sections = sectionDef.sub_sections || []
    sectionDef._incomingRefs = []
    sectionDef._parentSections = []
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
      from: from,
      to: to
    }
    edges.push(edge)
  }

  function addNode(def, x, y) {
    if (nodesMap[def._qualifiedName]) {
      return nodesMap[def._qualifiedName]
    }

    const node = {
      def: def,
      x: x || 0,
      y: y || 0
    }
    node.index = nodes.push(node)
    nodesMap[def._qualifiedName] = node

    if (def.m_def === 'Section') {
      def._parentSections.forEach(parentSection => {
        const parentIndex = addNode(parentSection)
        addEdge(node.index, parentIndex, {})
      })
      def.quantities
        .filter(quantity => quantity.type.type_kind === 'reference')
        .forEach(reference => {
          const referencedSectionDef = resolveRef(reference.type_data)
          const index = addNode(referencedSectionDef)
          addEdge(node.index, index, reference)
        })
    } else if (def.m_def === 'Quantity') {
      const sectionIndex = addNode(def._section)
      addEdge(node.index, sectionIndex, {})
    }

    if (def.categories) {
      def.categories.forEach(category => {
        const categoryIndex = addNode(category)
        addEdge(node.index, categoryIndex, {})
      })
    }

    return node
  }

  addNode(def)

  return {
    nodes: nodes,
    edges: edges
  }
}

export default metainfo
