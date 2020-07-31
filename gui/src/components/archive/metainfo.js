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
  const dx = 100
  const dy = 75

  function addEdge(from, to, def) {
    const edge = {
      def: def,
      source: from,
      target: to,
      value: 1
    }
    edges.push(edge)
  }

  function addNode(def, more, id) {
    id = id || (d => d._qualifiedName)
    const {recursive, x, y, i} = more || {}

    const key = id(def)
    if (nodesMap[key]) {
      return nodesMap[key]
    }

    const node = {
      id: key,
      def: def,
      x: x, y: y, i: i
    }

    nodes.push(node)
    nodesMap[key] = node

    if (recursive) {
      if (def.m_def === 'Section') {
          def._parentSections.forEach(parentSection => {
            const parent = addNode(parentSection, {
              recursive: true,
              x: x - dx, y: y, i: i + 1})
            addEdge(node, parent, {})
          })
          const references = def.quantities.filter(quantity => quantity.type.type_kind === 'reference')
          const layoutMiddle = (references.length - 1) * dx / 2
          references.forEach((reference, i) => {
              const referencedSectionDef = resolveRef(reference.type.type_data)
              const referenced = addNode(
                referencedSectionDef,
                {x: x + i * dx - layoutMiddle, y: y + dy, i: i},
                () => reference._qualifiedName)
              addEdge(node, referenced, reference)
            })
      } else if (def.m_def === 'Quantity') {
        const section = addNode(def._section, {
          recursive: true,
          x: x - dx, y: y, i: i + 1})
        addEdge(node, section, {})
      }

      if (def.categories) {
        const layoutMiddle = (def.categories.length - 1) * dx / 2
        def.categories.forEach((categoryDef, i) => {
          const category = addNode(resolveRef(categoryDef), {
            recursive: true,
            x: x + i * dx - layoutMiddle, y: y - dy, i: i
          })
          addEdge(node, category, {})
        })
      }
    }

    return node
  }

  addNode(def, {recursive: true, x: 0, y: 0, i: 0})

  return {
    nodes: nodes,
    links: edges
  }
}

export default metainfo
