/*
 * Copyright The NOMAD Authors.
 *
 * This file is part of NOMAD. See https://nomad-lab.eu for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
import metainfo from '../../metainfo'

export const SectionMDef = 'nomad.metainfo.metainfo.Section'
export const QuantityMDef = 'nomad.metainfo.metainfo.Quantity'
export const SubSectionMDef = 'nomad.metainfo.metainfo.SubSection'
export const CategoryMDef = 'nomad.metainfo.metainfo.Category'

export const defs = []
export const defsByName = {}
export const packageDefs = {}
export const packagePrefixes = {}

const addDef = def => {
  const defsForName = defsByName[def.name] || []
  defsByName[def.name] = defsForName
  if (!defsForName.find(existing => existing === def)) {
    defsForName.push(def)
  }
  defs.push(def)
  def.more = def.more || {}
}

function sortDefs(defs) {
  return defs.sort((a, b) => a.name.localeCompare(b.name))
}

metainfo.packages.forEach(pkg => {
  packageDefs[pkg.name] = pkg
  const prefix = pkg.name.split('.')[0]
  packagePrefixes[prefix] = packagePrefixes[prefix] || {}
  packagePrefixes[prefix][pkg.name] = pkg

  pkg._sections = {}
  pkg.category_definitions = pkg.category_definitions || []
  pkg.section_definitions = pkg.section_definitions || []
  pkg.category_definitions.forEach(categoryDef => {
    categoryDef._qualifiedName = `${pkg.name}.${categoryDef.name}`
    categoryDef._package = pkg
    addDef(categoryDef)
  })

  const initSection = sectionDef => {
    sectionDef.base_sections = sectionDef.base_sections || []
    sectionDef.quantities = sectionDef.quantities || []
    sectionDef.sub_sections = sectionDef.sub_sections || []
    sectionDef.inner_section_definitions = sectionDef.inner_section_definitions || []
    sectionDef._incomingRefs = sectionDef._incomingRefs || []
    sectionDef._parentSections = sectionDef._parentSections || []
    sectionDef._parentSubSections = sectionDef._parentSubSections || []
  }

  const computeAllProperties = sectionDef => {
    const results = {}
    function addProperties(sectionDef) {
      sectionDef.quantities.concat(sectionDef.sub_sections).forEach(
        property => {
          results[property.name] = property
        })
    }
    const baseSections = []
    function addBaseSection(sectionDef) {
      initSection(sectionDef)
      if (!sectionDef.extends_base_section) {
        sectionDef.base_sections.map(baseSectionRef => resolveRef(baseSectionRef)).forEach(addBaseSection)
      }
      baseSections.push(sectionDef)
    }
    addBaseSection(sectionDef)
    baseSections.forEach(addProperties)
    return Object.keys(results).map(key => results[key])
  }

  const addSectionDef = (sectionDef, parentDef) => {
    sectionDef._parent = parentDef
    pkg._sections[sectionDef.name] = sectionDef
    initSection(sectionDef)
    sectionDef._qualifiedName = parentDef ? `${parentDef._qualifiedName || parentDef.name}.${sectionDef.name}` : sectionDef.name
    sectionDef._package = pkg

    sectionDef.inner_section_definitions.forEach(innerSectionDef => addSectionDef(innerSectionDef, sectionDef))

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
      addDef(sectionDef)
    }

    sectionDef._allProperties = computeAllProperties(sectionDef)
    sectionDef._properties = {}
    const addProperty = property => {
      sectionDef._properties[property.name] = property
      if (!sectionDef.extends_base_section) {
        property._section = sectionDef
        property._qualifiedName = `${sectionDef._qualifiedName}.${property.name}`
      }
      property._package = pkg
      addDef(property)
    }
    sectionDef._allProperties.forEach(property => {
      addProperty(property)
      if (property.m_def === QuantityMDef) {
        property.shape = property.shape || []
        if (isReference(property)) {
          const referencedSection = resolveRef(property.type.type_data)
          referencedSection._incomingRefs = referencedSection._incomingRefs || []
          referencedSection._incomingRefs.push(property)
        }
      } else if (property.m_def === SubSectionMDef) {
        property._subSection = resolveRef(property.sub_section)
        const subSectionsSectionDef = property._subSection
        initSection(subSectionsSectionDef)
        subSectionsSectionDef._parentSections.push(sectionDef)
        subSectionsSectionDef._parentSubSections.push(property)
        property._section = sectionDef
      }
    })
  }

  pkg.section_definitions.forEach(sectionDef => addSectionDef(sectionDef, pkg))
})

export const rootSections = sortDefs(defs.filter(def => (
  def.m_def === SectionMDef && !def.extends_base_section && def._parentSections.length === 0)
))

export async function resolveRefAsync(ref, data, fetchArchive) {
  if (!data) {
    return resolveRef(ref)
  }

  if (!ref.includes('#')) {
    return resolveRef(ref, data)
  }

  const [url] = ref.split('#')
  if (url === '') {
    return resolveRef(ref, data)
  }

  if (data?.resources[url]) {
    return resolveRef(ref, data)
  }

  const uploadId = data?.metadata?.upload_id
  if (!(url.startsWith('../upload/archive/') && uploadId)) {
    return null
  }

  const archive = await fetchArchive(`uploads/${uploadId}/${url.slice('../upload/'.length)}`)
  if (!data.resources) {
    data.resource = {}
  }
  data.resources[url] = archive

  resolveRef(ref, data)
}

/**
 * Resolves the given string reference into the actual data.
 * @param {string} ref Reference.
 * @param {object} data Archive.
 */
export function resolveRef(ref, data) {
  if (!ref) {
    return null
  }

  const resolve = (ref, context) => {
    const parts = ref.split('#')
    if (parts.length === 2 && parts[0] !== '') {
      const url = parts[0]
      ref = parts[1]
      data = data?.resources[url]
      if (!data) {
        return null
      }
    } else {
      if (parts.length === 2) {
        ref = parts[1]
      } else {
        ref = parts[0]
      }
    }

    try {
      context = data || metainfo
      const segments = ref.split('/').filter(segment => segment !== '')
      const reducer = (current, segment) => {
        return isNaN(segment) ? current[segment] : current[parseInt(segment)]
      }
      return segments.reduce(reducer, context)
    } catch (e) {
      console.log('could not resolve: ' + ref)
      throw e
    }
  }
  if (Array.isArray(ref)) {
    return ref.map(x => resolve(x, data))
  }
  return resolve(ref, data)
}

/**
 * Converts a reference given in the /section/<index>/subsection format (used
 * in the metainfo) to the /section:<index>/subsection format (used by the
 * archive browser).
 * @param {*} ref The reference to convert.
 */
export function refPath(ref) {
  try {
    const segments = ref.split('/').filter(segment => segment !== '')
    const reducer = (current, segment) => {
      return isNaN(segment) ? `${current}/${segment}` : `${current}:${segment}`
    }
    return segments.reduce(reducer)
  } catch (e) {
    console.log('could not convert the path: ' + ref)
    throw e
  }
}

export function metainfoDef(qualifiedName) {
  const defQualifiedNameSegments = qualifiedName.split('.')
  const packageName = defQualifiedNameSegments.slice(0, -1).join('.')
  const sectionName = defQualifiedNameSegments[defQualifiedNameSegments.length - 1]
  return packageDefs[packageName]._sections[sectionName]
}

export function isReference(property) {
  return property.type && property.type.type_kind === 'reference'
}

export function path(nameOrDef) {
  let def
  if (typeof nameOrDef === 'string') {
    def = defsByName[nameOrDef] && defsByName[nameOrDef].find(def => true)
  } else {
    def = nameOrDef
  }

  if (!def) {
    return null
  }

  if (def.m_def === SubSectionMDef) {
    def = resolveRef(def.sub_section)
  }

  if (def.m_def === CategoryMDef) {
    return `${def._package.name.split('.')[0]}/category_definitions@${def._qualifiedName}`
  }

  const path = []
  while (def) {
    const parentSection = def
    const parentSubSection = def._parentSubSections && def._parentSubSections.filter(
      // Filter for direct recursions in the possible section containment.
      // This only catches direct connections where a sub section uses its parent
      // section as the sub section definition
      subSection => parentSection !== subSection._section)[0]
    if (parentSubSection) {
      def = parentSubSection
    }
    path.push(def.name)
    if (def.m_def === SubSectionMDef) {
      def = def._section
    } else {
      def = def._parentSections && def._parentSections[0]
    }

    while (def && def.extends_base_section) {
      def = resolveRef(def.base_sections[0])
    }
  }
  return path.reverse().join('/')
}

/**
 * @param {*} def A section definition.
 * @returns True, if sections of the given section def are editable.
 */
export function isEditable(def) {
  return !!def._allProperties.find(prop => prop.m_annotations?.eln) || !!def.m_annotations?.eln
}

export function removeSubSection(section, subSectionDef, index) {
  if (subSectionDef.repeats) {
    section[subSectionDef.name].splice(index, 1)
  } else {
    section[subSectionDef.name] = undefined
    delete section[subSectionDef.name]
  }
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
      x: x,
      y: y,
      i: i
    }

    nodes.push(node)
    nodesMap[key] = node

    if (recursive) {
      if (def.m_def === SectionMDef) {
        def._parentSections.forEach(parentSection => {
          const parent = addNode(parentSection, {
            recursive: true,
            x: x - dx,
            y: y,
            i: i + 1})
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
      } else if (def.m_def === QuantityMDef) {
        const section = addNode(def._section, {
          recursive: true,
          x: x - dx,
          y: y,
          i: i + 1})
        addEdge(node, section, {})
      }

      if (def.categories) {
        const layoutMiddle = (def.categories.length - 1) * dx / 2
        def.categories.forEach((categoryDef, i) => {
          const category = addNode(resolveRef(categoryDef), {
            recursive: true,
            x: x + i * dx - layoutMiddle,
            y: y - dy,
            i: i
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
