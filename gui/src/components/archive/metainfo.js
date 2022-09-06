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

import React, { useCallback, useContext, useEffect, useState } from 'react'
import PropTypes from 'prop-types'
import metainfoData from '../../metainfo'
import { useApi } from '../api'

const metainfoContext = React.createContext()

export const GlobalMetainfo = React.memo(function GlobalMetainfo({children}) {
  const [globalMetainfo, setGlobalMetainfo] = useState(metainfoData._metainfo)
  useEffect(() => {
    if (!globalMetainfo) {
      createGlobalMetainfo().then(metainfo => {
        setGlobalMetainfo(metainfo)
      })
    }
  }, [setGlobalMetainfo, globalMetainfo])

  const {api} = useApi()
  const [allCustomMetainfos, setAllCustomMetainfos] = useState()

  const fetchAllCustomMetainfos = useCallback(async (refresh, query) => {
    if (allCustomMetainfos && !refresh) {
      return allCustomMetainfos
    }

    // TODO paginate?
    // TODO Only grab new ones?
    query = query || {}
    const response = await api.post(`entries/archive/query`, {
      owner: 'visible',
      query: {
        ...query,
        quantities: 'definitions.section_definitions',
        processed: true
      },
      required: {
        definitions: '*',
        metadata: {
          mainfile: '*',
          entry_name: '*'
        }
      }
    })
    for (const data of response.data) {
      const archive = data.archive
      // TODO we should not createMetainfo all the time? There needs to be a
      // register/cache based on hashes or something
      archive._metainfo = await createMetainfo(archive, globalMetainfo, {api: api, archive: archive})
    }
    setAllCustomMetainfos(response.data)
    return response.data
  }, [globalMetainfo, api, allCustomMetainfos, setAllCustomMetainfos])

  if (globalMetainfo) {
    globalMetainfo.fetchAllCustomMetainfos = fetchAllCustomMetainfos
  }
  return <metainfoContext.Provider value={globalMetainfo}>
    {children}
  </metainfoContext.Provider>
})
GlobalMetainfo.propTypes = {
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
}

export function useGlobalMetainfo() {
  return useContext(metainfoContext)
}

export function useMetainfo(data) {
  const [metainfo, setMetainfo] = useState()
  const globalMetainfo = useGlobalMetainfo()
  useEffect(() => {
    if (data) {
      createMetainfo(data, globalMetainfo).then(setMetainfo)
    } else if (globalMetainfo) {
      setMetainfo(globalMetainfo)
    }
  }, [data, globalMetainfo, setMetainfo])

  return metainfo
}

export const PackageMDef = 'nomad.metainfo.metainfo.Package'
export const SectionMDef = 'nomad.metainfo.metainfo.Section'
export const QuantityMDef = 'nomad.metainfo.metainfo.Quantity'
export const SubSectionMDef = 'nomad.metainfo.metainfo.SubSection'
export const CategoryMDef = 'nomad.metainfo.metainfo.Category'
export const AttributeMDef = 'nomad.metainfo.metainfo.Attribute'

export async function createGlobalMetainfo() {
  return createMetainfo(metainfoData)
}

export async function createMetainfo(data, parentMetainfo, context) {
  if (!(data.packages || data.definitions) && parentMetainfo) {
    return parentMetainfo
  }

  if (data._metainfo) {
    return data._metainfo
  }

  const metainfo = new Metainfo(data, parentMetainfo, context)
  if (data.packages) {
    await metainfo._addPackages(data.packages)
  }
  if (data.definitions) {
    const entryId = data?.metadata?.entry_id
    const mainfile = data?.metadata?.mainfile
    const uploadId = data?.metadata?.upload_id
    const url = mainfile && uploadId && `../uploads/${uploadId}/raw/${mainfile}#definitions`
    await metainfo._addPackages([data.definitions], entryId ? `entry_id:${entryId}` : null, url)
  }
  data._metainfo = metainfo
  return metainfo
}

/**
 * Represents and manages schema data.
 *
 * It is responsible to create derived definition properties that are necessary for the
 * browser's function and to allow to resolve references to definitions.
 *
 * It allows to access schema data in two ways: via references and via qualified names
 * (e.g. based on packages name + section name).
 *
 * For package-based access, metainfo instances can be linked. A metainfo can have a parent.
 * Here metainfo instances will first try to resolve a qualified names with itself and it's
 * parent.
 *
 * For reference-based access, contexts can be used. A metainfo instance represents a single
 * root section (e.g. Package or Environment). URL fragments can be resolved here. The
 * rest of the reference URL can be used load respective resources with the context.
 */
class Metainfo {
  constructor(data, parent, context) {
    this._context = context
    this._parent = parent
    this._data = data
    this._defs = new Set()

    this._packageDefs = {}
    this._defsByNameCache = null
    this._packagePrefixCache = null
    this._rootSectionsCache = null
  }

  /**
   * @returns All definitions as an Array.
   */
  getDefs() {
    if (this._parent) {
      return this._parent.getDefs().concat([...this._defs])
    } else {
      return [...this._defs]
    }
  }

  /**
   * @returns An object with all definition names as keys. Packages are excluded.
   *   The values are arrays with all the definitions that share the name.
   */
  getDefsByName() {
    if (this._defsByNameCache) {
      return this._defsByNameCache
    }
    this._defsByNameCache = this.getDefs().filter(def => def.m_def !== PackageMDef).reduce((result, def) => {
      result[def.name] = result[def.name] || []
      result[def.name].push(def)
      return result
    }, {})
    return this._defsByNameCache
  }

  /**
   * @returns An object with all package prefixes (the part of the name before the first ".")
   *   The values are arrays with all the packages that share the prefix.
   */
  getPackagePrefixes() {
    if (this._packagePrefixCache) {
      return this._packagePrefixCache
    }
    this._packagePrefixCache = this.getDefs()
      .filter(def => def.m_def === PackageMDef)
      .reduce((results, pkg) => {
        const packageName = pkg.name
        if (packageName !== '*') {
          const prefix = packageName.split('.')[0]
          results[prefix] = results[prefix] || {}
          results[prefix][packageName] = pkg
        }
        return results
      }, {})
    return this._packagePrefixCache
  }

  /**
   * @returns An array with all section definitions where not sub-section is using them
   *   as section definition.
   */
  getRootSectionDefinitions() {
    if (this._rootSectionsCache) {
      return this._rootSectionsCache
    }
    const sortDefs = defs => defs.sort((a, b) => a.name.localeCompare(b.name))
    this._rootSectionsCache = sortDefs(this.getDefs().filter(def => (
      def.m_def === SectionMDef && !def.extends_base_section && def._parentSections.length === 0)
    ))
    return this._rootSectionsCache
  }

  getEntryArchiveDefinition() {
    // TODO this is a super wasteful implementation
    const entryArchiveDefinition = this.getRootSectionDefinitions().find(def => def.name === 'EntryArchive')
    return entryArchiveDefinition
  }

  async _addDef(def) {
    // only add if not already added
    if (this._defs.has(def)) {
      return
    }
    def.more = def.more || {}
    def.categories = await this.resolveDefinitionList(def.categories || [])
    this._defs.add(def)
    this._defsByNameCache = null
  }

  async _initSection(sectionDef) {
    if (sectionDef._allBaseSections !== undefined) {
      return sectionDef
    }

    sectionDef.base_sections = await this.resolveDefinitionList(sectionDef.base_sections || [])
    sectionDef.quantities = sectionDef.quantities || []
    sectionDef.sub_sections = sectionDef.sub_sections || []
    sectionDef.inner_section_definitions = sectionDef.inner_section_definitions || []
    sectionDef._allBaseSections = []
    sectionDef._allInheritingSections = []
    sectionDef._incomingRefs = []
    sectionDef._parentSections = []
    sectionDef._parentSubSections = []

    if (!sectionDef.extends_base_section) {
      for (const baseSection of sectionDef.base_sections) {
        await this._initSection(baseSection)
        sectionDef._allBaseSections.push(baseSection)
        baseSection._allBaseSections.forEach(baseBaseSection => sectionDef._allBaseSections.push(baseBaseSection))
        if (!baseSection._allInheritingSections.includes(sectionDef)) {
          baseSection._allInheritingSections.push(sectionDef)
          sectionDef._allInheritingSections.forEach(inheritingInheritingSection => {
            if (!baseSection._allInheritingSections.includes(inheritingInheritingSection)) {
              baseSection._allInheritingSections.push(inheritingInheritingSection)
            }
          })
        }
      }
    }

    return sectionDef
  }

  async _getAllProperties(sectionDef) {
    const results = {}
    function createAddProperties(inherited) {
      return (sectionDef) => {
        function createAddProperty(m_def) {
          return (property) => {
            const propertyToAdd = inherited ? {} : property
            if (inherited) {
              Object.assign(propertyToAdd, property)
              propertyToAdd._inherited = true
            } else {
              if (results[property.name]) {
                propertyToAdd._overwritten = true
              }
            }
            property.m_def = m_def
            results[property.name] = propertyToAdd
          }
        }
        sectionDef.quantities.forEach(createAddProperty(QuantityMDef))
        sectionDef.sub_sections.forEach(createAddProperty(SubSectionMDef))
      }
    }
    sectionDef = await this._initSection(sectionDef)
    sectionDef._allBaseSections.forEach(createAddProperties(true))
    createAddProperties(false)(sectionDef)
    return Object.keys(results).map(key => results[key])
  }

  async _addSection(pkg, sectionDef, parentDef, parentProperty, parentIndex) {
    this._rootSectionsCache = null
    sectionDef.m_def = SectionMDef
    sectionDef._parent = parentDef
    sectionDef._parentProperty = parentProperty
    sectionDef._parentIndex = parentIndex
    pkg._sections[sectionDef.name] = sectionDef
    await this._initSection(sectionDef)
    sectionDef._qualifiedName = parentDef ? `${parentDef._qualifiedName || parentDef._unique_id || parentDef.name}.${sectionDef.name}` : sectionDef.name
    if (parentDef?._url) {
      sectionDef._url = `${parentDef._url}/section_definitions/${parentIndex}`
    }
    sectionDef._package = pkg

    let index = 0
    for (const innerSectionDef of sectionDef.inner_section_definitions) {
      await this._addSection(pkg, innerSectionDef, sectionDef, 'inner_section_definitions', index)
      index++
    }

    const addPropertiesFromSections = async sections => {
      for (const ref of sections) {
        const extendingSectionDef = await this.resolveDefinition(ref)
        if (extendingSectionDef.quantities) {
          sectionDef.quantities.push(...extendingSectionDef.quantities)
        }
        if (extendingSectionDef.sub_sections) {
          sectionDef.sub_sections.push(...extendingSectionDef.sub_sections)
        }
      }
    }
    sectionDef.extending_sections = sectionDef.extending_sections || []
    await addPropertiesFromSections(sectionDef.extending_sections)
    if (!sectionDef.extends_base_section) {
      await this._addDef(sectionDef)
    }

    sectionDef._allProperties = await this._getAllProperties(sectionDef)
    sectionDef._properties = {}
    const addProperty = async property => {
      sectionDef._properties[property.name] = property
      if (!sectionDef.extends_base_section) {
        property._section = sectionDef
        property._qualifiedName = `${sectionDef._qualifiedName}.${property.name}`
      }
      property._package = pkg
      await this._addDef(property)
    }
    for (const property of sectionDef._allProperties) {
      await addProperty(property)
      if (property.m_def === QuantityMDef) {
        property.shape = property.shape || []
        if (isReference(property)) {
          const referencedSection = await this.resolveDefinition(property.type.type_data)
          referencedSection._incomingRefs = referencedSection._incomingRefs || []
          referencedSection._incomingRefs.push(property)
          property.type._referencedSection = referencedSection
        }
      } else if (property.m_def === SubSectionMDef) {
        property.sub_section = await this.resolveDefinition(property.sub_section)
        const subSectionsSectionDef = property.sub_section
        await this._initSection(subSectionsSectionDef)
        subSectionsSectionDef._parentSections.push(sectionDef)
        subSectionsSectionDef._parentSubSections.push(property)
        property._section = sectionDef
      }
    }
  }

  async _addPackage(pkg, unique_id, url) {
    this._packagePrefixCache = null
    pkg.m_def = PackageMDef
    if (unique_id) {
      pkg._unique_id = unique_id
    }
    pkg._url = url
    const packageName = pkg.name || '*'
    this._packageDefs[packageName] = pkg
    await this._addDef(pkg)

    pkg._sections = {}
    pkg.category_definitions = pkg.category_definitions || []
    pkg.section_definitions = pkg.section_definitions || []
    for (const categoryDef of pkg.category_definitions) {
      categoryDef.m_def = CategoryMDef
      categoryDef._qualifiedName = `${pkg._unique_id || pkg.name}.${categoryDef.name}`
      categoryDef._package = pkg
      await this._addDef(categoryDef)
    }

    let index = 0
    for (const sectionDef of pkg.section_definitions) {
      await this._addSection(pkg, sectionDef, pkg, 'section_definitions', index++)
    }
  }

  async _addPackages(packages, unique_id, url) {
    for (const pkg of packages) {
      await this._addPackage(pkg, unique_id, url)
    }
  }

  path(nameOrDef) {
    let def
    if (typeof nameOrDef === 'string') {
      const defsByName = this.getDefsByName()
      def = defsByName[nameOrDef] && defsByName[nameOrDef].find(def => true)
    } else {
      def = nameOrDef
    }

    if (!def) {
      return null
    }

    if (def.m_def === SubSectionMDef) {
      def = def.sub_section
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
        def = def.base_sections[0]
      }
    }
    return path.reverse().join('/')
  }

  async resolveDefinitionList(references, context) {
    const result = []
    for (const reference of references) {
      result.push(await this.resolveDefinition(reference, context))
    }
    return result
  }

  async resolveDefinition(reference, context) {
    context = context || this._context
    if (typeof reference !== 'string') {
      // already resolved
      return reference
    }
    if (reference.match(/.+#.*/)) {
      if (!context) {
        console.error('cannot resolve definition without context', reference)
        return null
      }
      return resolveRefAsync(reference, this._data, context, async archive => {
        if (!archive._metainfo) {
          archive._metainfo = await createMetainfo(archive, context.metainfo, context)
        }
      })
    }

    const resolved = await this._parent?.resolveDefinition(reference, context)
    if (resolved) {
      return resolved
    }

    if (reference.includes('#') || reference.includes('/')) {
      return resolveRef(reference, this._data)
    } else {
      const qualifiedName = reference
      const defQualifiedNameSegments = qualifiedName.split('.')
      const packageName = defQualifiedNameSegments.slice(0, -1).join('.')
      const sectionName = defQualifiedNameSegments[defQualifiedNameSegments.length - 1]
      return this._packageDefs[packageName]?._sections?.[sectionName]
    }
  }
}

export async function resolveRefAsync(reference, data, context, adaptArchive) {
  if (!data) {
    data = context?.archive
  }

  if (!data) {
    return resolveRef(reference)
  }

  if (!reference.includes('#') || reference.startsWith('#')) {
    return resolveRef(reference, data)
  }

  let resourceUrl = reference.slice(0, reference.indexOf('#'))
  reference = reference.slice(reference.indexOf('#'))
  if (!context.resources[resourceUrl]) {
    try {
      let apiUrl
      let uploadId = context.uploadId

      const uploadIdMatch = resourceUrl.match(/\.\.\/uploads\/([a-zA-Z0-9_-]+)\//)
      if (uploadIdMatch) {
        uploadId = uploadIdMatch[1]
        resourceUrl = '../upload/' + resourceUrl.slice(uploadIdMatch[0].length)
      }

      if (resourceUrl.startsWith('../upload/archive')) {
        apiUrl = `uploads/${uploadId}/${resourceUrl.slice('../upload/'.length)}`
      } else if (resourceUrl.startsWith('../upload/raw')) {
        const mainfile = resourceUrl.slice('../upload/raw/'.length)
        const queryBody = ({
          owner: 'visible',
          query: {
            upload_id: uploadId,
            'mainfile': mainfile
          },
          required: {
            include: ['entry_id']
          }
        })
        const queryResponse = await context.api.post(`/entries/query`, queryBody)
        if (!queryResponse.data[0]) {
          return null
        }
        const entryId = queryResponse.data[0].entry_id
        apiUrl = `uploads/${uploadId}/archive/${entryId}`
      } else {
        console.error(`Reference solutions for urls like ${resourceUrl} is not yet implemented.`)
        return null
      }

      const response = await context.api.get(apiUrl)
      const archive = response.data.archive
      context.resources[resourceUrl] = archive
    } catch {
      return null
    }
  }
  data = context.resources[resourceUrl]
  if (adaptArchive) {
    await adaptArchive(data)
  }
  if (!data) {
    return null
  }

  return resolveRef(reference, data)
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
      data = data?.resources?.[url]
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
      context = data || metainfoData
      const segments = ref.split('/').filter(segment => segment !== '')
      const reducer = (current, segment) => {
        return isNaN(segment)
          ? current?.[segment] || current?._properties?.[segment] || current?.sub_section?._properties?.[segment]
          : current?.repeats && current?.m_def
            ? current
            : current?.[parseInt(segment)]
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
    if (ref.startsWith('#')) {
      ref = ref.slice(1)
    }
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

export function isReference(property) {
  return property.type && property.type.type_kind === 'reference'
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
 * @param {*} definition The section definition to create a reference for.
 * @returns The reference fragment for the given section definition.
 */
export function getSectionReference(definition) {
  if (!definition._parent) {
    return ''
  }
  const ref = `${getSectionReference(definition._parent)}/${definition._parentProperty}`
  if (!isNaN(definition._parentIndex)) {
    return `${ref}/${definition._parentIndex}`
  }
  return ref
}

/**
 * Allows to traverse a given section through all its sub-sections.
 * @param {Object} section The section to traverse
 * @param {Object} definition The definition of the section
 * @param {str} path The archive path to the section
 * @param {function} callback The callback that is called on each section
 */
export function traverse(section, definition, path, callback, depthFirst) {
  callback(section, definition, path)
  for (const subSectionDef of definition._allProperties.filter(prop => prop.m_def === SubSectionMDef)) {
    let subSections = []
    if (!subSectionDef.repeats) {
      const subSection = section[subSectionDef.name]
      if (subSection) {
        subSections = [subSection]
      }
    } else {
      subSections = section[subSectionDef.name] || []
    }

    subSections.forEach((subSection, index) => {
      let childPath = `${path}/${subSectionDef.name}`
      if (subSectionDef.repeats) {
        childPath = `${childPath}/${index}`
      }
      traverse(subSection, subSectionDef.sub_section, childPath, callback)
    })
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
          const category = addNode(categoryDef, {
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
