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
    const addProperty = property => { sectionDef._properties[property.name] = property }
    sectionDef.quantities.forEach(quantitiy => {
      addProperty(quantitiy)
      quantitiy.shape = quantitiy.shape || []
    })
    sectionDef.sub_sections.forEach(addProperty)
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

export default metainfo
