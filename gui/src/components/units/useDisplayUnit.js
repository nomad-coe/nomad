import {useErrors} from "../errors"
import {useMemo} from "react"
import {Unit} from "./Unit"
import {useUnitContext} from "./UnitContext"
import {getFieldProps} from "../editQuantity/StringEditQuantity"

function getUnitSystem(def) {
  let unit_system
  if (def?.m_def === "nomad.metainfo.metainfo.Quantity") {
    unit_system = def?._parent?.m_annotations?.display?.[0]?.unit_system
    if (!unit_system && def?._section?._parentSections?.[0]) {
      unit_system = getUnitSystem(def._section._parentSections[0])
    }
  } else if (def?.m_def === "nomad.metainfo.metainfo.Section") {
    unit_system = def?.m_annotations?.display?.[0]?.unit_system
    if (!unit_system && def?._parentSections?.[0]) {
      unit_system = getUnitSystem(def._parentSections[0])
    }
  }
  return unit_system
}

export function useDisplayUnit(quantityDef) {
  const {units, unitSystems} = useUnitContext()
  const {raiseError} = useErrors()
  const defaultUnit = useMemo(() => quantityDef.unit && new Unit(quantityDef.unit), [quantityDef])
  const dimension = defaultUnit && defaultUnit.dimension(false)
  const {defaultDisplayUnit: deprecatedDefaultDisplayUnit} = getFieldProps(quantityDef)
  const defaultDisplayUnit = quantityDef?.m_annotations?.display?.[0]?.unit || deprecatedDefaultDisplayUnit

  const displayUnitObj = useMemo(() => {
    if (!dimension) return
    let defaultDisplayUnitObj
    const section_default_unit_system = getUnitSystem(quantityDef)

    // TODO: If we enable the new 'Schema' scope in the unit context, we should
    // prioritize those values there. But for now we just read unit info from
    // the schema.

    // If a default display unit has been defined, use it instead
    if (defaultDisplayUnit) {
      try {
        defaultDisplayUnitObj = new Unit(defaultDisplayUnit)
      } catch (e) {
        raiseError(`The provided defaultDisplayUnit for ${quantityDef.name} field is not valid.`)
      }
      if (defaultDisplayUnitObj.dimension(true) !== defaultUnit.dimension(true)) {
        raiseError(`The provided defaultDisplayUnit for ${quantityDef.name} has incorrect dimensionality for this field.`)
      }
    // If a unit system has been defined, use it
    } else if (section_default_unit_system && unitSystems[section_default_unit_system]) {
      defaultDisplayUnitObj = new Unit(unitSystems[section_default_unit_system].units[dimension].definition)
    // Fallback option is the unit system in the global unit scope
    } else {
      defaultDisplayUnitObj = new Unit(defaultUnit).toSystem(units)
    }

    return defaultDisplayUnitObj
  }, [defaultDisplayUnit, defaultUnit, dimension, quantityDef, raiseError, unitSystems, units])

  return displayUnitObj
}
