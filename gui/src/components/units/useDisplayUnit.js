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

export function useDisplayUnit(quantity, quantityDef) {
  const {units, isReset, unitSystems} = useUnitContext()
  const {raiseError} = useErrors()
  const defaultUnit = useMemo(() => quantityDef.unit && new Unit(quantityDef.unit), [quantityDef])
  const dimension = defaultUnit && defaultUnit.dimension(false)
  const {defaultDisplayUnit: deprecatedDefaultDisplayUnit} = getFieldProps(quantityDef)
  const defaultDisplayUnit = quantityDef?.m_annotations?.display?.[0]?.unit || deprecatedDefaultDisplayUnit

  const displayUnitObj = useMemo(() => {
    if (!dimension) {
      return
    }
    let defaultDisplayUnitObj
    const section_default_unit_system = getUnitSystem(quantityDef)

    if (isReset && (section_default_unit_system || defaultDisplayUnit)) {
      if (section_default_unit_system && unitSystems[section_default_unit_system]) {
        defaultDisplayUnitObj = new Unit(unitSystems[section_default_unit_system].units[dimension].definition)
      }
      if (defaultDisplayUnit) {
        try {
          defaultDisplayUnitObj = new Unit(defaultDisplayUnit)
        } catch (e) {
          raiseError(`The provided defaultDisplayUnit for ${quantityDef.name} field is not valid.`)
        }
        if (defaultDisplayUnitObj.dimension(false) !== dimension) {
          raiseError(`The provided defaultDisplayUnit for ${quantityDef.name} has incorrect dimensionality for this field.`)
        }
      }
    } else {
      defaultDisplayUnitObj = new Unit(defaultUnit).toSystem(units)
    }
    return defaultDisplayUnitObj
  }, [defaultDisplayUnit, defaultUnit, dimension, isReset, quantityDef, raiseError, unitSystems, units])

  return {displayUnit: displayUnitObj}
}
