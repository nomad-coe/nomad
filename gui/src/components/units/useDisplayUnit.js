import {useErrors} from "../errors"
import {useMemo} from "react"
import {Unit} from "./Unit"
import {useUnitContext} from "./UnitContext"
import {getFieldProps} from "../editQuantity/StringEditQuantity"

export function useDisplayUnit(quantityDef) {
  const {units} = useUnitContext()
  const {raiseError} = useErrors()
  const defaultUnit = useMemo(() => quantityDef.unit && new Unit(quantityDef.unit), [quantityDef])
  const dimension = defaultUnit && defaultUnit.dimension(false)
  const {defaultDisplayUnit: deprecatedDefaultDisplayUnit} = getFieldProps(quantityDef)
  const defaultDisplayUnit = quantityDef?.m_annotations?.display?.[0]?.unit || deprecatedDefaultDisplayUnit

  const displayUnitObj = useMemo(() => {
    if (!dimension) return
    let defaultDisplayUnitObj

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
    // Use the global unit system defined in the schema
    } else {
      defaultDisplayUnitObj = new Unit(defaultUnit).toSystem(units)
    }

    return defaultDisplayUnitObj
  }, [defaultDisplayUnit, defaultUnit, dimension, quantityDef, raiseError, units])

  return displayUnitObj
}
