import React, { useState, useCallback, useMemo } from 'react'
import PropTypes from 'prop-types'
import {
  FormControlLabel,
  FormControl,
  RadioGroup,
  Radio
} from '@material-ui/core'
import { makeStyles } from '@material-ui/core/styles'
import Structure from '../visualization/Structure'
import BrillouinZone from '../visualization/BrillouinZone'
import BandStructure from '../visualization/BandStructure'
import Spectra from '../visualization/Spectra'
import DOS from '../visualization/DOS'
import { Quantity } from '../units/Quantity'
import { useUnitContext } from '../units/UnitContext'
import { electronicRange } from '../../config'
import EnergyVolumeCurve from '../visualization/EnergyVolumeCurve'

const useOverviewAtomsStyles = makeStyles({
  root: {
    width: '30rem',
    height: '20rem',
    margin: 'auto'
  }
})
export const OverviewAtoms = React.memo(({def, section}) => {
  const style = useOverviewAtomsStyles()
  const system = useMemo(() => ({
    'species': section.species,
    'cell': section.lattice_vectors
      ? new Quantity(section.lattice_vectors, 'meter').to('angstrom').value()
      : undefined,
    'positions': new Quantity(section.positions, 'meter').to('angstrom').value(),
    'pbc': section.periodic
  }), [section])

  return <Structure
    className={style.root}
    data={system}
  ></Structure>
})
OverviewAtoms.propTypes = ({
  def: PropTypes.object,
  section: PropTypes.object,
  units: PropTypes.object
})

const useOverviewDOSStyles = makeStyles({
  root: {
    width: '20rem',
    height: '40rem',
    margin: 'auto'
  }
})
export const OverviewDOSElectronic = React.memo(({def, section, units}) => {
  const style = useOverviewDOSStyles()
  const data = useMemo(() => ([{
    version: 'old',
    energies: section.energies,
    densities: section.total.map(dos => dos.value),
    energy_highest_occupied: section.band_gap
      ? Math.max(...section.band_gap.map(x => x.energy_highest_occupied))
      : undefined
  }]), [section])

  return <DOS
    className={style.root}
    data={data}
    units={units}
    type="electronic"
  />
})
OverviewDOSElectronic.propTypes = ({
  def: PropTypes.object,
  section: PropTypes.object,
  units: PropTypes.object
})

export const OverviewDOSPhonon = React.memo(({def, section, units}) => {
  const style = useOverviewDOSStyles()
  const data = useMemo(() => ([{
    energies: section.energies,
    densities: section.total.map(dos => dos.value)
  }]), [section])

  return <DOS
    className={style.root}
    data={data}
    units={units}
    type="vibrational"
  />
})
OverviewDOSPhonon.propTypes = ({
  def: PropTypes.object,
  section: PropTypes.object,
  units: PropTypes.object
})

const useOverviewBandstructurePhononStyles = makeStyles({
  root: {
    width: '30rem',
    height: '30rem',
    margin: 'auto'
  }
})
export const OverviewBandstructurePhonon = React.memo(({def, section, units}) => {
  const style = useOverviewBandstructurePhononStyles()
  const data = useMemo(() => ([{
    segment: section.segment,
    reciprocal_cell: section.reciprocal_cell
  }]), [section])

  return <BandStructure
    className={style.root}
    data={data}
    units={units}
    type='vibrational'
  />
})
OverviewBandstructurePhonon.propTypes = ({
  def: PropTypes.object,
  section: PropTypes.object,
  units: PropTypes.object
})

const useOverviewBandstructureElectronicStyles = makeStyles({
  root: {
    width: '30rem',
    height: '30rem',
    margin: 'auto'
  },
  radio: {
    display: 'flex',
    justifyContent: 'center'
  }
})
export const OverviewBandstructureElectronic = React.memo(({def, section, units}) => {
  const [mode, setMode] = useState('bs')
  const toggleMode = useCallback((event) => {
    setMode(event.target.value)
  }, [setMode])

  const style = useOverviewBandstructureElectronicStyles()
  const data = useMemo(() => ([{
    segment: section.segment,
    reciprocal_cell: section.reciprocal_cell,
    energy_highest_occupied: section.band_gap
      ? Math.max(...section.band_gap.map(x => x.energy_highest_occupied))
      : undefined
  }]), [section])
  const layout = useMemo(() => ({
    yaxis: {
      autorange: false,
      range: new Quantity(electronicRange, 'electron_volt').toSystem(units).value()
    }
  }), [units])

  const bs = <BandStructure
    className={style.root}
    data={data}
    layout={layout}
    units={units}
  />

  return data[0].reciprocal_cell
    ? <>
      {mode === 'bs'
        ? bs
        : <BrillouinZone
          className={style.root}
          data={section}
        />
      }
      <FormControl component="fieldset" className={style.radio}>
        <RadioGroup
          row
          aria-label="position"
          name="position"
          defaultValue="bs"
          onChange={toggleMode}
          className={style.radio}
        >
          <FormControlLabel
            value="bs"
            control={<Radio color="primary" />}
            label="Band structure"
            labelPlacement="end"
          />
          <FormControlLabel
            value="bz"
            control={<Radio color="primary" />}
            label="Brillouin zone"
            labelPlacement="end"
          />
        </RadioGroup>
      </FormControl>
    </>
    : bs
})
OverviewBandstructureElectronic.propTypes = ({
  def: PropTypes.object,
  section: PropTypes.object,
  units: PropTypes.object
})

const useOverviewSpectraStyles = makeStyles({
  root: {
    width: '30rem',
    height: '15rem',
    margin: 'auto'
  }
})
export const OverviewSpectra = React.memo(({def, section, units}) => {
  const style = useOverviewSpectraStyles()

  return <Spectra
    className={style.root}
    data={[section]}
    layout={{yaxis: {autorange: true}}}
    units={units}
  />
})
OverviewSpectra.propTypes = ({
  def: PropTypes.object,
  section: PropTypes.object,
  units: PropTypes.object
})

const useOverviewEquationOfStateStyles = makeStyles({
  root: {
    width: '30rem',
    height: '25rem',
    margin: 'auto'
  }
})
export const OverviewEquationOfState = React.memo(({def, section, units}) => {
  const style = useOverviewEquationOfStateStyles()
  const data = [{
    volumes: section.volumes,
    energies: section.energies,
    name: 'raw'
  }]
  section.eos_fit.forEach(fit => {
    data.push({
      volumes: section.volumes,
      energies: fit.fitted_energies,
      name: fit.function_name
    })
  })

  return <EnergyVolumeCurve
    className={style.root}
    data={{data: data}}
    units={units}
  />
})
OverviewEquationOfState.propTypes = ({
  def: PropTypes.object,
  section: PropTypes.object,
  units: PropTypes.object
})

export const Overview = React.memo((props) => {
  const {def} = props
  const {units} = useUnitContext()
  const path = window.location.href.split('/').pop().split(':')[0]

  if (def.name === 'BandStructure' && path === 'band_structure_electronic') {
    return <OverviewBandstructureElectronic {...props} units={units}/>
  } else if (def.name === 'BandStructure' && path === 'band_structure_phonon') {
    return <OverviewBandstructurePhonon {...props} units={units}/>
  } else if (def.name === 'Atoms' && path === 'atoms') {
    return <OverviewAtoms {...props} units={units} />
  } else if (def.name === 'Dos' && path === 'dos_electronic') {
    return <OverviewDOSElectronic {...props} units={units} />
  } else if (def.name === 'Dos' && path === 'dos_phonon') {
    return <OverviewDOSPhonon {...props} units={units}/>
  } else if (def.name === 'Spectrum') {
    return <OverviewSpectra {...props} units={units}/>
  } else if (def.name === 'EquationOfState') {
    return <OverviewEquationOfState {...props} units={units}/>
  }
  return null
})
Overview.propTypes = ({
  def: PropTypes.object,
  section: PropTypes.object
})
