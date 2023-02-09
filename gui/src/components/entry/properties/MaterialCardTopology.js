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
import React, { useCallback, useEffect, useState, useMemo } from 'react'
import PropTypes from 'prop-types'
import { isEmpty } from 'lodash'
import { PropertyCard } from './PropertyCard'
import { Tab, Tabs, Typography, Box } from '@material-ui/core'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import { getTopology } from '../../visualization/Structure'
import StructureNGL from '../../visualization/StructureNGL'
// import { DownloadSystemButton } from '../../visualization/StructureBase'
import { FloatableNoReparent } from '../../visualization/Floatable'
import NoData from '../../visualization/NoData'
import TreeView from '@material-ui/lab/TreeView'
import TreeItem from '@material-ui/lab/TreeItem'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ChevronRightIcon from '@material-ui/icons/ChevronRight'
import { QuantityTable, QuantityRow, QuantityCell } from '../../Quantity'

/**
 * For displaying the most descriptive chemical formula that is present in an
 * entry.
*/
const nElementMap = {
  1: 'unary',
  2: 'binary',
  3: 'ternary',
  4: 'quaternary',
  5: 'quinary',
  6: 'sexinary',
  7: 'septenary',
  8: 'octanary',
  9: 'nonary',
  10: 'decimary'
}

/**
 * Displays information about the material and it's structure.
 */
const useMaterialCardStyles = makeStyles((theme) => ({
  root: {
    height: '650px'
  },
  card: {
    height: '100%',
    width: '100%'
  }
}))
const MaterialCardTopology = React.memo(({index, archive}) => {
  const styles = useMaterialCardStyles()
  const [topologyTree, topologyMap] = useMemo(() => getTopology(index, archive), [index, archive])
  const [selected, setSelected] = useState(topologyTree.system_id)
  const [tab, setTab] = useState(0)
  const [structure, setStructure] = useState()
  const [float, setFloat] = useState(false)

  // Used to resolve the path to the structure for the given topology item
  const resolveStructure = useCallback((top, topologyMap) => {
    return top.atoms_ref ? top.atoms_ref.slice(0, -6) : (top.atoms
      ? top.system_id
      : resolveStructure(topologyMap[top.parent_system], topologyMap))
  }, [])

  // When archive is loaded, this effect handles changes in visualizing the
  // selected system
  useEffect(() => {
    if (!archive) return
    const top = topologyMap[selected]
    const structure = resolveStructure(top, topologyMap)
    setStructure(structure)
  }, [archive, selected, topologyMap, resolveStructure])

  // Handle tab change
  const handleTabChange = useCallback((event, value) => {
    setTab(value)
  }, [])

  // Handle float change
  const handleFloat = useCallback(() => {
    setFloat(old => !old)
  }, [])

  return <div className={styles.root}>
    <PropertyCard title="Material" className={styles.card}>
      <FloatableNoReparent float={float} onFloat={handleFloat}>
        <Box display="flex" flexDirection="column" height="100%" width="100%">
          <Box display="flex" flex="1 1 auto" minHeight={0}>
            <Box display="flex" flexDirection="row" width="100%">
              <Box flex="1 1 33%">
                <Topology
                  topologyTree={topologyTree}
                  topologyMap={topologyMap}
                  selected={selected}
                  onSelect={setSelected}
                />
              </Box>
              <Box flex="0 0 8px"/>
              <Box flex="1 1 67%">
                <StructureNGL
                  data={structure}
                  topologyTree={topologyTree}
                  topologyMap={topologyMap}
                  entryId={index.entry_id}
                  selection={selected}
                  onFullscreen={handleFloat}
                />
              </Box>
            </Box>
          </Box>
          <Box flex="0 0 auto" minHeight={0}>
            <MaterialTabs value={tab} onChange={handleTabChange} node={topologyMap[selected]} />
          </Box>
        </Box>
      </FloatableNoReparent>
    </PropertyCard>
  </div>
})

MaterialCardTopology.propTypes = {
  index: PropTypes.object.isRequired,
  archive: PropTypes.object
}
const useTopologyStyles = makeStyles((theme) => ({
  root: {
    height: '100%',
    display: 'flex',
    flexDirection: 'column'
  },
  tree: {
    borderTop: `1px solid ${theme.palette.divider}`,
    overflow: 'auto'
  },
  spacer: {
    flexGrow: 1
  }
}))

/**
 * Displays the topology in an interactive tree.
 */
const Topology = React.memo(({topologyTree, topologyMap, selected, onSelect}) => {
  const styles = useTopologyStyles()

  // Add selection handlers and gather a list of all nodes.
  const nodes = useMemo(() => {
    const nodes = []
    Object.values(topologyMap).forEach(node => {
      node.onClick = () => onSelect(node.system_id)
      nodes.push(node.system_id)
    })
    return nodes
  }, [onSelect, topologyMap])

  const [expanded, setExpanded] = useState(nodes)
  const description = topologyMap?.[selected].description

  // Handle node expansion
  const handleToggle = useCallback((event, nodeIds) => {
    setExpanded(nodeIds)
  }, [])

  return <div className={styles.root}>
    <TreeView
      selected={selected}
      expanded={expanded}
      className={styles.tree}
      onNodeToggle={handleToggle}
      defaultCollapseIcon={<ExpandMoreIcon />}
      defaultExpandIcon={<ChevronRightIcon />}
    >
      <TopologyItem node={topologyTree} level={0} selected={selected}/>
    </TreeView>
    <div className={styles.spacer} />
    {description && <Box padding={1} paddingBottom={0}>
      <Typography variant="body1">Description:</Typography>
      <Typography variant="caption" >{description}</Typography>
    </Box>}
  </div>
})
Topology.propTypes = {
  topologyTree: PropTypes.object,
  topologyMap: PropTypes.object,
  selected: PropTypes.string,
  onSelect: PropTypes.func
}

const useTopologyItemStyles = makeStyles({
  root: (props) => ({
    '&:hover > $content': {
      backgroundColor: props.theme.palette.action.hover
    },
    '&$selected > $content': {
      backgroundColor: props.theme.palette.primary.main,
      color: 'white'
    },
    '&:focus > $content $label, &:hover > $content $label, &$selected > $content $label': {
      backgroundColor: 'transparent'
    }
  }),
  content: (props) => ({
    borderBottom: `1px solid ${props.theme.palette.divider}`,
    fontWeight: props.theme.typography.fontWeightMedium,
    '$expanded > &': {
      fontWeight: props.theme.typography.fontWeightRegular
    },
    paddingLeft: props.theme.spacing(props.level),
    boxSizing: 'border-box'
  }),
  group: (props) => ({marginLeft: 0}),
  // These empty classes are needed for reset
  expanded: (props) => ({}),
  selected: (props) => ({}),
  label: (props) => ({
    fontWeight: 'inherit',
    color: 'inherit'
  }),
  row: (props) => ({
    display: 'flex',
    alignItems: 'center',
    flexDirection: 'row',
    paddingRight: props.theme.spacing(0.5)
  }),
  column: (props) => ({
    display: 'flex',
    alignItems: 'flex-start',
    justifyContent: 'center',
    flexDirection: 'column',
    padding: props.theme.spacing(0.25, 0),
    flexGrow: 1
  }),
  nodeLabelPrimary: (props) => ({
    textTransform: 'uppercase',
    fontSize: 14,
    marginBottom: -props.theme.spacing(0.3)
  }),
  nodeLabelSecondary: (props) => ({
    fontSize: 12
  }),
  iconContainer: (props) => ({
    paddingLeft: props.theme.spacing(0.5)
  })
})

const TopologyItem = React.memo(({node, level, selected}) => {
  // The best way to dynamically manipulate the classes created by makeStyles
  // was to make each style dependent on a set of properties. Creating the
  // styles within the function does not seem to work.
  const theme = useTheme()
  const styleProps = useMemo(() => ({theme, level: level}), [theme, level])
  const classes = useTopologyItemStyles(styleProps)
  const isSelected = selected === node.system_id

  return <TreeItem
    onLabelClick={(event) => { node.onClick && node.onClick(); event.preventDefault() }}
    nodeId={node.system_id}
    label={
      <div className={classes.row}>
        <div className={classes.column}>
          <Typography variant="body1" className={classes.nodeLabelPrimary}>
            {`${node.label}`}
          </Typography>
          <Typography
            variant="body2"
            className={classes.nodeLabelSecondary}
            style={{
              color: isSelected ? 'white' : theme.palette.text.secondary
            }}
          >
            {`${node.structural_type}` + (!isEmpty(node?.child_systems) && node.structural_type === 'group'
              ? ` (${[...new Set(node.child_systems.map(x => x.structural_type))].join(', ')})`
              : ''
            )}
          </Typography>
        </div>
        {/* {node.atoms_ref &&
          <DownloadSystemButton path={node.atoms_ref} color={isSelected ? 'inherit' : 'default'}/>
        } */}
      </div>
    }
    classes={{
      root: classes.root,
      content: classes.content,
      expanded: classes.expanded,
      selected: classes.selected,
      iconContainer: classes.iconContainer,
      group: classes.group,
      label: classes.label
    }}
  >
    {node.child_systems && node.child_systems.map(node => {
      return <TopologyItem key={node.system_id} node={node} level={level + 1} selected={selected}/>
    })}
  </TreeItem>
})
TopologyItem.propTypes = {
  node: PropTypes.object,
  level: PropTypes.number,
  selected: PropTypes.string
}

/**
 * Displays the information that is present for a single node within tabs.
 */
const useMaterialTabsStyles = makeStyles(theme => ({
  noData: {
    height: 116.81 // The height of two QuantityRows
  }
}))
const MaterialTabs = React.memo(({value, onChange, node}) => {
  const styles = useMaterialTabsStyles()
  const nElems = node?.n_elements
  const nElemsLabel = nElementMap[nElems]
  const n_elements = `${nElems}${nElemsLabel ? ` (${nElemsLabel})` : ''}`
  const tabMap = useMemo(() => {
    const cellTab = {...(node?.cell || {})}
    const symmetryTab = {...(node?.symmetry || {})}
    const prototypeTab = {...(node?.prototype || {})}
    const hasTopology = !isEmpty(node)
    const hasCell = !isEmpty(cellTab)
    const hasSymmetry = !isEmpty(symmetryTab)
    const hasPrototype = !isEmpty(prototypeTab)
    return {
      0: {disabled: !hasTopology, label: 'Composition'},
      1: {disabled: !hasCell, label: 'Cell'},
      2: {disabled: !hasSymmetry, label: 'Symmetry'},
      3: {disabled: !hasPrototype, label: 'Prototype'}
    }
  }, [node])

  return <>
    <Tabs
      value={value}
      onChange={onChange}
      indicatorColor="primary"
      textColor="primary"
    >
      <Tab {...tabMap[0]}/>
      <Tab {...tabMap[1]}/>
      <Tab {...tabMap[2]}/>
      <Tab {...tabMap[3]}/>
    </Tabs>
    <MaterialTab value={value} index={0}>
      {!tabMap[0].disabled
      ? <QuantityTable fixed>
        <QuantityRow>
          <QuantityCell value={node?.chemical_formula_hill} quantity="results.material.topology.chemical_formula_hill"/>
          <QuantityCell value={node?.chemical_formula_iupac} quantity="results.material.topology.chemical_formula_iupac"/>
          <QuantityCell value={node?.structural_type} quantity="results.material.topology.structural_type"/>
          <QuantityCell value={node?.material_id} quantity="results.material.topology.material_id"/>
        </QuantityRow>
        <QuantityRow>
          <QuantityCell value={node?.label} quantity="results.material.topology.label"/>
          <QuantityCell value={node?.elements} quantity="results.material.topology.elements"/>
          <QuantityCell
            label="number of elements"
            quantity="results.material.topology.n_elements"
            value={n_elements}
            format={false}
          />
          <QuantityCell value={node?.n_atoms} quantity="results.material.topology.n_atoms"/>
        </QuantityRow>
      </QuantityTable>
      : <NoData className={styles.noData}/>}
    </MaterialTab>
    <MaterialTab value={value} index={1}>
    {!tabMap[1].disabled
      ? <QuantityTable>
        <QuantityRow>
          <QuantityCell value={node?.cell?.a} quantity="results.material.topology.cell.a"/>
          <QuantityCell value={node?.cell?.b} quantity="results.material.topology.cell.b"/>
          <QuantityCell value={node?.cell?.c} quantity="results.material.topology.cell.c"/>
          <QuantityCell value={node?.cell?.alpha} quantity="results.material.topology.cell.alpha"/>
          <QuantityCell value={node?.cell?.beta} quantity="results.material.topology.cell.beta"/>
          <QuantityCell value={node?.cell?.gamma} quantity="results.material.topology.cell.gamma"/>
        </QuantityRow>
        <QuantityRow>
          <QuantityCell value={node?.cell?.volume} quantity="results.material.topology.cell.volume"/>
          <QuantityCell value={node?.cell?.mass_density} quantity="results.material.topology.cell.mass_density"/>
          <QuantityCell value={node?.cell?.atomic_density} quantity="results.material.topology.cell.atomic_density"/>
        </QuantityRow>
      </QuantityTable>
      : <NoData className={styles.noData}/>}
    </MaterialTab>
    <MaterialTab value={value} index={2}>
      {!tabMap[2].disabled
        ? <QuantityTable>
            <QuantityRow>
              <QuantityCell value={node?.symmetry?.crystal_system} quantity="results.material.topology.symmetry.crystal_system"/>
              <QuantityCell value={node?.symmetry?.bravais_lattice} quantity="results.material.topology.symmetry.bravais_lattice"/>
              <QuantityCell value={node?.symmetry?.strukturbericht_designation} quantity="results.material.topology.symmetry.strukturbericht_designation"/>
              <QuantityCell value={node?.symmetry?.space_group_symbol} quantity="results.material.topology.symmetry.space_group_symbol"/>
            </QuantityRow>
            <QuantityRow>
              <QuantityCell value={node?.symmetry?.space_group_number} quantity="results.material.topology.symmetry.space_group_number"/>
              <QuantityCell value={node?.symmetry?.point_group} quantity="results.material.topology.symmetry.point_group"/>
              <QuantityCell value={node?.symmetry?.hall_number} quantity="results.material.topology.symmetry.hall_number"/>
              <QuantityCell value={node?.symmetry?.hall_symbol} quantity="results.material.topology.symmetry.hall_symbol"/>
            </QuantityRow>
          </QuantityTable>
        : <NoData className={styles.noData}/>
      }
    </MaterialTab>
    <MaterialTab value={value} index={3}>
    {!tabMap[3].disabled
      ? <QuantityTable>
          <QuantityRow>
            <QuantityCell value={node?.prototype?.name} quantity="results.material.topology.prototype.name"/>
            <QuantityCell value={node?.prototype?.aflow_id} quantity="results.material.topology.prototype.aflow_id"/>
            <QuantityCell value={node?.prototype?.formula} quantity="results.material.topology.prototype.formula"/>
          </QuantityRow>
        </QuantityTable>
      : <NoData className={styles.noData}/>}
    </MaterialTab>
  </>
})
MaterialTabs.propTypes = {
  node: PropTypes.object,
  value: PropTypes.any,
  onChange: PropTypes.func
}

const MaterialTab = React.memo(({value, index, children}) => {
  return <div
    role="tabpanel"
    hidden={value !== index}
    id={`full-width-tabpanel-${index}`}
    aria-labelledby={`full-width-tab-${index}`}
  >
    {children}
  </div>
})
MaterialTab.propTypes = {
  index: PropTypes.number.isRequired,
  value: PropTypes.number.isRequired,
  children: PropTypes.node
}

export default MaterialCardTopology
