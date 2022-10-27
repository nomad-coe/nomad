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
import React, { useEffect, useMemo, useState } from 'react'
import PropTypes from 'prop-types'
import { Typography, makeStyles, Box, Grid, Divider } from '@material-ui/core'
import { resolveNomadUrlNoThrow } from '../../utils'
import Quantity from '../Quantity'
import ElectronicPropertiesCard from '../entry/properties/ElectronicPropertiesCard'
import OptoelectronicPropertiesCard from '../entry/properties/OptoelectronicPropertiesCard'
import MaterialCard from '../entry/properties/MaterialCard'
import MaterialCardTopology from '../entry/properties/MaterialCardTopology'
import NexusCard from './properties/NexusCard'
import VibrationalPropertiesCard from '../entry/properties/VibrationalPropertiesCard'
import MechanicalPropertiesCard from '../entry/properties/MechanicalPropertiesCard'
import ThermodynamicPropertiesCard from '../entry/properties/ThermodynamicPropertiesCard'
import StructuralProperties from '../entry/properties/StructuralPropertiesCard'
import DynamicalProperties from '../entry/properties/DynamicalPropertiesCard'
import GeometryOptimizationCard from '../entry/properties/GeometryOptimizationCard'
import SpectroscopyCard from './properties/SpectroscopyCard'
import { MethodMetadata } from './EntryDetails'
import Page from '../Page'
import { SourceApiCall, SourceApiDialogButton, SourceDialogDivider } from '../buttons/SourceDialogButton'
import { useEntryPageContext } from './EntryPageContext'
import SectionCard from './properties/SectionCard'
import { useMetainfoDef, traverse } from '../archive/metainfo'
import {
  ArchiveSaveButton, ArchiveDeleteButton, ArchiveReloadButton, ArchiveReUploadButton
} from '../archive/ArchiveBrowser'
import { useErrors } from '../errors'
import DefinitionsCard from './properties/DefinitionsCard'
import { ErrorHandler } from '../ErrorHandler'

function MetadataSection({title, children}) {
  return <Box marginTop={2} marginBottom={2}>
    {title && <Typography component="div">
      <Box fontSize="h6.fontSize" marginBottom={1}>
        {title}
      </Box>
    </Typography>}
    {children}
  </Box>
}

MetadataSection.propTypes = {
  title: PropTypes.string,
  children: PropTypes.any
}

const useStyles = makeStyles(theme => ({
  root: {
    marginBottom: theme.spacing(4)
  },
  leftColumn: {
    maxWidth: '32%',
    flexBasis: '32%',
    flexGrow: 0,
    paddingRight: theme.spacing(3)
  },
  rightColumn: {
    maxWidth: '67.99%',
    flexBasis: '67.99%',
    flexGrow: 0,
    '& > div': {
      marginBottom: theme.spacing(2)
    }
  },
  editActions: {
    marginBottom: `${theme.spacing(1)}px !important`
  },
  divider: {
    marginTop: theme.spacing(1),
    marginBottom: theme.spacing(1)
  }
}))

const overviewArchiveFilter = Object.freeze({
  // Optimization: for the overview page, we only need to fetch some of the archive data
  metadata: '*',
  data: '*',
  definitions: '*',
  results: {
    material: '*',
    method: '*',
    properties: {
      structures: '*',
      structural: '*',
      dynamical: '*',
      electronic: 'include-resolved',
      mechanical: 'include-resolved',
      spectroscopy: 'include-resolved',
      vibrational: 'include-resolved',
      thermodynamic: 'include-resolved',
      // For geometry optimizations we require only the energies.
      // Trajectory, optimized structure, etc. are unnecessary.
      geometry_optimization: {
        energies: 'include-resolved'
      }
    }
  }
})

/**
 * Shows an informative overview about the selected entry.
 */
const OverviewView = React.memo((props) => {
  const {
    url,
    metadata,
    metadataApiData,
    exists,
    editable,
    archive,
    archiveApiData} = useEntryPageContext(overviewArchiveFilter, false)
  const index = metadata
  const [sections, setSections] = useState([])
  const {raiseError} = useErrors()
  const dataMetainfoDefUrl = resolveNomadUrlNoThrow(archive?.data?.m_def, url)
  const dataMetainfoDef = useMetainfoDef(dataMetainfoDefUrl)

  const properties = useMemo(() => {
    return new Set(index?.results
      ? index.results.properties.available_properties
      : []
    )
  }, [index])

  useEffect(() => {
    if (!dataMetainfoDef || !archive?.data) {
      return
    }
    const getSections = async () => {
      const sections = []
      traverse(archive.data, dataMetainfoDef, 'data', (section, sectionDef, path) => {
        if (path === 'data' || sectionDef.m_annotations?.eln?.[0]?.overview) {
          sections.push({
            archivePath: path,
            sectionDef: sectionDef,
            section: section
          })
        }
      })
      return sections
    }
    getSections().then(setSections).catch(raiseError)
  }, [archive, dataMetainfoDef, setSections, raiseError])

  const classes = useStyles()

  if (!exists) {
    return <Page>
      <Typography>
        This entry does not exist.
      </Typography>
    </Page>
  }

  if (!index) {
    return null
  }

  return <Page limitedWidth>
    <Grid container spacing={0} className={classes.root}>
      <Grid item xs={4} className={classes.leftColumn}>
        <MetadataSection title="Metadata">
          <MethodMetadata data={index} />
        </MetadataSection>
        <Divider className={classes.divider} />
        <MetadataSection>
          <Quantity flex>
            <Quantity quantity='comment' data={index} />
            <Quantity quantity='references' data={index}/>
            <Quantity quantity='authors' data={index}/>
            <Quantity quantity="datasets" data={index}/>
          </Quantity>
        </MetadataSection>
        <Divider className={classes.divider}/>
        <MetadataSection>
          <Quantity column style={{maxWidth: 350}}>
            <Quantity quantity="mainfile" data={index}/>
            <Quantity quantity="entry_id" data={index}/>
            <Quantity quantity="results.material.material_id" data={index}/>
            <Quantity quantity="upload_id" data={index}/>
            <Quantity quantity="upload_create_time" data={index}/>
            <Quantity quantity="raw_id" data={index}/>
            <Quantity quantity="external_id" data={index}/>
            <Quantity quantity="last_processing_time" data={index}/>
            <Quantity quantity="last_processing_version" data={index}/>
          </Quantity>
        </MetadataSection>
        <SourceApiDialogButton label="API" maxWidth="lg" fullWidth buttonProps={{variant: 'contained', size: 'small'}}>
          {metadataApiData && <SourceApiCall
            {...metadataApiData}
            description="The basic metadata shown on this page is retrieved from the *entry metadata* API."
          />}
          <SourceDialogDivider />
          {archiveApiData && <SourceApiCall
            {...archiveApiData}
            description="The detailed property information is retrieved from the *entry archive* API. Only a specific parts of the archive are *required*."
          />}
        </SourceApiDialogButton>
      </Grid>

      <Grid item xs={8} className={classes.rightColumn}>
        {editable && (
          <Box textAlign="right" className={classes.editActions} display={'flex'} justifyContent={'flex-end'}>
            <ArchiveReloadButton />
            <ArchiveSaveButton />
            <ArchiveReUploadButton />
            <ArchiveDeleteButton />
          </Box>
        )}
        {sections
          .map((section, index) => (
            <ErrorHandler key={index} message="Could not render section card.">
              <SectionCard
                {...section}
                archivePath={section.archivePath.replace(/\./g, '/')}
                readOnly={!editable}
              />
            </ErrorHandler>
          ))
        }
        <ErrorHandler message="Could not render definitions card.">
          <DefinitionsCard index={index} archive={archive}/>
        </ErrorHandler>
        <ErrorHandler message="Could not render Nexus card.">
          <NexusCard index={index}/>
        </ErrorHandler>
        <ErrorHandler message="Could not render materials card.">
          {index?.results?.material?.topology
            ? <MaterialCardTopology index={index} archive={archive} properties={properties}/>
            : <MaterialCard index={index} archive={archive} properties={properties}/>}
        </ErrorHandler>
        <ErrorHandler message="Could not render electronic properties.">
          <ElectronicPropertiesCard index={index} archive={archive} properties={properties}/>
        </ErrorHandler>
        <ErrorHandler message="Could not render opto-electronic properties.">
          <OptoelectronicPropertiesCard index={index} archive={archive} properties={properties}/>
        </ErrorHandler>
        <ErrorHandler message="Could not render vibrational properties.">
          <VibrationalPropertiesCard index={index} archive={archive} properties={properties}/>
        </ErrorHandler>
        <ErrorHandler message="Could not render mechanical properties.">
          <MechanicalPropertiesCard index={index} archive={archive} properties={properties}/>
        </ErrorHandler>
        <ErrorHandler message="Could not render thermodynamic properties.">
          <ThermodynamicPropertiesCard index={index} archive={archive} properties={properties}/>
        </ErrorHandler>
        <ErrorHandler message="Could not render structural properties.">
          <StructuralProperties index={index} archive={archive} properties={properties}/>
        </ErrorHandler>
        <ErrorHandler message="Could not render dynamical properties.">
          <DynamicalProperties index={index} archive={archive} properties={properties}/>
        </ErrorHandler>
        <ErrorHandler message="Could not render geometry optimization.">
          <GeometryOptimizationCard index={index} archive={archive} properties={properties}/>
        </ErrorHandler>
        <ErrorHandler message="Could not render spectroscopy.">
          <SpectroscopyCard index={index} archive={archive} properties={properties}/>
        </ErrorHandler>
      </Grid>
    </Grid>
  </Page>
})

OverviewView.propTypes = {
}

OverviewView.whyDidYouRender = true

export default OverviewView
