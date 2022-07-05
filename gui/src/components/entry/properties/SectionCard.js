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
import React, {useCallback, useMemo, useRef, useState} from 'react'
import PropTypes from 'prop-types'
import {PropertyCard} from './PropertyCard'
import SectionEditor from '../../archive/SectionEditor'
import { useEntryPageContext } from '../EntryPageContext'
import {Box, IconButton, Typography, makeStyles} from '@material-ui/core'
import CodeIcon from '@material-ui/icons/Code'
import MoreIcon from '@material-ui/icons/MoreVert'
import { ArchiveButton } from '../../nav/Routes'
import {SectionPlots} from '../../archive/ArchiveBrowser'
import Quantity, {QuantityCell, QuantityRow, QuantityTable} from '../../Quantity'
import {Editor} from '@tinymce/tinymce-react'
import { pluralize } from '../../../utils'
import {Matrix} from '../../archive/visualizations'

const useStyles = makeStyles(theme => ({
  subSection: {
    width: '100%',
    display: 'flex',
    flexDirection: 'column',
    justifyContent: 'center'
  },
  arrayAction: {
    padding: 0
  },
  subSectionAction: {
  },
  row: {
    display: 'flex',
    justifyContent: 'space-between'
  }
}))

const PropertyPreview = React.memo(({quantityDef, section}) => {
  const classes = useStyles()
  const {entryId, uploadId} = useEntryPageContext()
  const maxPreviewLength = 5
  if (!quantityDef.type) {
    return ''
  }

  const shape = quantityDef.shape || []
  if (shape.length === 0) {
    if (quantityDef.m_annotations?.eln?.[0]?.component === 'RichTextEditQuantity') {
      return <QuantityCell >
        <Quantity data={section} quantity={quantityDef}>
          <Editor
            init={{
              height: '600px',
              width: '100%',
              plugins: [
                'advlist autolink lists link image charmap print preview anchor',
                'searchreplace visualblocks code',
                'insertdatetime media table paste code help wordcount'
              ],
              skin: 'nomadPreview',
              resize: false,
              toolbar: false,
              menubar: false,
              branding: false,
              statusbar: false
            }}
            initialValue={section[quantityDef.name]}
            disabled={true}
          />
        </Quantity>
      </QuantityCell>
    }
    return <QuantityCell >
      <Quantity data={section} quantity={quantityDef} />
    </QuantityCell>
  } else if (shape.length === 1) {
    const values = section[quantityDef.name] || []
    return <QuantityCell>
      <Quantity data={section} quantity={quantityDef} >
        <Box width='100%'>
          {values.slice(0, maxPreviewLength).map((value, index) =>
            <Quantity key={index} data={{[quantityDef.name]: value}} quantity={quantityDef} noLabel/>
          )}
          {values.length >= maxPreviewLength &&
            <Box className={classes.row}>
              <Typography noWrap variant={'body1'}>
                {`and ${values.length - maxPreviewLength} more ${pluralize('item', values.length - maxPreviewLength, false)}`}
              </Typography>
              <ArchiveButton
                className={classes.arrayAction}
                component={IconButton}
                entryId={entryId} uploadId={uploadId}
                path={`data/${quantityDef.name}`}
              >
                <MoreIcon/>
              </ArchiveButton>
            </Box>
          }
        </Box>
      </Quantity>
    </QuantityCell>
  } else {
    return <QuantityCell>
      <Quantity data={section} quantity={quantityDef} >
        <Box width={'100%'} textAlign='center'>
          <Matrix
            values={section[quantityDef.name]}
            shape={shape}
            invert={false}
            type={quantityDef.type.type_data}
            key={`matrix:${quantityDef.name}`}
          />
          <Typography noWrap variant="caption">
            ({quantityDef.shape.map((dimension, index) => <span key={index}>
              {index > 0 && <span>&nbsp;&times;&nbsp;</span>}{String(dimension)}
            </span>)}&nbsp;)
          </Typography>
        </Box>
      </Quantity>
    </QuantityCell>
  }
})
PropertyPreview.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired
}

const SectionPreview = React.memo(({sectionDef, section}) => {
  const rootRef = useRef()

  const filterHiddenProperties = useCallback((property) => {
    const hiddenPropertyNames = sectionDef?.m_annotations?.eln?.[0]?.hide || []
    return !hiddenPropertyNames.includes(property.name)
  }, [sectionDef])

  const properties = useMemo(() => {
    return sectionDef._allProperties
      .filter(filterHiddenProperties)
  }, [filterHiddenProperties, sectionDef._allProperties])

  return (
    <div ref={rootRef}>
      {<QuantityTable data={section}>
        {properties.map(property => (
          <QuantityRow key={property.name} >
            <PropertyPreview quantityDef={property} section={section}/>
          </QuantityRow>
        ))}
      </QuantityTable>}
    </div>
  )
})
SectionPreview.propTypes = {
  sectionDef: PropTypes.object.isRequired,
  section: PropTypes.object
}

const SectionCard = React.memo(({archivePath, sectionDef, section, readOnly, ...props}) => {
  const {entryId, uploadId} = useEntryPageContext()
  const [showJson, setShowJson] = useState(false)

  const actions = <React.Fragment>
    <IconButton onClick={() => setShowJson(value => !value)}>
      <CodeIcon />
    </IconButton>
    <ArchiveButton
      component={IconButton}
      entryId={entryId} uploadId={uploadId}
      path={archivePath}
    >
      <MoreIcon/>
    </ArchiveButton>
  </React.Fragment>

  if (!sectionDef) {
    console.error('SectionCard: section definition is not available')
    return ''
  }

  return <PropertyCard title={sectionDef.name} action={actions}>
    {sectionDef._qualifiedName === 'nomad.parsing.tabular.Table' && (
      <Box margin={2}>
        <Typography>Data from a table with {section.row_refs?.length || '...'} rows.</Typography>
      </Box>
    )}
    <Box margin={2}>
      {readOnly
        ? (
          <SectionPreview
            sectionDef={sectionDef}
            section={section}
            archivePath={archivePath}
            {...props}
          />
        ) : (
          <SectionEditor
            sectionDef={sectionDef}
            section={section}
            showJson={showJson}
            {...props}
          />
        )}
      {sectionDef.m_annotations?.plot && <SectionPlots sectionDef={sectionDef} section={section}/>}
    </Box>
  </PropertyCard>
})

SectionCard.propTypes = {
  archivePath: PropTypes.string.isRequired,
  sectionDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  readOnly: PropTypes.bool
}

export default SectionCard
