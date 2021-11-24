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

import React from 'react'
import PropTypes from 'prop-types'
import {
  Divider,
  Typography,
  Link,
  AccordionActions,
  Button,
  Grid,
  makeStyles
} from '@material-ui/core'
import MuiAccordion from '@material-ui/core/Accordion'
import MuiAccordionSummary from '@material-ui/core/AccordionSummary'
import MuiAccordionDetails from '@material-ui/core/AccordionDetails'
import { styled } from '@material-ui/core/styles'
import ArrowForwardIosSharpIcon from '@material-ui/icons/ArrowForwardIosSharp'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import Markdown from '../Markdown'
import { StringParam, useQueryParam } from 'use-query-params'
import AccessIcon from '../../images/AIT_ico_bd_link_external_big.svg'
import WatchIcon from '../../images/AIT_ico_bd_youtube.svg'
import PdfIcon from '../../images/AIT_ico_bd_link_pdf.svg'
import DoiIcon from '../../images/AIT_ico_bd_link_doi.svg'

const NotebookButton = React.memo(function NotebookButton({...props}) {
  // TODO this button has to be implemented and send users to the respective AItoolkit
  // notebook on north.
  return <Button {...props} />
})
NotebookButton.propTypes = {
  metadata: PropTypes.Object,
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
}

const useStyles = makeStyles(theme => ({
  tutorialTitleGrid: {
    marginRight: '40px'
  },
  tutorialTitleText: {
    fontWeight: theme.typography.fontWeightMedium,
    fontSize: '28px',
    color: '#2A3C67',
    lineHeight: '30px'
  },
  fieldText: {
    color: '#2A3C67'
  },
  linkAuthors: {
    color: '#2A3C67',
    cursor: 'pointer',
    lineHeight: '20px',
    fontSize: '16px'
  },
  tutorialDescriptionGrid: {
    marginLeft: '50px'
  },
  tutorialDescriptionText: {
    color: '#2A3C67',
    fontSize: '18px'
  },
  keywordsGrid: {
    marginLeft: '80px'
  },
  linkKeywords: {
    border: '1.5px solid #00DFE0',
    lineHeight: '35px',
    color: '#2A3C67',
    cursor: 'pointer',
    fontStyle: 'normal',
    fontSize: '16px'
  },
  tutorialActions: {
    marginLeft: '50px'
  },
  tutorialResources: {
    marginTop: '-17px',
    marginLeft: '-6px'
  }
}))

const Accordion = styled((props) => (
  <MuiAccordion {...props} />
))(({ theme }) => ({
  backgroundColor: theme.palette.background.default,
  borderBottom: '13px solid #7FEFEF'
}))

const AccordionSummary = styled((props) => (
  <MuiAccordionSummary
    expandIcon={<ArrowForwardIosSharpIcon sx={{ fontSize: '0.9rem' }} />}
    {...props}
  />
))(({ theme }) => ({
  flexDirection: 'row-reverse',
  '& .MuiAccordionSummary-expandIconWrapper.Mui-expanded': {
    transform: 'rotate(90deg)'
  },
  '& .MuiAccordionSummary-content': {
    marginLeft: theme.spacing(1),
    marginTop: '20px',
    marginBottom: '20px'
  }
}))

const AccordionDetails = styled(MuiAccordionDetails)(({ theme }) => ({
  padding: theme.spacing(2)
}))

function AccordionsList(props) {
  const styles = useStyles()
  const [expanded, setExpanded] = useQueryParam('expanded', StringParam)
  return (
    props.tutorials_list.map(tutorial => (
      props.filter(tutorial) &&
      <div key={tutorial.title} >
        <Accordion
          key={tutorial.key}
          expanded={expanded === tutorial.key}
          onChange={() => setExpanded(expanded === tutorial.key ? null : tutorial.key)}
          elevation={0}
        >
          <AccordionSummary expandIcon={<ExpandMoreIcon />}>
            <Grid container spacing={1} >
              <Grid item xs={7} className={styles.tutorialTitleGrid} >
                <Typography className={styles.tutorialTitleText}>
                  {tutorial.title}
                </Typography>
              </Grid>
              <Grid item xs={4}>
                <Typography className={styles.fieldText}>
                  {<b>Authors: </b> }
                  {tutorial.authors
                    .map(name => {
                      const label = name.split(',').reverse().join(' ')
                      return <Link
                        className={styles.linkAuthors}
                        key={name}
                        onClick={() => {
                          props.setQueryParameters({
                            ...props.emptyQuery,
                            author: props.queryParameters.author === name ? null : name
                          })
                          props.setAuthor(name)
                        }}
                      >
                        {label}
                      </Link>
                    }).reduce((prev, curr) => [prev, ' | ', curr])
                  }
                </Typography>
              </Grid>
            </Grid>
          </AccordionSummary>
          <AccordionDetails >
            <Grid container spacing={4}>
              <Grid item xs={6} className={styles.tutorialDescriptionGrid}>
                <Markdown className={styles.tutorialDescriptionText}>
                  {tutorial.description}
                </Markdown>
              </Grid>
              <Grid item xs={4} className={styles.keywordsGrid}>
                <Typography className={styles.fieldText}>
                  <b>AI methods</b>:
                </Typography>
                <Typography>
                  {tutorial.labels.ai_methods
                    .map(method => (
                      <Link
                        className={styles.linkKeywords}
                        key={method}
                        onClick={() => {
                          props.setQueryParameters({
                            ...props.emptyQuery,
                            method: props.queryParameters.method === method ? null : method
                          })
                          props.setMethod(method)
                        }}
                      >
                        {method}
                      </Link>
                    )).reduce((prev, curr) => [prev, '    ', curr])
                  }
                </Typography>
                <Typography className={styles.fieldText}>
                  <b>System</b>:
                </Typography>
                <Typography>
                  {tutorial.labels.application_system
                    .map(system => (
                      <Link
                        className={styles.linkKeywords}
                        key={system}
                        onClick={() => {
                          props.setQueryParameters({
                            ...props.emptyQuery,
                            system: props.queryParameters.system === system ? null : system
                          })
                          props.setSystem(system)
                        }}
                      >
                        {system}
                      </Link>
                    )).reduce((prev, curr) => [prev, '    ', curr])
                  }
                </Typography>
              </Grid>
            </Grid>
          </AccordionDetails>
          <AccordionActions>
            <Grid container spacing={4}>
              <Grid item xs={7} className={styles.tutorialActions}>
                <Grid container spacing={0}>
                  <Grid item xs={5}>
                    <NotebookButton
                      metadata={tutorial}
                      startIcon={<img alt='Access icon' src={AccessIcon}></img>}
                    >
                      <Typography className={styles.fieldText} >
                        <b>Access tutorial</b>
                      </Typography>
                    </NotebookButton>
                  </Grid>
                  <Grid item xs={5} >
                    <div>
                      { tutorial.link_video && <Button
                        width='10px'
                        color='#2A3C67'
                        href={tutorial.link_video}
                        target="tutorial"
                        startIcon={<img alt='Watch icon' src={WatchIcon}></img>}
                      >
                        <Typography className={styles.fieldText} >
                          <b>Watch video</b>
                        </Typography>
                      </Button>}
                    </div>
                  </Grid>
                </Grid>
              </Grid>
              {tutorial.link_paper &&
              <Grid item xs={4} className={styles.tutorialResources}>
                <Grid item xs={12}>
                  <Typography className={styles.fieldText}>
                    <b>Additional Resources</b>:
                  </Typography>
                </Grid>
                <Grid container spacing={0}>
                  <Grid item xs={2} >
                    <div>
                      {tutorial.link_paper && <Button
                        color='#2A3C67'
                        href={tutorial.link_paper}
                        target="tutorial"
                        startIcon={<img alt='DOI icon' src={DoiIcon}></img>}>
                      </Button>}
                    </div>
                  </Grid>
                  <Grid item xs={2}>
                    <div>
                      {tutorial.link_paper && <Button
                        color='#2A3C67'
                        href={tutorial.link_paper}
                        target="tutorial"
                        startIcon={<img alt='PDF icon' src={PdfIcon}/>}>
                      </Button>}
                    </div>
                  </Grid>
                </Grid>
              </Grid>
              }
            </Grid>
          </AccordionActions>
          <Divider />
        </Accordion>
      </div>
    )))
}

export default AccordionsList
