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
import React, { useMemo } from 'react'
import { Box, Divider, Typography, makeStyles, Link, AccordionActions, Button, Grid, TextField } from '@material-ui/core'
import MUIAccordion from '@material-ui/core/Accordion'
import MUIAccordionSummary from '@material-ui/core/AccordionSummary'
import MUIAccordionDetails from '@material-ui/core/AccordionDetails'

import { withStyles } from '@material-ui/core/styles'
import tutorials from '../../toolkitMetadata'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import Markdown from '../Markdown'
import { StringParam, useQueryParams, useQueryParam } from 'use-query-params'
import Autocomplete from '@material-ui/lab/Autocomplete'
import TutorialsIcon from './assets/AIT_ico_bp_tutorial.svg'
import AccessIcon from './assets/AIT_ico_bd_link_external_big.svg'
import WatchIcon from './assets/AIT_ico_bd_youtube.svg'
import PdfIcon from './assets/AIT_ico_bd_link_pdf.svg'
import DoiIcon from './assets/AIT_ico_bd_link_doi.svg'

const useStyles = makeStyles(theme => ({

  root: {
    margin: theme.spacing(3),
    width: '100%',
    marginLeft: 'auto',
    marginRight: 'auto',
    maxWidth: '1024px'
  },
  sectionIcon: {
    marginTop: theme.spacing(3)
  },
  sectionTitle: {
    marginBottom: theme.spacing(1),
    marginLeft: theme.spacing(2),
    marginTop: '105px'
  },
  title: {
    fontWeight: 'bold',
    color: '#2A3C67',
    fontSize: 30,
    marginLeft: '-10px',
    fontFamily: 'TitilliumBold'

  },
  deck: {
    color: '#2A3C67',
    fontSize: 15,
    marginTop: '20px',
    lineHeight: '24px',
    marginLeft: '-10px',
    fontFamily: 'TitilliumRegular'
  },
  icon: {
    height: '350px',
    marginTop: '-20px',
    marginLeft: '100px'
  },
  filter: {
    fontWeight: 'bold',
    color: '#2A3C67',
    fontSize: 15,
    marginTop: '60px',
    marginLeft: '0px',
    fontFamily: 'TitilliumBold'
  },
  autocomplete: {
    height: 'auto',
    color: '#2A3C67',
    border: '3px solid rgba(127, 239, 239, 1)',
    borderRadius: '10px 10px 10px 10px',
    marginTop: '10px',
    marginLeft: '0px'
  },
  tutorialsList: {
    marginTop: '50px'
  },
  tutorialTitleGrid: {
    marginRight: '40px'
  },
  tutorialTitleText: {
    fontSize: 22,
    lineHeight: '30px',
    color: '#2A3C67',
    fontFamily: 'TitilliumBold'
  },
  authorsGrid: {
    marginLeft: '150px',
    marginRight: '30px'
  },
  fieldText: {
    color: '#2A3C67',
    fontWeight: 'TitilliumBold'
  },
  linkAuthors: {
    color: '#2A3C67',
    cursor: 'pointer',
    fontFamily: 'TitilliumRegular',
    lineHeight: '20px'
  },
  tutorialDescriptionGrid: {
    marginTop: '-30px',
    marginLeft: '50px'
  },
  tutorialDescriptionText: {
    fontFamily: 'TitilliumRegular',
    color: '#2A3C67'
  },
  keyworksGrid: {
    marginTop: '-30px',
    marginLeft: '63px'
  },
  linkKeywords: {
    border: '1.5px solid rgba(127, 239, 239, 1)',
    lineHeight: '35px',
    color: '#2A3C67',
    cursor: 'pointer',
    fontStyle: 'normal',
    fontFamily: 'TitilliumRegular'
  },
  tutorialActions: {
    marginLeft: '50px'
  },
  tutorialResources: {
    marginTop: '-17px',
    marginLeft: '-6px'
  }
}))

const Accordion = withStyles({
  root: {
    borderTop: '10px solid rgba(127, 239, 239, 1)',
    scrollbarGutter: 'false',
    width: '100%',
    display: 'block',
    '&:not(:last-child)': {
    },
    '&:before': {
      display: 'none'
    },
    '&$expanded': {
      margin: 'auto'
    }
  },
  heading: {
    fontSize: 35,
    flexBasis: '33.33%',
    flexShrink: 0,
    marginLeft: '0px'
  },
  secondaryHeading: {
    fontSize: 10
  },
  expanded: {}
})(MUIAccordion)

const AccordionSummary = withStyles({
  root: {
    // flexDirection: 'column',
    // height: '150px',
    width: '100%',
    display: 'block'
  }
  // content: {
  //   marginBottom: 0,
  //   flexGrow: 1
  // },
  // expandIcon: {
  //   marginRight: '10px',
  //   paddingTop: '10px'
  // }
})(MUIAccordionSummary)

const AccordionDetails = withStyles({
  root: {
    width: '100%',
    display: 'block',
    marginTop: '30px'

  }
  // content: {
  //   marginBottom: 0,
  //   flexGrow: 1
  // },
  // expandIcon: {
  //   marginRight: '10px',
  //   paddingTop: '20px'
  // }
})(MUIAccordionDetails)

export default function AIToolkitPage() {
  const classes = useStyles()
  const [expanded, setExpanded] = useQueryParam('expanded', StringParam)
  const [queryParameters, setQueryParameters] = useQueryParams({
    author: StringParam, keyword: StringParam, method: StringParam, filterString: StringParam
  })
  const emptyQuery = {
    author: null,
    keyword: null,
    method: null,
    filterString: null
  }

  const filter = tutorial => {
    const {author, keyword, method} = queryParameters
    if (author && tutorial.authors.indexOf(author) === -1) {
      return false
    }
    if (keyword && tutorial.labels.application_keyword.indexOf(keyword) === -1) {
      return false
    }
    if (method && tutorial.labels.data_analytics_method.indexOf(method) === -1) {
      return false
    }
    return true
  }

  const tutorials_list = tutorials.tutorials.filter(tutorial => tutorial.labels.application_section[0] === 'Tutorials for artificial-intelligence methods')

  const {authors, keywords, methods} = useMemo(() => {
    const authors = {}
    const keywords = {}
    const methods = {}
    tutorials_list.forEach(tutorial => {
      tutorial.key = tutorial.title.replace(/\W/gm, '_').toLowerCase()
      tutorial.authors.forEach(i => { authors[i] = i })
      tutorial.labels.application_keyword.forEach(i => { keywords[i] = i })
      tutorial.labels.data_analytics_method.forEach(i => { methods[i] = i })
    }
    )
    return {
      authors: Object.keys(authors).sort(),
      keywords: Object.keys(keywords).sort(),
      methods: Object.keys(methods).sort()
    }
  }, [tutorials_list])

  return <Grid container spacing={1} className={classes.root}>
    <Grid container spacing={0}>
      <Grid item xs={6} className={classes.sectionTitle} >
        <Box className={classes.title}>
          {
            'Learn from tutorials'
          }
        </Box>
        <Box className={classes.deck}>
          {
            'We develop and implement methods that identify correlations and structure in big data of materials. This will enable scientists and engineers to decide which materials are useful for specific applications or which new materials should be the focus of future studies. The following tutorials are designed to get started with the AI Toolkit.'
          }
        </Box>
      </Grid>
      <Grid item xs={4} className={classes.sectionIcon}>
        <img src={TutorialsIcon} className={classes.icon}/>
      </Grid>
    </Grid>
    <Grid container spacing={0}>
      <Grid item xs={12} >
        <Box className={classes.filter} >
          {
            'Filter Tutorials'
          }
        </Box>
      </Grid>
      <Grid item xs={2}>
        <Autocomplete
          id="combo-box-demo"
          options={authors}
          className={classes.autocomplete}
          getOptionLabel={option => option}
          style={{height: '50px', width: '150px'}}
          renderInput={params => (
            <TextField {...params} label="Author" InputProps={{...params.InputProps, disableUnderline: true}} fullWidth />
          )}
          value={queryParameters.author}
          onChange={(_, value) => setQueryParameters({...emptyQuery, author: value})}
        />
      </Grid>
      <Grid item xs={2}>
        <Autocomplete
          id="combo-box-demo"
          options={keywords}
          style={{height: '50px', width: '150px'}}
          className={classes.autocomplete}
          getOptionLabel={option => option}
          renderInput={params => (
            <TextField {...params} label="Keyword" InputProps={{...params.InputProps, disableUnderline: true}} fullWidth />
          )}
          value={queryParameters.keyword}
          onChange={(_, value) => setQueryParameters({...emptyQuery, keyword: value})}
        />
      </Grid>
      <Grid item xs={2}>
        <Autocomplete
          id="combo-box-demo"
          options={methods}
          style={{height: '50px', width: '150px'}}
          className={classes.autocomplete}
          renderInput={params => (
            <TextField {...params} label="Method" InputProps={{...params.InputProps, disableUnderline: true}} fullWidth />
          )}
          value={queryParameters.method}
          onChange={(_, value) => setQueryParameters({...emptyQuery, method: value})}
        />
      </Grid>
    </Grid>

    <Grid container spacing={1} className={classes.tutorialsList}>
      <Grid item xs={12}>
        {tutorials_list.map(tutorial => (
          <div key={tutorial.title} >
            <Accordion
              key={tutorial.key}
              disabled={!filter(tutorial)}
              expanded={expanded === tutorial.key}
              onChange={() => setExpanded(expanded === tutorial.key ? null : tutorial.key)}
              className={classes.tutorial}
              elevation={0}
            >
              <AccordionSummary expandIcon={<ExpandMoreIcon />}>
                <Grid container spacing={1} >
                  <Grid item xs={7} className={classes.tutorialTitleGrid} >
                    <Typography className={classes.tutorialTitleText}>
                      {tutorial.title}
                    </Typography>
                  </Grid>
                  <Grid item xs={4} classes={classes.authorsGrid} >
                    <Typography className={classes.fieldText}>
                      {<b>Authors: </b> }
                      {tutorial.authors
                        .map(name => {
                          const label = name.split(',').reverse().join(' ')
                          return <Link
                            className={classes.linkAuthors}
                            key={name}
                            onClick={() => setQueryParameters({
                              ...emptyQuery,
                              author: queryParameters.author === name ? null : name
                            })}
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
                  <Grid item xs={6} className={classes.tutorialDescriptionGrid}>
                    <Markdown className={classes.tutorialDescriptionText}>
                      {tutorial.description}
                    </Markdown>
                  </Grid>
                  <Grid item xs={4} className={classes.keyworksGrid}>
                    <Typography className={classes.fieldText}>
                      <b>Keywords</b>:
                    </Typography>
                    <Typography>
                      {tutorial.labels.application_keyword
                        .map(keyword => (
                          <Link
                            className={classes.linkKeywords}
                            key={keyword}
                            onClick={() => setQueryParameters({
                              ...emptyQuery,
                              keyword: queryParameters.keyword === keyword ? null : keyword
                            })}
                          >
                            {keyword}
                          </Link>
                        )).reduce((prev, curr) => [prev, '    ', curr])
                      }
                    </Typography>
                    <Typography className={classes.fieldText}>
                      <b>Methods</b>:
                    </Typography>
                    <Typography>
                      {tutorial.labels.data_analytics_method
                        .map(method => (
                          <Link
                            className={classes.linkKeywords}
                            key={method}
                            onClick={() => setQueryParameters({
                              ...emptyQuery,
                              method: queryParameters.method === method ? null : method
                            })}
                          >
                            {method}
                          </Link>
                        )).reduce((prev, curr) => [prev, '    ', curr])
                      }
                    </Typography>
                  </Grid>
                </Grid>
              </AccordionDetails>

              <AccordionActions>
                <Grid container spacing={4}>
                  <Grid item xs={7} className={classes.tutorialActions}>
                    <Grid container spacing={0}>
                      <Grid item xs={5}>
                        <Button href={tutorial.link} target="tutorial" startIcon={<img src={AccessIcon}></img>}>
                          <Typography className={classes.fieldText} >
                            <b>Access tutorial</b>
                          </Typography>

                        </Button>
                      </Grid>

                      <Grid item xs={5} >

                        <Button width='10px' color='#2A3C67' href={tutorial.link_public} target="tutorial" startIcon={<img src={WatchIcon}></img>}>
                          <Typography className={classes.fieldText} >
                            <b>Watch video</b>
                          </Typography>                        </Button>
                      </Grid>
                    </Grid>

                  </Grid>
                  <Grid item xs={4} className={classes.tutorialResources}>
                    <Grid container spacing={0}>
                      <Grid item xs={12}>
                        <Typography className={classes.fieldText}>
                          <b>Additional Resources</b>:
                        </Typography>
                      </Grid>
                      <Grid item xs={2} >
                        <Button color='#2A3C67' href={tutorial.link} target="tutorial" startIcon={<img src={DoiIcon}></img>}>
                        </Button>
                      </Grid>
                      <Grid item xs={2}>
                        <Button color='#2A3C67' href={tutorial.link} target="tutorial" startIcon={<img src={PdfIcon}></img>}>
                        </Button>
                      </Grid>
                    </Grid>
                  </Grid>
                </Grid>
              </AccordionActions>
              <Divider />
            </Accordion>
          </div>
        ))}
      </Grid>
    </Grid>
  </Grid>
}
