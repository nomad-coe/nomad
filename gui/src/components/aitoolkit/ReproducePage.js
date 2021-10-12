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
import {useStylesList} from './styles.js'
import { Box, Button, Grid, TextField, Divider } from '@material-ui/core'
import tutorials from '../../toolkitMetadata'
import { StringParam, useQueryParams } from 'use-query-params'
import Autocomplete from '@material-ui/lab/Autocomplete'
import TutorialsIcon from './assets/AIT_ico_bb_tutorial.svg'
import ArrowIcon from './assets/AIT_ico_bd_link_go_to.svg'
import ReproduceIcon from './assets/AIT_ico_bp_replicate.svg'
import AccordionsList from './AccordionsList'

export default function AIToolkitPage() {
  const classes = useStylesList()
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

  const tutorials_list_advanced = tutorials.tutorials.filter(tutorial => tutorial.labels.category[0] === 'advanced_tutorial')

  const {authors, keywords, methods} = useMemo(() => {
    const authors = {}
    const keywords = {}
    const methods = {}
    tutorials_list_advanced.forEach(tutorial => {
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
  }, [tutorials_list_advanced])

  return <Grid container spacing={1} className={classes.root}>
    <Grid container spacing={0} className={classes.Heading}>
      <Grid item xs={6} className={classes.sectionTitle} >
        <Box className={classes.title}>
          {
            'Reproduce published results'
          }
        </Box>
        <Box className={classes.deck}>
          {
            'We develop and implement methods that identify correlations and structure in big data of materials. This will enable scientists and engineers to decide which materials are useful for specific applications or which new materials should be the focus of future studies. The following tutorials are designed to reproduce result published in scientific journals.'
          }
        </Box>
      </Grid>
      <Grid item xs={4} className={classes.sectionIcon}>
        <img src={ReproduceIcon} className={classes.icon}/>
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
        <Divider disableGutters className={classes.tutorialsDivider}></Divider>
      </Grid>
      <Grid item xs={12}>
        <AccordionsList tutorials_list={tutorials_list_advanced}
          author={authors}
          keyword={keywords}
          method={methods}
          filter={filter}
          setQueryParameters={setQueryParameters}
          queryParameters={queryParameters}
          emptyQuery={queryParameters} />
      </Grid>

    </Grid>
    <Grid item xs={6} className={classes.sectionTitle} >
      <Box className={classes.titleSecondary}>
        {
          'Explore the basics of AI'
        }
      </Box>
      <Box className={classes.deck}>
        {
          'Recent applications of artificial-intelligence in science build on top of solid methodologies, that have been developing over last decades. You can learn the basics of AI to better understand its latest applications in materials science.'
        }
      </Box>
      <Button width='10px' color='#2A3C67' href={'tutorials'} target="tutorial" className={classes.bottomButton} endIcon={<img src={ArrowIcon}></img>}>
        <Box className={classes.fieldText} >
          AI tutorials
        </Box>
      </Button>
    </Grid>
    <Grid item xs={4} className={classes.sectionIcon}>
      <img src={TutorialsIcon} className={classes.bottomIcon}/>
    </Grid>

  </Grid>
}
