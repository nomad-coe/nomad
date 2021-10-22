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

import React, { useMemo, useState } from 'react'
import { useStyles } from './TutorialsPage'
import {
  Button,
  Grid,
  TextField,
  Typography,
  Divider,
  IconButton
} from '@material-ui/core'
import { Link } from 'react-router-dom'
import tutorials from '../../toolkitMetadata'
import { StringParam, useQueryParams } from 'use-query-params'
import Autocomplete from '@material-ui/lab/Autocomplete'
import TutorialsIcon from '../../images/AIT_ico_bb_tutorial.svg'
import ArrowIcon from '../../images/AIT_ico_bd_link_go_to.svg'
import ReproduceIcon from '../../images/AIT_ico_bp_replicate.svg'
import AccordionsList from './AccordionsList'
import FigureAI from '../../images/AIT_illu_AIT.svg'

export default function ReproducePage() {
  const styles = useStyles()
  const [queryParameters, setQueryParameters] = useQueryParams({
    author: StringParam, keyword: StringParam, method: StringParam, filterString: StringParam
  })

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

  const [valAuthor, setAuthor] = useState({})
  const [valKeyword, setKeyword] = useState({})
  const [valMethod, setMethod] = useState({})

  return <Grid container spacing={1} className={styles.root}>
    <Grid container spacing={0} className={styles.Heading}>
      <Grid item xs={6} className={styles.sectionTitle} >
        <Typography className={styles.title}>
          Reproduce published results
        </Typography>
        <Typography className={styles.deck}>
          We develop and implement methods that identify correlations and
          structure in big data of materials. This will enable scientists and
          engineers to decide which materials are useful for specific
          applications or which new materials should be the focus of future
          studies. The following notebooks are designed to reproduce result
          published in scientific journals.
        </Typography>
      </Grid>
      <Grid item xs={4} className={styles.sectionIcon}>
        <img alt='Reproduce icon' src={ReproduceIcon} className={styles.icon}/>
      </Grid>
    </Grid>
    <Grid container spacing={0}>
      <Grid item xs={12} >
        <Typography className={styles.filter} >
          Filter Tutorials
        </Typography>
      </Grid>
      <Grid item xs={3}>
        <Autocomplete
          id="combo-box-demo"
          options={authors}
          className={styles.autocomplete}
          getOptionLabel={option => option}
          renderInput={params => (
            <TextField
              {...params}
              label="Author"
              InputProps={{...params.InputProps, disableUnderline: true}}
              fullWidth
            />
          )}
          value={valAuthor}
          onChange={(_, value) => {
            setQueryParameters({author: value})
            setAuthor(value)
          }}
        />
      </Grid>
      <Grid item xs={3}>
        <Autocomplete
          id="combo-box-demo"
          options={keywords}
          className={styles.autocomplete}
          getOptionLabel={option => option}
          renderInput={params => (
            <TextField
              {...params}
              label="Keyword"
              InputProps={{...params.InputProps, disableUnderline: true}}
              fullWidth
            />
          )}
          value={valKeyword}
          onChange={(_, value) => {
            setQueryParameters({keyword: value})
            setKeyword(value)
          }}
        />
      </Grid>
      <Grid item xs={3}>
        <Autocomplete
          id="combo-box-demo"
          options={methods}
          className={styles.autocomplete}
          renderInput={params => (
            <TextField
              {...params}
              label="Method"
              InputProps={{...params.InputProps, disableUnderline: true}}
              fullWidth
            />
          )}
          value={valMethod}
          onChange={(_, value) => {
            setQueryParameters({method: value})
            setMethod(value)
          }}
        />
      </Grid>
    </Grid>
    <Grid container spacing={1} className={styles.tutorialsList}>

      <Grid item xs={12}>
        <Divider
          style={{
            backgroundColor: 'rgba(127, 239, 239, 1)',
            height: '13px',
            borderRadius: '4px',
            marginBottom: '-8px'
          }}
        />
      </Grid>
      <Grid item xs={12}>
        <AccordionsList tutorials_list={tutorials_list_advanced}
          author={authors}
          setAuthor = {setAuthor}
          keyword={keywords}
          setKeyword={setKeyword}
          method={methods}
          setMethod={setMethod}
          filter={filter}
          setQueryParameters={setQueryParameters}
          queryParameters={queryParameters}
          emptyQuery={queryParameters} />
      </Grid>

    </Grid>
    <Grid item xs={6} className={styles.sectionTitle} >
      <Typography className={styles.titleSecondary}>
        Explore the basics of AI
      </Typography>
      <Typography className={styles.deck}>
        Recent applications of artificial intelligence in science build on top
        of solid methodologies that have been being developed over the last decades.
        Exploring the following tutorials, you can learn the basics of AI to
        better understand their latest applications in materials science.
      </Typography>
      <Grid container spacing={1}>
        <Grid item xs={4}>
          <IconButton
            component={Link}
            to="aitoolkit"
            style={{marginRight: '0px', marginTop: '20px'}}
          >
            <img alt='AI toolkit logo' src={FigureAI} style={{width: '120px'}}/>
          </IconButton>
        </Grid>
        <Grid item xs={8}>
          <Button
            width='10px'
            color='#2A3C67'
            component={Link}
            to="tutorials"
            className={styles.bottomButton}
            endIcon={<img alt='Arrow icon' src={ArrowIcon}/>}
          >
            <Typography className={styles.fieldText} >
              AI tutorials
            </Typography>
          </Button>
        </Grid>
      </Grid>
    </Grid>
    <Grid item xs={4} className={styles.sectionIcon}>
      <IconButton
        component={Link}
        to="tutorials"
        className={styles.bottomIcon}
      >
        <img alt='Tutorials icon' src={TutorialsIcon} style={{width: '300px'}}/>
      </IconButton>

    </Grid>
  </Grid>
}
