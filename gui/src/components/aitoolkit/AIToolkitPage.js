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
import { Grid, Box } from '@material-ui/core'
import IconQuery from './assets/AIT_ico_bb_query.svg'
import IconReplicate from './assets/AIT_ico_bb_replicate.svg'
import IconTutorial from './assets/AIT_ico_bb_tutorial.svg'
import IconWork from './assets/AIT_ico_bb_work.svg'
import {useStylesLanding} from './styles.js'

export default function AIToolkitPage() {
  const classes = useStylesLanding()

  return <Grid container spacing={2} className={classes.root}>

    <Grid item xs={12} >
      <Box className={classes.title}>
         Powerful Tools for Materials Science  </Box>
    </Grid>
    <Grid item xs={12}>
      <Box className={classes.deck}>
      Find new Patterns and Information in Materials Science Big Data  </Box>
    </Grid>
    <Grid container spacing={1} className={classes.boxIcons}>
      <Grid item xs={3}> <img src={IconQuery} />
      </Grid>
      <Grid item xs={3}> <img src={IconTutorial} />
      </Grid>
      <Grid item xs={3}> <img src={IconReplicate} />
      </Grid>
      <Grid item xs={3}> <img src={IconWork} />
      </Grid>
      <Grid item xs={3}>
        <button className={classes.button}
          type="button"
        > Query the archive</button>
      </Grid>
      <Grid item xs={3} >
        <button className={classes.button}
          type="button"
          onClick={(e) => {
            e.preventDefault()
            window.location.href = 'tutorials'
          }}
        > View tutorials</button>
      </Grid>
      <Grid item xs={3}>
        <button className={classes.button}
          type="button"
          onClick={(e) => {
            e.preventDefault()
            window.location.href = 'reproduce'
          }}
        > Reproduce published results</button>
      </Grid>
      <Grid item xs={3}>
        <button className={classes.button}
          type="button"
        > Get to work</button>
      </Grid>
    </Grid>
  </Grid>
}
