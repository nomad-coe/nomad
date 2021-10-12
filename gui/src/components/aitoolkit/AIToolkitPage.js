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
import { Grid, Box, Button } from '@material-ui/core'
import IconQuery from './assets/AIT_ico_bb_query.svg'
import IconReplicate from './assets/AIT_ico_bb_replicate.svg'
import IconTutorial from './assets/AIT_ico_bb_tutorial.svg'
import IconWork from './assets/AIT_ico_bb_work.svg'
import {useStylesLanding} from './styles.js'
import ArrowIcon from './assets/AIT_ico_bd_link_go_to.svg'

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
      <Grid item xs={3} style={{zIndex: 1}}> <img src={IconQuery} />
      </Grid>
      <Grid item xs={3} style={{zIndex: 1}}> <img src={IconTutorial} />
      </Grid>
      <Grid item xs={3} style={{zIndex: 1}}> <img src={IconReplicate} />
      </Grid>
      <Grid item xs={3} style={{zIndex: 1}}> <img src={IconWork} />
      </Grid>
      <Grid item xs={3}>
        <Button width='10px' color='#2A3C67' href='reproduce' className={classes.button} endIcon={<img src={ArrowIcon}></img>}>
          <Box className={classes.fieldText} >
          Query the Archive
          </Box>
        </Button>
      </Grid>
      <Grid item xs={3} >
        <Button width='10px' color='#2A3C67' href='tutorials' className={classes.button} endIcon={<img src={ArrowIcon}></img>}>
          <Box className={classes.fieldText} >
          View tutorials
          </Box>
        </Button>
      </Grid>
      <Grid item xs={3}>
        <Button width='10px' color='#2A3C67' href='reproduce' className={classes.button} endIcon={<img src={ArrowIcon}></img>}>
          <Box className={classes.fieldText} >
          Reproduce published results
          </Box>
        </Button>
      </Grid>
      <Grid item xs={3}>
        <Button width='10px' color='#2A3C67' className={classes.button} endIcon={<img src={ArrowIcon}></img>}>
          <Box className={classes.fieldText} >
          Get to work
          </Box>
        </Button>

      </Grid>
    </Grid>
  </Grid>
}
