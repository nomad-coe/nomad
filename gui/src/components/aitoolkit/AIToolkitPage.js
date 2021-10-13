/* eslint-disable react/no-unescaped-entities */
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
import { Grid, Box, Button, Typography, Popover, IconButton } from '@material-ui/core'
import IconQuery from './assets/AIT_ico_bb_query.svg'
import IconReplicate from './assets/AIT_ico_bb_replicate.svg'
import IconTutorial from './assets/AIT_ico_bb_tutorial.svg'
import IconWork from './assets/AIT_ico_bb_work.svg'
import IconQuery2 from './assets/AIT_ico_bp_query.svg'
import IconReplicate2 from './assets/AIT_ico_bp_replicate.svg'
import IconTutorial2 from './assets/AIT_ico_bp_tutorial.svg'
import IconWork2 from './assets/AIT_ico_bp_work.svg'
import {useStylesLanding} from './styles.js'
import ArrowIcon from './assets/AIT_ico_bd_link_go_to.svg'
import FigureAI from './assets/AIT_illu_AIT.svg'
import YoutubeEmbed from './YoutubeEmbed'
import InfoIcon from './assets/AIT_ico_bd_info_circle.svg'

export default function AIToolkitPage() {
  const [anchorEl, setAnchorEl] = React.useState(null)

  const handleClick = (event) => {
    setAnchorEl(event.currentTarget)
  }

  const handleClose = () => {
    setAnchorEl(null)
  }

  const open = Boolean(anchorEl)
  const id = open ? 'simple-popover' : undefined
  const classes = useStylesLanding()

  return <Grid container spacing={2} className={classes.root}>
    <Grid container className={classes.background}>
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
          <Button color='#2A3C67' href='https://analytics-toolkit.nomad-coe.eu/public/user-redirect/notebooks/tutorials/query_nomad_archive.ipynb' className={classes.button} endIcon={<img src={ArrowIcon}></img>}>
            <Box className={classes.fieldText} > Query the Archive </Box>
          </Button>
        </Grid>
        <Grid item xs={3} >
          <Button color='#2A3C67' href='tutorials' className={classes.button} endIcon={<img src={ArrowIcon}></img>}>
            <Box className={classes.fieldText} > View tutorials </Box>
          </Button>
        </Grid>
        <Grid item xs={3}>
          <Button href='reproduce' className={classes.button} endIcon={<img src={ArrowIcon}></img>}>
            <Box className={classes.fieldText} > Published results </Box>
          </Button>
        </Grid>
        <Grid item xs={3}>
          <Grid container >
            <Grid item xs={11}>
              <Button color='#2A3C67' href='https://analytics-toolkit.nomad-coe.eu/hub/user-redirect/notebooks' className={classes.button} endIcon={<img src={ArrowIcon}></img>}>
                <Box className={classes.fieldText} > Get to work </Box>
              </Button>
            </Grid>
            <Grid item xs={1} style={{marginTop: '-35px', marginLeft: '-20px'}}>
              <IconButton aria-describedby={id} variant="contained" onClick={handleClick}>
                <img src={InfoIcon} ></img>
              </IconButton>
              <Popover
                id={id}
                open={open}
                anchorEl={anchorEl}
                onClose={handleClose}
                anchorOrigin={{
                  vertical: 'bottom',
                  horizontal: 'left'
                }}
              >
                <Typography style={{width: '400px', height: '500px', marginRight: '15px', marginLeft: '15px', marginTop: '10px', fontFamily: 'TitilliumRegular'}}>
                  Clicking on the 'Get to work' button will redirect to a personal space that is available to each NOMAD user.
                  After logging in, you can access a 'tutorials' and 'work' directory.
                  The 'tutorials' directory contains all notebooks available in the AI toolkit, while the 'work' directory offers some space to save personal work.
                  When you are in the 'work' directory, click on the 'new' icon on the top right and then select 'Python 3'.
                  This will create a Jupyter notebook that is stored in the AI toolkit and can be reaccessed and iteratively modified by the user.
                  All packages installed in the AI toolkit are also available in the 'work' directory, thus making possible to employ the same code syntax that is used in each tutorial presented in the AI toolkit.
                  It is possible, for example, to deploy any of the methodologies described in tutorials on a novel dataset.
                  You can also upload your own data with the 'Upload' button on the top right, or even download a novel dataset from the NOMAD Archive.
                  Make sure to learn how to download NOMAD data, that is explained in the tutorial accessible from the 'Query the Archive' button.
                </Typography>
              </Popover>
            </Grid>
          </Grid>
        </Grid>
      </Grid>
    </Grid>
    <Grid container spacing={1} className={classes.body}>
      <Grid item xs={8} >
        <Typography className={classes.highlightedText}>What is the NOMAD Artificial Intelligence Toolkit? </Typography>

        <Typography className={classes.bodyText}>
          The preparation, synthesis, and characterization of new materials is a complex and costly aspect of materials design.
          The number of possible materials is practically infinite, about 200,000 materials are “known” to exist.
          But the basic properties (e.g., optical gap, elasticity constants, plasticity, piezoelectric tensors, conductivity, etc.) have been determined for very few of them.
          NOMAD develops and provides a big set of tools- the Artificial Intelligence Toolkit - using the latest machine learning and artificial intelligence approaches that make it possible to sort all available material data, to identify correlations and structures, and to detect trends and anomalies.
          Thus, the Artificial Intelligence Toolkit enables scientists and engineers to decide which materials are useful for specific applications or which new materials should be the focus of future studies. </Typography>
      </Grid>

      <Grid item xs={4}>
        <img src={FigureAI} style={{marginTop: '150px'}}></img>
      </Grid>
      <Grid item xs={8} >
        <Typography className={classes.highlightedText}>How to get started </Typography>
        <div className="App" style={{marginTop: '30px'}}>
          <YoutubeEmbed embedId="v_Ie5TPXrd0" />
        </div>

        <Typography className={classes.bodyText}>
          NOMAD Artificial Toolkit is very accessible.
          Watch this video and learn more about its features. Our philosophy is... </Typography>
      </Grid>
      <div className="App" style={{marginTop: '30px'}}>
        <YoutubeEmbed embedId="7R4EHsSRork" />
      </div>

      <Grid item xs={8} >
        <Typography className={classes.highlightedText}>Access the tutorials </Typography>

        <Typography className={classes.bodyText}>Ready to start? Click on one of the options below. If you're new, we suggest starting with the tutorials. </Typography>
      </Grid>
    </Grid>

    <Grid container spacing={1} className={classes.boxIconsBottom}>
      <Grid item xs={3}> <img src={IconQuery2} className={classes.iconsBottom}/>
      </Grid>
      <Grid item xs={3}> <img src={IconTutorial2} className={classes.iconsBottom}/>
      </Grid>
      <Grid item xs={3}> <img src={IconReplicate2} className={classes.iconsBottom}/>
      </Grid>
      <Grid item xs={3}> <img src={IconWork2} className={classes.iconsBottom}/>
      </Grid>
      <Grid item xs={3}>
        <Button color='#2A3C67' href='https://analytics-toolkit.nomad-coe.eu/public/user-redirect/notebooks/tutorials/query_nomad_archive.ipynb' className={classes.buttonBottom} endIcon={<img src={ArrowIcon}></img>}>
          <Box className={classes.fieldText} >
          Query the Archive
          </Box>
        </Button>
      </Grid>
      <Grid item xs={3} >
        <Button color='#2A3C67' href='tutorials' className={classes.buttonBottom} endIcon={<img src={ArrowIcon}></img>}>
          <Box className={classes.fieldText} >
          View tutorials
          </Box>
        </Button>
      </Grid>
      <Grid item xs={3}>
        <Button href='reproduce' className={classes.buttonBottom} endIcon={<img src={ArrowIcon}></img>}>
          <Box className={classes.fieldText} >
          Published results
          </Box>
        </Button>
      </Grid>
      <Grid item xs={3}>
        <Button href='https://analytics-toolkit.nomad-coe.eu/hub/user-redirect/notebooks' color='#2A3C67' className={classes.buttonBottom} endIcon={<img src={ArrowIcon}></img>}>
          <Box className={classes.fieldText} >
          Get to work
          </Box>
        </Button>
      </Grid>
    </Grid>
  </Grid>
}
