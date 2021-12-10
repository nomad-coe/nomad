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
import {
  Grid,
  Box,
  Button,
  Typography,
  Popover,
  IconButton,
  makeStyles
} from '@material-ui/core'
import { Link } from 'react-router-dom'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import IconQuery from '../../images/AIT_ico_bb_query.svg'
import IconReplicate from '../../images/AIT_ico_bb_replicate.svg'
import IconTutorial from '../../images/AIT_ico_bb_tutorial.svg'
import IconWork from '../../images/AIT_ico_bb_work.svg'
import IconQuery2 from '../../images/AIT_ico_bp_query.svg'
import IconReplicate2 from '../../images/AIT_ico_bp_replicate.svg'
import IconTutorial2 from '../../images/AIT_ico_bp_tutorial.svg'
import IconWork2 from '../../images/AIT_ico_bp_work.svg'
import ArrowIcon from '../../images/AIT_ico_bd_link_go_to.svg'
import FigureAI from '../../images/AIT_illu_AIT.svg'
import YouTubeEmbed from '../YouTubeEmbed'
import ScrollButton from '../buttons/ScrollButton'
import InfoIcon from '../../images/AIT_ico_bd_info_circle.svg'
import Background from '../../images/AIT_bg_title.jpg'
import ImgNatRev from '../../images/AIT_slide_nature.png'

const useStyles = makeStyles(theme => ({
  root: {
    margin: theme.spacing(3),
    width: '100%',
    marginLeft: 'auto',
    marginRight: 'auto',
    maxWidth: '1920px',
    marginBottom: '150px'
  },
  background: {
    position: 'relative',
    backgroundImage: `url(${Background})`,
    height: '840px',
    marginTop: '-50px',
    zIndex: 0
  },
  boxIcons: {
    width: '1000px',
    margin: 'auto',
    marginTop: '-150px'
  },
  title: {
    fontWeight: theme.typography.fontWeightMedium,
    fontSize: '55px',
    margin: 'auto',
    textAlign: 'center',
    align: 'center',
    marginTop: '100px',
    width: '650px',
    height: '140px',
    letterSpacing: 0,
    lineHeight: '62px',
    color: 'white'
  },
  deck: {
    letterSpacing: 0,
    marginTop: '-170px',
    wordSpacing: '5px',
    lineHeight: '42px',
    color: 'white',
    fontSize: '35px',
    margin: 'auto',
    textAlign: 'center',
    align: 'center',
    left: '736px',
    top: '270px',
    width: '550px',
    height: '140px'
  },
  topIcon: {
    width: '180px',
    marginBottom: '-10px'
  },
  topButton: {
    fontWeight: theme.typography.fontWeightMedium,
    backgroundColor: 'white',
    borderRadius: '30px',
    textTransform: 'none',
    marginTop: '-38px',
    fontSize: '20px',
    lineHeight: '20px',
    color: '#2A3C67',
    textAlign: 'center',
    align: 'center',
    paddingTop: '15px',
    paddingBottom: '15px',
    paddingRight: '35px',
    paddingLeft: '15px'
  },
  arrowGrid: {
    marginTop: '-32px',
    marginLeft: '-45px'
  },
  scrollButton: {
    position: 'absolute',
    bottom: '70px',
    left: theme.spacing(2),
    color: '#2A3C67',
    backgroundColor: theme.palette.background.paper,
    '&:hover, &.Mui-focusVisible': {
      color: '#2A3C67',
      backgroundColor: theme.palette.background.paper
    }
  },
  body: {
    width: '1052px',
    margin: theme.spacing(3),
    marginLeft: 'auto',
    marginRight: 'auto'
  },
  highlightedText: {
    fontWeight: theme.typography.fontWeightMedium,
    color: '#00DFE0',
    fontSize: '35px',
    width: '518px',
    lineHeight: '42px',
    marginTop: '80px'
  },
  bodyText: {
    color: '#2A3C67',
    lineHeight: '30px',
    fontSize: '22px',
    width: '607px',
    marginTop: '40px'
  },
  boxIconsBottom: {
    width: '1000px',
    margin: 'auto',
    marginTop: '50px'
  },
  buttonBottom: {
    fontWeight: theme.typography.fontWeightMedium,
    backgroundColor: '#F3F2F5',
    fontSize: '20px',
    lineHeight: '20px',
    color: '#2A3C67',
    textAlign: 'center',
    align: 'center',
    borderRadius: '30px',
    width: '207px',
    height: '70px',
    textTransform: 'none'
  },
  iconsBottom: {
    width: '200px'
  },
  toolTipText: {
    width: '400px',
    height: '520px',
    marginRight: '15px',
    marginLeft: '15px',
    marginTop: '10px'
  }
}))

export default function AIToolkitPage() {
  const styles = useStyles()

  const [anchorEl, setAnchorEl] = React.useState(null)
  const handleClick = (event) => {
    setAnchorEl(event.currentTarget)
  }
  const handleClose = () => {
    setAnchorEl(null)
  }
  const open = Boolean(anchorEl)
  const id = open ? 'simple-popover' : undefined

  return <Grid container spacing={2} className={styles.root}>
    <Grid container className={styles.background}>
      <ScrollButton scrollAmount={840} className={styles.scrollButton}>
        <ExpandMoreIcon/>
      </ScrollButton>
      <Grid item xs={12} >
        <Typography className={styles.title}>
          Artificial-Intelligence Tools for Materials Science
        </Typography>
      </Grid>
      <Grid item xs={12}>
        <Typography className={styles.deck}>
          Find new Patterns and Trends in Materials Science Big Data
        </Typography>
      </Grid>
      <Grid container spacing={1} className={styles.boxIcons}>
        <Grid item xs={3} >
          <IconButton href='https://analytics-toolkit.nomad-coe.eu/public/user-redirect/notebooks/tutorials/query_nomad_archive.ipynb'>
            <img
              src={IconQuery}
              className={styles.topIcon}
              style={{zIndex: 2, position: 'relative'}}
              alt='Query Archive icon'
            />
          </IconButton>
          <Grid container spacing={0}>
            <Grid item xs={11}>
              <Typography
                className={styles.topButton}
                style={{zIndex: 1, position: 'relative'}}
              >
                Query the Archive
              </Typography>
            </Grid>
            <Grid item xs={1} className={styles.arrowGrid}>
              <IconButton href='https://analytics-toolkit.nomad-coe.eu/public/user-redirect/notebooks/tutorials/query_nomad_archive.ipynb'>
                <img
                  src={ArrowIcon}
                  style={{width: '20px', zIndex: 4, position: 'relative'}}
                  alt = 'Arrow Icon'
                />
              </IconButton>
            </Grid>
          </Grid>
        </Grid>
        <Grid item xs={3} >
          <IconButton to="tutorials" component={Link}>
            <img
              alt='Tutorials icon'
              src={IconTutorial}
              className={styles.topIcon}
              style={{zIndex: 2, position: 'relative'}}
            />
          </IconButton>
          <Grid container spacing={0}>
            <Grid item xs={11}>
              <Typography
                className={styles.topButton}
                style={{zIndex: 1, position: 'relative'}}
              >
                View tutorials
              </Typography>
            </Grid>
            <Grid item xs={1} className={styles.arrowGrid}>
              <IconButton to="tutorials" component={Link}>
                <img
                  alt='Arrow Icon'
                  src={ArrowIcon}
                  style={{width: '20px', zIndex: 4, position: 'relative'}}
                />
              </IconButton>
            </Grid>
          </Grid>
        </Grid>
        <Grid item xs={3} >
          <IconButton to="reproduce" component={Link}>
            <img
              alt='Reproduce icon'
              src={IconReplicate}
              className={styles.topIcon}
              style={{zIndex: 2, position: 'relative'}}
            />
          </IconButton>
          <Grid container spacing={0}>
            <Grid item xs={11}>
              <Typography
                className={styles.topButton}
                style={{zIndex: 1, position: 'relative'}}
              >
                Published results
              </Typography>
            </Grid>
            <Grid item xs={1} className={styles.arrowGrid}>
              <IconButton to="reproduce" component={Link}>
                <img
                  alt='Arrow icon'
                  src={ArrowIcon}
                  style={{width: '20px', zIndex: 4, position: 'relative'}}
                />
              </IconButton>
            </Grid>
          </Grid>
        </Grid>
        <Grid item xs={3} >
          <Grid container >
            <Grid item xs={11}>
              <IconButton href="https://analytics-toolkit.nomad-coe.eu/hub/user-redirect/notebooks">
                <img
                  alt='Get to work icon'
                  src={IconWork}
                  className={styles.topIcon}
                  style={{zIndex: 2, position: 'relative'}}
                />
              </IconButton>
              <Grid container spacing={0}>
                <Grid item xs={11}>
                  <Typography
                    className={styles.topButton}
                    style={{zIndex: 1, position: 'relative'}}
                  >
                Get to work
                  </Typography>
                </Grid>
                <Grid item xs={1} className={styles.arrowGrid}>
                  <IconButton href="https://analytics-toolkit.nomad-coe.eu/hub/user-redirect/notebooks">
                    <img
                      alt='Arrow cion'
                      src={ArrowIcon}
                      style={{width: '20px', zIndex: 4, position: 'relative'}} />
                  </IconButton>
                </Grid>
              </Grid>
            </Grid>
            <Grid item xs={1} style={{marginTop: '161px', marginLeft: '-20px'}}>
              <IconButton aria-describedby={id} variant="contained" onClick={handleClick}>
                <img alt='Info icon' src={InfoIcon} />
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
                <Typography className={styles.toolTipText}>
                  By clicking on the 'Get to work' button you will access a
                  personal space that is available to each NOMAD user.  After
                  logging in, you will see a 'tutorials' and 'work' directory.
                  The 'tutorials' directory contains all notebooks available in
                  the AI toolkit, while the 'work' directory offers some space
                  to save personal work.  When you are in the 'work' directory,
                  click on the 'new' icon on the top right and then select
                  'Python 3'. This will create a Jupyter notebook that is
                  stored in the AI toolkit and can be reaccessed and iteratively
                  modified by the user.  All packages and software installed in the AI
                  toolkit are also available in the 'work' directory, that
                  makes it possible to employ the same code syntax used in
                  each tutorial contained in the AI toolkit for your own project.
                  For example, you can deploy any of the methodologies described in
                  tutorials on a different dataset. You can also upload your own data with
                  the 'Upload' button (via the menu bar on top, under the 'Publish' menu),
                  or directly access datasets in the NOMAD Archive. Make sure to learn
                  how to access the data in the NOMAD Archive, which is explained in the
                   tutorial accessible from the 'Query the Archive' button.
                </Typography>
              </Popover>
            </Grid>
          </Grid>
        </Grid>
      </Grid>
    </Grid>
    <Grid container spacing={1} className={styles.body}>
      <Grid item xs={8} >
        <Typography className={styles.highlightedText}>
          What is the NOMAD Artificial-Intelligence Toolkit?
        </Typography>
        <Typography className={styles.bodyText}>
          The preparation, synthesis, and characterization of new materials is a
          complex and costly aspect of materials design. The number of possible
          materials is practically infinite, about 200,000 materials are “known”
          to exist. But the basic properties (e.g., optical gap, elasticity
          constants, plasticity, piezoelectric tensors, conductivity, etc.) have
          been determined for very few of them.  NOMAD develops and provides a
          big set of tools - the Artificial-Intelligence Toolkit - using the
          latest artificial-intelligence approaches (including machine-learning,
          compressed sensing, and data mining) that make it possible to sort all
          available material data, to identify correlations and structures, and
          to detect trends and anomalies. Thus, the Artificial Intelligence Toolkit
          enables scientists and engineers to decide which materials are useful for
          specific applications or which new materials should be the focus of future
          studies.
        </Typography>
      </Grid>
      <Grid item xs={4}>
        <img alt='AI toolkit logo' src={FigureAI} style={{marginTop: '150px'}}/>

      </Grid>
      <Grid item xs={8} >
        <Typography className={styles.highlightedText}>
          How to get started
        </Typography>
        <Typography className={styles.bodyText}>
          Introduction to the scope of the NOMAD Artificial-Intelligence toolkit.
        </Typography>
        <div className="App" style={{marginTop: '30px'}}>
          <YouTubeEmbed embedId="v_Ie5TPXrd0" />
        </div>
        <Grid item xs={8}>
          <Typography className={styles.bodyText}>
          The NOMAD Artificial-Intelligence Toolkit is very accessible. Watch this video and
          learn more about its features.
          </Typography>
        </Grid>
      </Grid>
      <div className="App" style={{marginTop: '30px'}}>
        <YouTubeEmbed embedId="7R4EHsSRork" />
      </div>
      <Grid item xs={8} >
        <Typography className={styles.highlightedText}>
          Read about us!
        </Typography>
        <Typography className={styles.bodyText}>
          By clicking on the image below, you will access a Nature Reviews paper
          which gives an introduction to the NOMAD Artificial-Intelligence Toolkit.
        </Typography>
        <IconButton
          href='https://www.nature.com/articles/s42254-021-00373-8'
          style={{marginRight: '0px', marginTop: '20px'}}
        >
          <img alt='Nature logo' src={ImgNatRev}
            style={{width: '550px',
              marginTop: '15px',
              marginLeft: '-10px' }}
          />
        </IconButton>
      </Grid>

      <Grid item xs={8} >
        <Typography className={styles.highlightedText}>Access the tutorials </Typography>
        <Typography className={styles.bodyText}>
          Ready to start? Click on one of the options below. If you're new, we
          suggest starting with the tutorials.
        </Typography>
      </Grid>
    </Grid>

    <Grid container spacing={1} className={styles.boxIconsBottom}>
      <Grid item xs={3}>
        <IconButton href='https://analytics-toolkit.nomad-coe.eu/public/user-redirect/notebooks/tutorials/query_nomad_archive.ipynb'>
          <img alt='Query the Archive logo' src={IconQuery2} className={styles.iconsBottom}/>
        </IconButton>
      </Grid>
      <Grid item xs={3}>
        <img alt='Tutorials logo' src={IconTutorial2} className={styles.iconsBottom}/>
      </Grid>
      <Grid item xs={3}>
        <img alt='Reroduce logo' src={IconReplicate2} className={styles.iconsBottom}/>
      </Grid>
      <Grid item xs={3}>
        <img alt='Get to work log' src={IconWork2} className={styles.iconsBottom}/>
      </Grid>
      <Grid item xs={3}>
        <Button
          href='https://analytics-toolkit.nomad-coe.eu/public/user-redirect/notebooks/tutorials/query_nomad_archive.ipynb'
          className={styles.buttonBottom}
          endIcon={<img alt='Arrow icon' src={ArrowIcon}/>}
        >
          <Box className={styles.fieldText} >
            Query the Archive
          </Box>
        </Button>
      </Grid>
      <Grid item xs={3} >
        <Button
          component={Link}
          to="tutorials"
          className={styles.buttonBottom}
          endIcon={<img alt='Arrow icon' src={ArrowIcon}/>}
        >
          <Box className={styles.fieldText} >
            View tutorials
          </Box>
        </Button>
      </Grid>
      <Grid item xs={3}>
        <Button
          component={Link}
          to="reproduce"
          className={styles.buttonBottom}
          endIcon={<img alt='Arrow icon' src={ArrowIcon}/>}
        >
          <Box className={styles.fieldText}>
            Published results
          </Box>
        </Button>
      </Grid>
      <Grid item xs={3}>
        <Button
          href='https://analytics-toolkit.nomad-coe.eu/hub/user-redirect/notebooks'
          className={styles.buttonBottom}
          endIcon={<img alt='Arrow icon' src={ArrowIcon}/>}
        >
          <Box className={styles.fieldText}>
            Get to work
          </Box>
        </Button>
      </Grid>
    </Grid>
  </Grid>
}
