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
  Button,
  Grid,
  Typography,
  Divider,
  IconButton,
  makeStyles,
  AccordionActions

} from '@material-ui/core'
import { Link } from 'react-router-dom'
import { StringParam, useQueryParam } from 'use-query-params'
import TutorialsIcon from '../../images/AIT_ico_bp_tutorial.svg'
import ArrowIcon from '../../images/AIT_ico_bd_link_go_to.svg'
import FigureAI from '../../images/AIT_illu_AIT.svg'
import tutorials from '../../courseAI.json'
import { aitoolkitEnabled } from '../../config'
import MuiAccordion from '@material-ui/core/Accordion'
import MuiAccordionSummary from '@material-ui/core/AccordionSummary'
import MuiAccordionDetails from '@material-ui/core/AccordionDetails'
import { styled } from '@material-ui/core/styles'
import ArrowForwardIosSharpIcon from '@material-ui/icons/ArrowForwardIosSharp'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import Markdown from '../Markdown'
import AccessIcon from '../../images/AIT_ico_bd_link_external_big.svg'
import WatchIcon from '../../images/AIT_ico_bd_youtube.svg'
import PdfIcon from '../../images/AIT_ico_bd_link_pdf.svg'
import DoiIcon from '../../images/AIT_ico_bd_link_doi.svg'

export const useStyles = makeStyles(theme => ({
  root: {
    margin: theme.spacing(3),
    width: '100%',
    marginLeft: 'auto',
    marginRight: 'auto',
    maxWidth: '1052px',
    marginBottom: '150px'
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
    color: '#2A3C67',
    fontSize: '35px',
    marginLeft: '-10px',
    fontWeight: theme.typography.fontWeightMedium,
    marginTop: '-70px'
  },
  deck: {
    color: '#2A3C67',
    fontSize: '22px',
    marginTop: '20px',
    lineHeight: '30px',
    marginLeft: '-10px',
    width: '518px'
  },
  icon: {
    height: '371px',
    marginTop: '-20px',
    marginLeft: '100px'
  },
  filter: {
    fontWeight: theme.typography.fontWeightMedium,
    color: '#2A3C67',
    fontSize: '20px',
    marginTop: '60px',
    marginLeft: '0px'
  },
  autocomplete: {
    height: 'auto',
    color: '#2A3C67',
    border: '3px solid #00DFE0',
    borderRadius: '10px 10px 10px 10px',
    marginTop: '10px',
    marginLeft: '0px',
    width: '240px'
  },
  tutorialsList: {
    marginTop: '50px'
  },
  textLevel: {
    textAlign: 'left',
    color: '#2A3C67',
    fontSize: '22px',
    height: '22px',
    marginTop: '-16px'
  },
  titleSecondary: {
    fontWeight: 'bold',
    color: '#00DFE0',
    fontSize: '35px',
    marginLeft: '-10px'
  },
  bottomButton: {
    color: '#F3F2F5',
    backgroundColor: '#F3F2F5',
    borderRadius: '30px',
    width: '242px',
    height: '70px',
    textAlign: 'center',
    align: 'center',
    marginTop: '40px',
    textTransform: 'none',
    fontSize: '12pt',
    lineHeight: '20px'
  },
  bottomButtonText: {
    color: '#2A3C67',
    fontWeight: theme.typography.fontWeightMedium
  },
  bottomIcon: {
    marginTop: '80px',
    marginLeft: '120px'
  },
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

export default function TutorialsPage() {
  const styles = useStyles()

  const tutorials_list = tutorials.tutorials

  tutorials_list.forEach(tutorial => {
    tutorial.key = tutorial.title.replace(/\W/gm, '_').toLowerCase()
  })

  function AccordionsList() {
    const [expanded, setExpanded] = useQueryParam('expanded', StringParam)
    return (
      tutorials_list.map(tutorial => (
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
                    {<b>Lecturer: </b> }
                    {tutorial.lecturer

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

                  {tutorial.notebook_name_1 && <><Typography className={styles.fieldText}>
                    <b>Notebook authors:</b>
                  </Typography><Typography>
                    {tutorial.authors_notebook_1
                      .map(name => {
                        const label = name.split(',').reverse().join(' ')
                        return label
                      }).reduce((prev, curr) => [prev, ' | ', curr])}
                  </Typography>
                  <Button
                    href={tutorial.link_notebook_1}
                    target="tutorial"
                    startIcon={<img alt='Access icon' src={AccessIcon}></img>}
                  >
                    <Typography className={styles.fieldText} >
                      <b>Access notebook</b>
                    </Typography>
                  </Button>
                  </>}
                  {tutorial.notebook_name_2 && <><Typography className={styles.fieldText}>
                    <b>Notebook authors:</b>
                  </Typography><Typography>
                    {tutorial.authors_notebook_2
                      .map(name => {
                        const label = name.split(',').reverse().join(' ')
                        return label
                      }).reduce((prev, curr) => [prev, ' | ', curr])}
                  </Typography>
                  <Button
                    href={tutorial.link_notebook_2}
                    target="tutorial"
                    startIcon={<img alt='Access icon' src={AccessIcon}></img>}
                  >
                    <Typography className={styles.fieldText} >
                      <b>Access notebook</b>
                    </Typography>
                  </Button>
                  </>}
                </Grid>
              </Grid>
            </AccordionDetails>
            <AccordionActions>
              <Grid container spacing={4}>
                <Grid item xs={7} className={styles.tutorialActions}>
                  <Grid container spacing={0}>
                    <Grid item xs={5}>
                      <Button
                        href={tutorial.link_video_1}
                        target="tutorial"
                        startIcon={<img alt='Watch icon' src={WatchIcon}></img>}
                      >
                        <Typography className={styles.fieldText} >
                          <b>Watch video 1</b>
                        </Typography>
                      </Button>
                    </Grid>
                    <Grid item xs={5} >
                      <div>
                        { tutorial.link_video_2 && <Button
                          width='10px'
                          color='#2A3C67'
                          href={tutorial.link_video_2}
                          target="tutorial"
                          startIcon={<img alt='Watch icon' src={WatchIcon}></img>}
                        >
                          <Typography className={styles.fieldText} >
                            <b>Watch video 2</b>
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

  return <Grid container spacing={1} className={styles.root}>
    <Grid container spacing={0} className={styles.Heading}>
      <Grid item xs={6} className={styles.sectionTitle}>
        <Grid container spacing={0}>
          <Grid item xs={4} style={{marginTop: '-100px', marginLeft: '-20px'}}>
            <IconButton
              {...(aitoolkitEnabled ? ({to: 'aitoolkit', component: Link}) : ({href: 'https://nomad-lab.eu/AIToolkit', component: 'a'}))}
            >
              <img alt='AI toolkit logo' src={FigureAI} style={{width: '120px'}}/>
            </IconButton>
          </Grid>
          <Grid item xs={8}>
            <Typography className={styles.title}>
              AI lectures
            </Typography>
          </Grid>
        </Grid>
        <Typography className={styles.deck}>
        We present a virtual course in artificial intelligence, by proposing a path through some of the tutorial notebooks from the AI-toolkit. Each topic is presented via a video lecture and accompanied by one or more notebooks that allow for putting your hand on the explained methods.        </Typography>
      </Grid>
      <Grid item xs={4} className={styles.sectionIcon}>
        <img alt='Tutorials icon' src={TutorialsIcon} className={styles.icon}/>
      </Grid>
    </Grid>

    <Grid container spacing={1} className={styles.tutorialsList}>
      <Grid item xs={12}>
        <Divider
          style={{
            backgroundColor: '#7FEFEF',
            height: '13px',
            borderRadius: '4px',
            marginBottom: '-8px'
          }}
        />
      </Grid>
      <Grid item xs={12}>
        <AccordionsList/>
      </Grid>

    </Grid>
    <Grid item xs={6} className={styles.sectionTitle}>
      <Typography className={styles.titleSecondary}>
        Next intermediate level
      </Typography>
      <Typography className={styles.deck}>
        If you are still curious about applications of artificial intelligence to materials science, you can find
        more lectures and more notebooks in our tutorials section.
      </Typography>
      <Grid container spacing={1}>
        <Grid item xs={4}>
          <IconButton
            component={Link}
            {...(aitoolkitEnabled ? ({to: 'aitoolkit', component: Link}) : ({href: 'https://nomad-lab.eu/AIToolkit', component: 'a'}))}
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
            <Typography className={styles.bottomButtonText}>
              More AI tutorials
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
