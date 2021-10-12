import React from 'react'
import {useStyles} from './styles.js'
import { Divider, Typography, Link, AccordionActions, Button, Grid } from '@material-ui/core'
import MuiAccordion from '@material-ui/core/Accordion'
import MuiAccordionSummary from '@material-ui/core/AccordionSummary'
import MuiAccordionDetails from '@material-ui/core/AccordionDetails'
import { styled } from '@material-ui/core/styles'
import ArrowForwardIosSharpIcon from '@material-ui/icons/ArrowForwardIosSharp'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import Markdown from '../Markdown'
import { StringParam, useQueryParam } from 'use-query-params'
import AccessIcon from './assets/AIT_ico_bd_link_external_big.svg'
import WatchIcon from './assets/AIT_ico_bd_youtube.svg'
import PdfIcon from './assets/AIT_ico_bd_link_pdf.svg'
import DoiIcon from './assets/AIT_ico_bd_link_doi.svg'

const Accordion = styled((props) => (
  <MuiAccordion disableGutters elevation={0} square {...props} />
))(({ theme }) => ({
  borderBottom: '13px solid rgba(127, 239, 239, 1)',
  '&:not(:last-child)': {
    borderBottom: 0
  },
  '&:before': {
    display: 'none'
  }
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
    marginLeft: theme.spacing(1)
  }
}))

const AccordionDetails = styled(MuiAccordionDetails)(({ theme }) => ({
  padding: theme.spacing(2)
}))

function AccordionsList(props) {
  const classes = useStyles()
  const [expanded, setExpanded] = useQueryParam('expanded', StringParam)
  return (
    props.tutorials_list.map(tutorial => (
      <div key={tutorial.title} >
        <Accordion
          key={tutorial.key}
          disabled={!props.filter(tutorial)}
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
                        onClick={() => props.setQueryParameters({
                          ...props.emptyQuery,
                          author: props.queryParameters.author === name ? null : name
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
                        onClick={() => props.setQueryParameters({
                          ...props.emptyQuery,
                          keyword: props.queryParameters.keyword === keyword ? null : keyword
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
                        onClick={() => props.setQueryParameters({
                          ...props.emptyQuery,
                          method: props.queryParameters.method === method ? null : method
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
                      </Typography>
                    </Button>
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
    )))
}

export default AccordionsList
