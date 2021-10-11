import React from 'react'
import { Divider, Typography, makeStyles, Link, AccordionActions, Button, Grid } from '@material-ui/core'
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

const useStyles = makeStyles(theme => ({

  root: {
    margin: theme.spacing(3),
    width: '100%',
    marginLeft: 'auto',
    marginRight: 'auto',
    maxWidth: '1052px'
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
    fontSize: '35px',
    marginLeft: '-10px',
    fontFamily: 'TitilliumBold',
    marginTop: '-70px'
  },
  deck: {
    color: '#2A3C67',
    fontSize: '22px',
    marginTop: '20px',
    lineHeight: '30px',
    marginLeft: '-10px',
    fontFamily: 'TitilliumRegular',
    width: '518px'
  },
  icon: {
    height: '371px',
    marginTop: '-20px',
    marginLeft: '100px'
  },
  filter: {
    fontWeight: 'bold',
    color: '#2A3C67',
    fontSize: '20px',
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
    fontSize: '28px',
    color: '#2A3C67',
    fontFamily: 'TitilliumBold'
  },
  authorsGrid: {
    marginLeft: '150px',
    marginRight: '30px'
  },
  fieldText: {
    color: '#2A3C67'
  },
  linkAuthors: {
    color: '#2A3C67',
    cursor: 'pointer',
    fontFamily: 'TitilliumRegular',
    lineHeight: '20px',
    fontSize: '16px'
  },
  tutorialDescriptionGrid: {
    marginLeft: '50px'
  },
  tutorialDescriptionText: {
    fontFamily: 'TitilliumRegular',
    color: '#2A3C67',
    fontSize: '18px'
  },
  keyworksGrid: {
    marginLeft: '80px'
  },
  linkKeywords: {
    border: '1.5px solid rgba(127, 239, 239, 1)',
    lineHeight: '35px',
    color: '#2A3C67',
    cursor: 'pointer',
    fontStyle: 'normal',
    fontFamily: 'TitilliumRegular',
    fontSize: '16px'
  },
  tutorialActions: {
    marginLeft: '50px'
  },
  tutorialResources: {
    marginTop: '-17px',
    marginLeft: '-6px'
  },
  titleSecondary: {
    fontWeight: 'bold',
    color: 'rgba(127, 239, 239, 1)',
    fontSize: '35px',
    marginLeft: '-10px',
    fontFamily: 'TitilliumRegular'
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
    lineHeight: '20px',
    fontFamily: 'TitilliumBold'
  },
  bottomIcon: {
    height: '300px',
    marginTop: '80px',
    marginLeft: '120px'
  }
}))

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
