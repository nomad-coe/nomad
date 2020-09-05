import React from 'react'
import { Card, CardHeader, CardContent, makeStyles } from '@material-ui/core'
import RawFiles from '../entry/RawFiles'

const useStyles = makeStyles(theme => ({
  root: {
    marginTop: theme.spacing(2)
  }
}))

export default function QCMSEntryCards(props) {
  const classes = useStyles()
  return (
    <Card className={classes.root}>
      <CardHeader title="Raw files" />
      <CardContent>
        <RawFiles {...props} />
      </CardContent>
    </Card>
  )
}
