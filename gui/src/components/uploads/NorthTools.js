import React from 'react'
import { Button, ButtonGroup, Grid } from '@material-ui/core'

const NorthTools = React.memo(function NorthTools(props) {
  return <Grid container spacing={2}>
    <Grid item>
      <ButtonGroup variant="contained">
        <Button color="primary">Open Jupyter</Button>
        <Button>Delete Container</Button>
      </ButtonGroup>
    </Grid>
    {['Hyperspy', 'Nionswift'].map((tool, index) => (
      <Grid key={index} item>
        <Button color="primary" variant="contained">Launch {tool}</Button>
      </Grid>
    ))}
  </Grid>
})
NorthTools.propTypes = {}

export default NorthTools
