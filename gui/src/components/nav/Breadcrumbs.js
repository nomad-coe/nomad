import React from 'react'
import { makeStyles } from '@material-ui/core/styles'
import {
  Link,
  Typography
} from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { Link as RouterLink } from 'react-router-dom'
import HelpDialog from '../Help'
import { useRoute } from './Routes'

const useStyles = makeStyles((theme) => ({
  root: {
    display: 'flex',
    flexDirection: 'row',
    alignItems: 'center'
  },
  submenuName: {
    fontSize: '1.30rem'
  },
  helpButton: {
    marginTop: '0.1rem',
    marginLeft: theme.spacing(1)
  }
}))

/**
 * Shows the current page title together with navigation breadcrumbs.
 */
function Breadcrumbs({className}) {
  const styles = useStyles()
  const {help, title, breadCrumbs} = useRoute()

  return (
    <div className={clsx(styles.root, className)}>
      <Typography variant="h6" color="primary" noWrap className={styles.submenuName}>
        {breadCrumbs && breadCrumbs.map((breadCrumb, index) => <React.Fragment key={index}>
          <Link component={RouterLink} to={breadCrumb.path}>{breadCrumb.title}</Link>&nbsp;›&nbsp;
        </React.Fragment>)}
        {title}
      </Typography>
      {help
        ? <HelpDialog
          color="primary"
          maxWidth="md"
          classes={{root: styles.helpButton}}
          size="small"
          {...help}
        />
        : ''
      }
    </div>
  )
}

Breadcrumbs.propTypes = {
  className: PropTypes.string
}

export default Breadcrumbs
