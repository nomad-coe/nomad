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

import React, { useState, useCallback } from 'react'
import { makeStyles } from '@material-ui/core/styles'
import { useRoute } from './Routes'
import { useCookies } from 'react-cookie'
import { MenuBar, MenuBarItem, MenuBarMenu } from './MenuBar'
import Consent from './Consent'
import { appBase, oasis, aitoolkitEnabled } from '../../config'

import BackupIcon from '@material-ui/icons/Backup'
import SearchIcon from '@material-ui/icons/Search'
import UserDataIcon from '@material-ui/icons/AccountCircle'
import AboutIcon from '@material-ui/icons/Home'
import ForumIcon from '@material-ui/icons/QuestionAnswer'
import FAQIcon from '@material-ui/icons/LiveHelp'
import MetainfoIcon from '@material-ui/icons/Info'
import DocIcon from '@material-ui/icons/Help'
import CodeIcon from '@material-ui/icons/Code'
import TermsIcon from '@material-ui/icons/Assignment'
import AnalyticsIcon from '@material-ui/icons/ShowChart'

const useStyles = makeStyles(theme => ({
  spacer: {
    flexGrow: 1
  }
}))

/**
 * Main menu showing the available navigation items.
 */
const MainMenu = React.memo(() => {
  const route = useRoute()
  const cookies = useCookies()[0]
  const styles = useStyles()
  const selected = (route?.navPath) || 'publish/uploads'
  const [consentOpen, setConsentOpen] = useState(cookies['terms-accepted'] !== 'true')

  const handleConsentAccept = useCallback(() => {
    setConsentOpen(false)
  }, [])

  return <MenuBar selected={selected}>
    <MenuBarMenu name="publish" label="Publish" route="/uploads" icon={<BackupIcon/>}>
      <MenuBarItem
        name="uploads" label="Upload" route="/uploads" isDefault
        tooltip="Upload and publish new data" icon={<SearchIcon />}
      />
      <MenuBarItem
        label="Your data" name="userdata" route="/userdata"
        tooltip="Manage your uploaded data" icon={<UserDataIcon />}
      />
    </MenuBarMenu>
    <MenuBarMenu name="explore" route="/search" icon={<SearchIcon/>}>
      <MenuBarItem
        label="Search Entries" name="searchentries" route="/search/entries"
        tooltip="Search individual entries"
      />
      <MenuBarItem
        label="Search Materials" name="searchmaterials" route="/search/materials"
        tooltip="Search materials"
      />
    </MenuBarMenu>
    <MenuBarMenu name="analyze" route="/metainfo" icon={<AnalyticsIcon/>}>
      {(!oasis && aitoolkitEnabled)
        ? <MenuBarItem
          label="AI Toolkit" name="aitoolkit" route="/aitoolkit"
          tooltip="NOMAD's Artificial Intelligence Toolkit tutorial Jupyter notebooks"
          icon={<MetainfoIcon />}
        />
        : <MenuBarItem
          label="AI Toolkit" name="aitoolkit"
          href="https://nomad-lab.eu/AIToolkit"
          tooltip="Visit the NOMAD Artificial Intelligence Analytics Toolkit"
        />
      }
      <MenuBarItem
        name="metainfo" label="The NOMAD Metainfo" route="/metainfo" tooltip="Browse the NOMAD Metainfo Schema"
      />
      <MenuBarItem
        name="apis" label="APIs" route="/apis" tooltip="The list of APIs offered by NOMAD"
      />
    </MenuBarMenu>
    <MenuBarMenu name="about" route="/" icon={<AboutIcon/>}>
      <MenuBarItem
        label="Information" name="info" route="/"
        tooltip="Overview of the NOMAD Repository and Archive"
      />
      <MenuBarItem
        name="forum"
        href="https://matsci.org/c/nomad/"
        tooltip="The NOMAD user/developer forum on matsci.org"
        icon={<ForumIcon/>}
      />
      <MenuBarItem
        label="FAQ" name="faq"
        href="https://nomad-lab.eu/repository-archive-faqs"
        tooltip="Frequently Asked Questions (FAQ)"
        icon={<FAQIcon/>}
      />
      <MenuBarItem
        name="Docs"
        href={`${appBase}/docs/index.html`}
        tooltip="The full user and developer documentation"
        icon={<DocIcon/>}
      />
      <MenuBarItem
        name="Sources"
        href="https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR"
        tooltip="NOMAD's main Gitlab project"
        icon={<CodeIcon/>}
      />
      <MenuBarItem
        name="Terms"
        onClick={() => setConsentOpen(true)}
        tooltip="The terms of service and cookie consent"
        icon={<TermsIcon/>}
      />
    </MenuBarMenu>
    <div className={styles.spacer}></div>
    <Consent open={consentOpen} onAccept={handleConsentAccept}/>
  </MenuBar>
})

export default MainMenu
