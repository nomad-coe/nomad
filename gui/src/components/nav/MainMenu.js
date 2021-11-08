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
import { useRoute } from './Routes'
import { MenuBar, MenuBarItem, MenuBarMenu } from './MenuBar'
import { appBase, oasis, aitoolkitEnabled, encyclopediaEnabled } from '../../config'

import BackupIcon from '@material-ui/icons/Backup'
import SearchIcon from '@material-ui/icons/Search'
import UserDataIcon from '@material-ui/icons/AccountCircle'
import AboutIcon from '@material-ui/icons/Home'
import ForumIcon from '@material-ui/icons/QuestionAnswer'
import FAQIcon from '@material-ui/icons/LiveHelp'
import EncyclopediaIcon from '@material-ui/icons/Language'
import MetainfoIcon from '@material-ui/icons/Info'
import DocIcon from '@material-ui/icons/Help'
import CodeIcon from '@material-ui/icons/Code'
import TermsIcon from '@material-ui/icons/Assignment'
import AnalyticsIcon from '@material-ui/icons/ShowChart'

export default function MainMenu() {
  const route = useRoute()
  const selected = (route?.navPath) || 'publish/uploads'

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
        name="search" route="/search"
        tooltip="Find and download data"
      />
      {encyclopediaEnabled && <MenuBarItem
        name="encyclopedia"
        href={`${appBase}/encyclopedia/#/search`}
        tooltip="Visit the NOMAD Materials Encyclopedia"
        icon={<EncyclopediaIcon/>}
      />}
    </MenuBarMenu>
    <MenuBarMenu name="analyze" route="/metainfo" icon={<AnalyticsIcon/>}>
      {(!oasis && aitoolkitEnabled)
        ? <MenuBarItem
          label="AI Toolkit" name="aitoolkit" route="/aitoolkit/main"
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
        name="metainfo" route="/metainfo" tooltip="Browse the NOMAD Archive schema"
      />
      <MenuBarItem
        name="apis" label="APIs" route="/apis" tooltip="The list of APIs offered by NOMAD"
      />
    </MenuBarMenu>
    <MenuBarMenu name="about" route="/" icon={<AboutIcon/>}>
      <MenuBarItem
        label="Information" name="info" route="/"
        tooltip="About the NOMAD Repository and Archive"
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
        name="terms"
        href="https://nomad-lab.eu/terms"
        tooltip="The terms of service"
        icon={<TermsIcon/>}
      />
    </MenuBarMenu>
  </MenuBar>
}
