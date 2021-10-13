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
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import Markdown from './Markdown'
import { email } from '../config'

class FAQ extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired
  }

  static styles = theme => ({
    root: {
      padding: theme.spacing(3)
    },
    container: {
      maxWidth: 1024,
      margin: 'auto',
      width: '100%'
    }
  })

  render() {
    const { classes } = this.props

    return (
      <div className={classes.root}><div className={classes.container}>
        <Markdown >{`
          # Frequently Asked Questions (FAQ)

          These are often repeated questions that cover the basic NOMAD use-cases. If you have
          further questions, please write use: [${email}](mailto:${email}).

          ## General

          ### How can I be sure that my data will be cited properly?

          Sharing means a change of culture. Making data open access is comparable to a
          publication where references to other work are common practice ever since.
          Likewise, using someone's data requires proper citation. We recommend the uploader
          to provide references to their data (publications, websites) and the users to
          cite these references together with the link to the repository.

          ### Can I change my email address or affiliation

          The NOMAD user management allows you to edit (click on your name after login)
          your name, email address, and affiliation. The email address and affiliation
          associated with your data might not change right away, but it will be updated
          eventually.

          ### I'd like to install NOMAD on my local computers to be only used by my group

          Local NOMAD deployments are not actively supported at the moment. However, all
          [NOMAD software](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR) is
          available under the Apache 2.0 Open-Source License and could in principle be
          installed and operated by you. We are planing to actively support local NOMAD
          installations and potential (partial) mirrors of NOMAD data under the name
          *NOMAD Oasis* in the future.

          ## Upload data, datasets, embargo, and DOIs

          ### What steps are necessary to publish data and reference it in a paper or share it otherwise?

          1. *upload* your data
          2. *publish* your data
          3. create a *dataset* and assign a *DOI*
          4. cite your data using a DOI or share URLs for datasets or individual entries

          ### Do I need to prepare my data before I upload?

          You need to create a file archive that contains all the input and
          output files you want to provide. You can either create a .zip (recommended)
          or .tar(.gz) file. Ideally, the input and output files of different code runs
          are in separate sub-directories. Otherwise, you can use any directory structure
          you want and also include additional files for each code run, including notes
          or other related artifacts.

          ### Only some of my files are shown by NOMAD, where is the rest?

          NOMAD scans the whole upload for files that it recognizes as main code output
          files (*mainfile*) of the supported codes. For each recognized *mainfile* NOMAD will
          create and show exactly one entry. NOMAD considers all files in the same directory
          of the *mainfile* to be *auxiliary* files of this entry. These files can include
          the input files, further output files, or other files in the same directory. Only
          *mainfiles* and *auxiliary* files can be downloaded publicly.

          If you know a file to be the main output file of a supported code, and it is
          still not recognized by NOMAD, let us know: [${email}](mailto:${email}).

          ### Some of my data is marked as *not processed* or *unavailable*, what does this mean?

          Codes and their input and output file formats evolve and NOMAD's parsers might
          not always be able to successfully parse all files. In these scenarios certain
          entry quantities are marked as *not processed* indicating that the entry could
          not be processed. This state might change in the future automatically, when
          we improve parsers and reprocess data.

          Some quantities only make sense in the context of certain codes or systems.
          For example, NOMAD cannot classify an entry's symmetry, if only a single
          molecule was simulated. In such cases, quantities are marked *unavailable*. These
          values might also change in the future as algorithms change or improve.

          ### Is my data visible to others

          Ultimately. NOMAD is a platform to share data, but we recognize that you might not want to
          share all data right away. There are two mechanisms to hide data for two
          different underlying motivations.

          First, not yet *published* data has restricted visibility. Once uploaded, data
          is only shown to you. Only after you publish your upload, it is visible to others.
          This allows you to delete wronly uploaded or processed data privately before
          publishing anything.

          Second, you can publish your uploads with an *embargo* period. This can last up to
          3 years. You can lift the embargo at anytime. Embargoed data is
          visible to and findable by others. This makes only some few metadata (e.g.
          chemical formula, system type, spacegroup, etc.) public, but the raw-file
          and archive contents remain hidden (except to you, and users you explicitly
          share the data with).
          You can already create datasets and assign DOIs for data with embargo, e.g.
          to put it into your unpublished paper.
          The embargo will last up to 36 month. Afterwards, your data will be made publicly
          available. You can also lift the embargo on entries at any time.

          ### How do I cite uploaded data in a paper?

          NOMAD allows to create *datasets* to curate data. A *dataset* is just a named
          collection of entries. A *dataset* does not have to correspond to an *upload*;
          you can freely create datasets. Each entry can belong to multiple datasets
          and datasets can overlap and form hierarchies.
          Datasets are created by editing entries and assigning a new dataset name to
          them. Entries can be edited from the [Upload](./uploads) and [Your data](./userdata)
          pages.

          NOMAD allows you to assign a DOI (Digital Object Identifier) to datasets. A DOI
          is a unique string that can be resolved to a URL that brings your readers to
          NOMAD where they can inspect and download your data.
          A DOI can be assigned after the dataset was created via the **dataset** tab on
          the [Your data](./userdata) page.

          ### Can I change datasets after assigning a DOI

          A DOI is a persistent identifier that is supposed to always reference the same
          data. Therefore, you cannot remove data (or whole datasets) after a DOI has
          been created. We only assign DOIs to datasets that exclusively contain published
          (but potentially embargoed) data.

          ### How can I share data that has an embargo?

          When you edit entries, you can mark other NOMAD users as *shared with*. Those
          users will see you embargoed data like any other data in NOMAD. For example,
          you could create a dataset, edit all entries in the dataset and set certain
          users as *shared with*. Take the URL from the dataset page and send it to the
          other users via email.

          ### How can I share credit with my co authors?

          When you edit the metadata, you can mark other NOMAD users as *co-authors* of
          an upload or an entry. The main author and the co-authors will comprise the
          respective authors list displayed for each entry.

          ### I want to upload data from a code that is not yet supported?

          If you are familiar with the input and output format of other relevant codes,
          write us an Email ([${email}](mailto:${email})) and we will figure out if and how
          to support this code in the future.
        `}</Markdown>
      </div></div>
    )
  }
}

export default withStyles(FAQ.styles)(FAQ)
