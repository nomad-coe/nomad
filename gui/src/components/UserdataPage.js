import React from 'react'
import { withApi } from './api'
import Search from './search/Search'

export const help = `
This page allows you to **inspect** and **manage** you own data. It is similar to the
*search page*, but it will only show data that was uploaded by you or is shared with you.

Besides giving you a filter for your data, this page also allows you to edit the *user metadata*
on all entries. User metadata is assigned to individual entries, but of course you can
edit many entries at once. User metadata comprises comments, references, datasets, co-authors,
and users to share data with.

#### Selecting entries for edit

You can use the search in combination with the entries, dataset, or uploads table to
help select sets of entries that you like to edit. If you select individual entries,
only these entries will be edit. If you select all entries, all entries fitting your
search (including those
  on further pages) will be edited. If you hit edit on datasets or uploads all entries
  in these dataset or upload will be edited.

#### The edit user metadata dialog

If you click an edit button, a dialog will appear. It is pre-filled with the values of
the first selected entry. However, the settings will be applied to all selected entries,
potentially overwriting existing data. Only those settings that you change will be applied.
For example, if you only set a new comment, only the comments on all entries will be overwritten.
If you add a co-author, you changed the co-authors and all co-authors on all entries will be
set with the new authors.

#### Datasets

The edit dialog also allows you to put entries into datasets. Just enter the name of
an existing or new dataset. We only consider datasets with data. You should never see
an empty dataset. Entries can be put into many datasets. There is no explicit hierarchy
between datasets, but of course you can effectively form super-sets by assigning entries
respectively.

To remove entries from datasets, simple set the list of datasets for those entries to a
list that is empty or does not contain the respective dataset anymore. If you remove
all entries from a dataset like this, it will effectively disappear. You can also
delete datasets from the dataset table with the delete button. The dataset will be removed
from the dataset list of all its current entries.

#### DOIs

On the dataset table, you can also click the DOI button to assign a DOI to a dataset.
This DOI can be used in publications to link to your dataset. If people resovle the
DOI, they will be redirected to a NOMAD view that shows the dataset and allows its download.

Once you assigned a DOI to a dataset, no entries can be removed or added to the dataset.
`

function UserdataPage() {
  return <Search
    ownerTypes={['user', 'shared', 'staging']}
    initialOwner="user"
    initialRequest={{order_by: 'upload_time', uploads_grouped: true}}
    initialResultTab="uploads"
    availableResultTabs={['uploads', 'datasets', 'entries']}
    resultListProps={{selectedColumnsKey: 'userEntries', selectedColumns: ['formula', 'upload_time', 'mainfile', 'published', 'co_authors', 'references', 'datasets']}}
  />
}

export default withApi(true, false, 'To manage you data, you must log in.')(UserdataPage)
