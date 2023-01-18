---
title: Using Git/GitLab
---

## Branches

We use *protected* and *feature* branches. You must not (even if you have
the rights) commit to it directly.

- `master`, a *protected* branch and also the *default* branch. Represents the latest stable
release (i.e. what the current official NOMAD runs on).
- `develop`, a quasi *protected* branch. It contains the latest, potentially
unreleased, features and fixes.
- *feature branches*, this is where you work. Typically they are automatically
named after issues: `<issue-number>-<issue-title>`.
- `vX.X.X`, tags for releases.


## Flow: From issue to merge

### Create issues

Everyone with an MPCDF GitLab account can create issues.

- Use descriptive short titles. Specific words over general works. Try to describe
the problem and do not assume causes.
- You do not need to greet us or say thx and goodbye. Let's keep it purely technically.
- For bugs: think about how to re-produce the problem and provide examples when applicable.
- For features: in addition to the feature descriptive, try to provide a use-case that
might help us understand the feature and its scope.
- You can label the issue if you know the label system.

### Labeling issues

There are three main categories for labels. Ideally each issue gets one of each category:

- The *state* labels (grey). Those are used to manage when and how the issue is addressed.
This should be given the NOMAD team member who is currently responsible to moderate the
development. If its a simple fix and you want to do it yourself, assign yourself and use
"bugfixes".

- The *component* label (purple). This denote the part of the NOMAD software that is
most likely effected. Multiple purple labels are possible.

- The *kind* label (key). Wether this is a bug, feature, refactoring, or documentation issue.

Unlabeled issues will get labeled by the NOMAD team as soon as possible.


### Working on an issue

Within the development team and during development meetings we decide who is acting on
an issue. The person is assigned to the issue on GitLab. It typically switches its grey *state* label
from *backlog* to *current* or *bugfixes*. The assigned person is responsible to further
discuss the problem with the issue author and involve more people if necessary.

To contribute code changes, you should create a branch and merge request (MR) from the
GitLab issue via the button offered by GitLab. This will create a branch of *develop* and
create a merge request that targets *develop*.

You can work on this branch and push as often as you like. Each push will cause a pipeline
to run. The solution can be discussed in the merge request.


### Code review and merge

When you are satisfied with your solution, and your CI/CD pipeline passes you can mark your MR as *ready*.
Make sure that the `delete` option is checked to automatically remove you branch after
merge. In most cases you should also check `squash` commits. This will replace all the
commits in the MR with a single *squash* commit. If you do not want to *squash* your
branch, make sure that you produce a reasonable and [clean version history](#clean-version-history).

To review GUI changes, you should deploy your branch to the dev-cluster via CI/CD actions.
Find someone on the NOMAD developer team to review your MR and request a review through
GitLab. The review should be performed shortly and should not stall the MR longer than
two full work days.

The reviewer will open *threads* that need to be solved by the MR author. If all
threads are resolved, you can re-request a review. The reviewer should eventually merge
the MR. Typically we squash MRs to keep the revision history short.
This will typically auto-close the issue.

## Changelog

We have an automatically generated changelog in the repository file `CHANGELOG.md`.
This changelog is produced from commit messages and to maintain this file, you
need to write commit messages accordingly.

To trigger a changelog entry, your commit needs to end with a so called *git trailer*
called `Changelog`. A typical commit message for a changelog entry should look like this:

```
A brief one line title of the change.

A longer *markdown* formatted description of the change. Keep in mind that gitlab
will automatically link the changelog entry with this commit and a respective merge
requests. You do not need to manually link to any gitlab resources.

This could span multiple paragraphs. However, keep it short. Documentation should
go into the actual documentation, but you should mention breaks in backward compatibility,
deprecation of features, etc.

Changelog: Fixed
```

The trailer value (`Fixed` in the example) has to be one of the following values:

- `Fixed`, for bugfixes.
- `Added`, for new features.
- `Changed`, for general improvements, e.g. updated documentation, refactoring,
improving performance, etc.

These categories are consistent with (keepachangelog.com)[https://keepachangelog.com/].
For more information about the changelog generation read the [gitlab documentation](https://docs.gitlab.com/ee/api/repositories.html#add-changelog-data-to-a-changelog-file).


## Clean version history

It is often necessary to consider code history to reason about potential problems in
our code. This can only be done, if we keep a "clean" history.

- Use descriptive commit messages. Use simple verbs (*added*, *removed*, *refactored*, etc.)
name features and changed components. [Include issue numbers](https://docs.gitlab.com/ee/user/project/issues/crosslinking_issues.html)
to create links in gitlab.

- Learn how to amend to avoid lists of small related commits.

- Learn how to rebase. Only merging feature-branches should create merge commits.

- Squash commits when merging.

- Some videos on more advanced git usage: https://youtube.be/Uszj_k0DGsg, https://youtu.be/qsTthZi23VE

### amend
While working on a feature, there are certain practices that will help us to create
a clean history with coherent commits, where each commit stands on its own.

```sh
  git commit --amend
```

If you committed something to your own feature branch and then realize by CI that you have
some tiny error in it that you need to fix, try to amend this fix to the last commit.
This will avoid unnecessary tiny commits and foster more coherent single commits. With *amend*
you are basically adding changes to the last commit, i.e. editing the last commit. If
you push, you need to force it `git push origin feature-branch --force-with-lease`. So be careful, and
only use this on your own branches.

### rebase
```sh
  git rebase <version-branch>
```

Lets assume you work on a bigger feature that takes more time. You might want to merge
the version branch into your feature branch from time to time to get the recent changes.
In these cases, use rebase and not merge. Rebase puts your branch commits in front of the
merged commits instead of creating a new commit with two ancestors. It basically moves the
point where you initially branched away from the version branch to the current position in
the version branch. This will avoid merges, merge commits, and generally leave us with a
more consistent history.  You can also rebase before creating a merge request, which basically
allows no-op merges. Ideally the only real merges that we ever have, are between
version branches.


### squash
```sh
  git merge --squash <other-branch>
```

When you need multiple branches to implement a feature and merge between them, try to
use *squash*. Squashing basically puts all commits of the merged branch into a single commit.
It basically allows you to have many commits and then squash them into one. This is useful
if these commits were made just to synchronize between workstations, due to
unexpected errors in CI/CD, because you needed a save point, etc. Again the goal is to have
coherent commits, where each commits makes sense on its own.


## Submodules


The main NOMAD GitLab-project (`nomad-fair`) uses Git-submodules to maintain its
parsers and other dependencies. All these submodules are places in the `/dependencies`
directory. There are helper scripts to install (`./dependencies.sh`) and
commit changes to all submodules (`./dependencies-git.sh`). After merging or checking out,
you have to make sure that the modules are updated to not accidentally commit old
submodule commits again. Usually you do the following to check if you really have a
clean working directory.

```sh
  git checkout something-with-changes
  git submodule update --init --recursive
  git status
```

We typically use the `master`/`main` branch on all dependencies. Of course feature branches
can be used on dependencies to manage work in progress.
