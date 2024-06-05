# How to release a new NOMAD version

## What is a release

NOMAD is a public service, a Git repository, a Python package, and a docker image.
What exactly is a NOMAD release? It is all of the following:

- a version tag on the main NOMAD [git project](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR),
e.g. [`v1.3.0`](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/-/tags/v1.3.0)
- a gitlab release based on a tag with potential release notes
- a version of the `nomad-lab` Python package released to pypi.org, e.g. `nomad-lab==1.3.0`.
- a docker image tag, e.g. `gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-fair:v1.3.0`
- the docker image tag `stable` points to the image with the latest release tag

## Steps towards a new release

- Silently create a new version tag in the `v1.3.0` format.
- Deploy the build from this tag to the public NOMAD deployments.
What deployments are updated might depend on the current needs. But usually
the production and test deployment should be updated.
- Release the Python package to the local gitlab registry. (This will update the
NORTH Jupyter image in the next nightly build and most likely effect plugins)
- Bump the `latest` docker image tag.
- For minor and major releases, encourage (Oasis) users to test the public services and the latest docker image for a short trial phase (e.g. 3 days). For patch releases this step should be
skipped.
- Create a gitlab release from the tag with potential release notes. Those notes
should also be added to the README.md. It is ok, if the updated README.md is not part of the
release itself.
- Bump the `stable` docker image tag.
- Publish the Python package to [pypi.org](https://pypi.org/)

## How to deal with hotfixes

This depends on the current `develop` branch and requires a judgement call. There are
two opposing scenarios:

1. The `develop` branch only contains minor fixes or fix/features that are not likely to effect
the released functionality. In this case, a new release with an increased patch version
is the right call.

2. The `develop` branch adds major refactorings and commits that likely effect the
released functionality. In this case, a `v1.3.0-hotfix` branch should be created.
After adding commits with the hotfix, the release process can be applied to the
hotfix branch in order to create a `v1.3.1` release that only contains the hotfixes and
not the changes on develop. After the `v1.3.1` release, the `v1.3.0-hotfix` branch is merged
back into develop. Hotfix branches should not live longer than a week.

## Major, minor, patch versions

- **patch**: No significant refactorings. Only new/updated features behind disabled feature switches.
Bugfixes. Might mark features as deprecated.

- **minor**: Might enabled new features by default. Can contain major refactorings (especially if they effect to plugin developers, data stewards etc.). Might finally deprecate features.
Should "basically" be backwards compatible.

- **major**: Breaking changes and will require data migration.

What is a *breaking change* and what does "basically" backwards compatible mean?
We develop experimental functionality and often need multiple iterations
to get a feature right. This also means that we technically introduce breaking changes
far more often than we can issue major releases. It is again a judgement call to decide on
major vs minor. The following things would generally not be considered *breaking* and would be considered *backwards compatible*:

- the breaking change is for a feature that is not enabled by default
- data migration is necessary for new functionality, but optional for existing functionality
- it is unlikely that plugins not developed by FAIRmat are effected
- it is unlikely that data beyond the central NOMAD deployments need to be migrated

## Release schedule

Patch releases should happen frequently and at least once every other month. Also
minor releases should be done semi regular. Important new features or at least
bi-annual FAIRmat events should trigger a minor release. Major releases require
more involved planning, data migration, and respective instructions and assistance to
NOMAD (Oasis) users. They are also political. Therefore, they do not a have a regular
schedule.

With a one `develop` branch Git strategy, there might be necessary exceptions to
regular patch releases. In general, new features should be protected by feature switches,
and should not be an issue. However, major refactorings that might effect multiple components are hard to hide behind a feature switch. In such cases, the release schedule might be
put on hold for another month or two.