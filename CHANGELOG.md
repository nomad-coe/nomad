## v1.1.6 (2022-12-23)

### fixed (1 change)

- [Fixed failing plot annotation valiation.](nomad-lab/nomad-FAIR@3cd5d0635566d9882e1c44cd922e113dfbeb1292) ([merge request](nomad-lab/nomad-FAIR!1005))

### feature (1 change)

- [Added a way to define formal models for metainfo annotations.](nomad-lab/nomad-FAIR@6a84b55e17dd5cd41cb5ba571f2451d528dbe6a0) ([merge request](nomad-lab/nomad-FAIR!990))

### added (1 change)

- [Replaced config object with pydantic models.](nomad-lab/nomad-FAIR@b9c1e023d95d14fd337b98a7ac895f505d4715a9) ([merge request](nomad-lab/nomad-FAIR!981))
## 1.1.7 (2023-02-20)

### Added (17 changes)

- [Added a config switch to enable disable north tools.](nomad-lab/nomad-FAIR@37d9b6db1e8044ae28476488e8c3ae1896e149ad) ([merge request](nomad-lab/nomad-FAIR!1099))
- [Refactor NORTH and Jupyterhub deployment. Add a new north API endpoint. Mount...](nomad-lab/nomad-FAIR@b0e28fa7fd3f47fcd538f33bd84bad659ee6549a) ([merge request](nomad-lab/nomad-FAIR!996))
- [Added categories to create entry build-in schemas.](nomad-lab/nomad-FAIR@6a2b96f58c29a09364523c1d140c13c6a9cf24b0) ([merge request](nomad-lab/nomad-FAIR!1085))
- [Added support for slice notation when choosing axis quantities in the plot annotation.](nomad-lab/nomad-FAIR@1fedce5e2a310fde2d304aa38262cbc79adeb289) ([merge request](nomad-lab/nomad-FAIR!1073))
- [New workflow visualizer](nomad-lab/nomad-FAIR@e02bfb3cc5fdcc5eb3beda6910a204e84fde7f5f) ([merge request](nomad-lab/nomad-FAIR!957))
- [Added pagination feature for long sections in the archive browser.](nomad-lab/nomad-FAIR@b13839fdaa1bdd111b2cd32e35a7887ce65db816) ([merge request](nomad-lab/nomad-FAIR!1081))
- [Added categories to the example uploads dialog.](nomad-lab/nomad-FAIR@994e272b478682f25fbdc5b8748ef09764a3ac8f) ([merge request](nomad-lab/nomad-FAIR!1080))
- [Added a new 3D structure viewer based on NGL.](nomad-lab/nomad-FAIR@91e8cd67247b5881b842df14014929098abce0df) ([merge request](nomad-lab/nomad-FAIR!776))
- [Dialog to create entry](nomad-lab/nomad-FAIR@5602071ffb659700360983e2ca7b27105de391b9) ([merge request](nomad-lab/nomad-FAIR!1032))
- [Added ELN Base Sections](nomad-lab/nomad-FAIR@4895124d43de7f3fe57b5c97281e3e86676cab4a) ([merge request](nomad-lab/nomad-FAIR!1037))
- [Ability to create and modify child classes from referenceEditQuantity](nomad-lab/nomad-FAIR@70ec0c443ea2cbf2481c162215283981b714b0b5) ([merge request](nomad-lab/nomad-FAIR!1019))
- [Added possibility for changing the page title to distinguish OASIS installations.](nomad-lab/nomad-FAIR@3604c83813251df4ec9c7daee021ba98cb7baca3) ([merge request](nomad-lab/nomad-FAIR!1056))
- [Added automatic push of develop to github.](nomad-lab/nomad-FAIR@e806087e77f8c6bcbeea9e42c218d0b377be3c6e) ([merge request](nomad-lab/nomad-FAIR!1047))
- [Added a to_ase function to simulation system Atoms section.](nomad-lab/nomad-FAIR@f97425cc835a08f4809a676af246c8ad8e2ea54f) ([merge request](nomad-lab/nomad-FAIR!1027))
- [Added a automatically generated changelog to our git.](nomad-lab/nomad-FAIR@fe795eff97145e0e28d61ca840cbb3069217d24a) ([merge request](nomad-lab/nomad-FAIR!1016))
- [Added a glossary to the documentation.](nomad-lab/nomad-FAIR@aa42834ed51b3d231a8edc58b00ea1f2dddc21f5) ([merge request](nomad-lab/nomad-FAIR!1024))
- [Added dropdown list to select from sections.](nomad-lab/nomad-FAIR@a332abed1318f60078463e2ee71ae68b2086ef1a) ([merge request](nomad-lab/nomad-FAIR!992))

### Fixed (16 changes)

- [Fixed wrong property inheritance base section order in UI metainfo.](nomad-lab/nomad-FAIR@9cbd6d9c37d21f06cf9fc3f5d9fbde746d7873b3) ([merge request](nomad-lab/nomad-FAIR!1098))
- [Fixed the build-in schema entry create to use the right definitions.](nomad-lab/nomad-FAIR@8db7bb1be8ab502f4aef178107f50d27297123a2) ([merge request](nomad-lab/nomad-FAIR!1095))
- [Fixed direct navigation of reference edit quantity by disabling it.](nomad-lab/nomad-FAIR@b6ba106989192407954d7846f13abfe13fde431d) ([merge request](nomad-lab/nomad-FAIR!1095))
- [Fixed bi-directional referencing in ELNs.](nomad-lab/nomad-FAIR@b5f8c7ac5f120943f5c235e3ee216424b26f0188) ([merge request](nomad-lab/nomad-FAIR!1095))
- [Resolve "quickfix: typo in scripts/build_sdist.sh"](nomad-lab/nomad-FAIR@0abcb2a53011529d628074203e4b6a67a5fbdf25) ([merge request](nomad-lab/nomad-FAIR!1090))
- [Fixed incorrect URL reported in the API call dialogs.](nomad-lab/nomad-FAIR@1f58dd2f3271c95829867933255bbb11c2540e61) ([merge request](nomad-lab/nomad-FAIR!1089))
- [Resolve the reference to a data defined together the schema in the same file](nomad-lab/nomad-FAIR@94a58aeaf49970b8bc9f8b1ef9605b790e125d31) ([merge request](nomad-lab/nomad-FAIR!1079))
- [ability to reference all with m_def not only data](nomad-lab/nomad-FAIR@d56a3645ed3515a2e63231d93ca013b4b6be3316) ([merge request](nomad-lab/nomad-FAIR!1079))
- [only fetch those entries which are processed](nomad-lab/nomad-FAIR@029e5ead244382330856fc541c4ad3c9d5d2e3c2) ([merge request](nomad-lab/nomad-FAIR!1079))
- [normalize path in url](nomad-lab/nomad-FAIR@d6fbcec4bc58b20910ea283dd3c6c1b741866b9d) ([merge request](nomad-lab/nomad-FAIR!1079))
- [URLEditQuantity is added as a string eln_component to annotations](nomad-lab/nomad-FAIR@a5a916b32358bfa4648ce7899e5dd0cf59186a40) ([merge request](nomad-lab/nomad-FAIR!1062))
- [Fixed issues with the developers setup of nomad.](nomad-lab/nomad-FAIR@d76fc7e061fb6a1662cfa13da6519dde51345ffa) ([merge request](nomad-lab/nomad-FAIR!1064))
- [Fixed endless recursion in metainfo browser.](nomad-lab/nomad-FAIR@1443063a8b3a8df36975204c12d23ba51b9962e2) ([merge request](nomad-lab/nomad-FAIR!1059))
- [Fixed the README batches and changed pipeline to push tags to github.](nomad-lab/nomad-FAIR@8c42e102387b2d0c5f0a55fd368b266306b705f1) ([merge request](nomad-lab/nomad-FAIR!1057))
- [Fixed issue with optimade filter menu performing premature API calls.](nomad-lab/nomad-FAIR@7c4ad789591ce6a30cb6236825798e440711ce21) ([merge request](nomad-lab/nomad-FAIR!1048))
- [Fixed generation of duplicated subsections on updating an entry that triggers tabular parser](nomad-lab/nomad-FAIR@e91176fe2a197f9d8450e29dc8b599f669ecad1d) ([merge request](nomad-lab/nomad-FAIR!1033))

### Changed (2 changes)

- [Small UX changes to the upload interface.](nomad-lab/nomad-FAIR@5e18b0e12d4ce4d608a5cd68efbeb02ee2d51f04) ([merge request](nomad-lab/nomad-FAIR!1085))
- [Refactored the overview UI code to optimise loading performance.](nomad-lab/nomad-FAIR@161051655e8b9406197247ed943e027fdd13771b) ([merge request](nomad-lab/nomad-FAIR!1075))

### fixed (2 changes)

- [Fix the select feature of the Datatable component](nomad-lab/nomad-FAIR@8fd52ed7bccc14c13c6c83d90e3a931df3cf9d25) ([merge request](nomad-lab/nomad-FAIR!1060))
- [Datafile quantity containing `tabular_parser` annotation can be referenced](nomad-lab/nomad-FAIR@8ad191bf7f37b9c1286200127f6f5b2f0cf3701a) ([merge request](nomad-lab/nomad-FAIR!1017))

### added (1 change)

- [Implement extended matching of h5 files.](nomad-lab/nomad-FAIR@96372aa26a0675a0bbae99b8e5059f5e0f1a5f3a) ([merge request](nomad-lab/nomad-FAIR!1023))

### feature (1 change)

- [Adds an ELN schema that allows to create entries that clone data from a labfolder project.](nomad-lab/nomad-FAIR@166a8aaa65f63bf368453365689e1dab577fc94d) ([merge request](nomad-lab/nomad-FAIR!948))
