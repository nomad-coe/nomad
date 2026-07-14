# AGENTS

## Datetime handling

- Prefer `nomad.common.now()` for current timestamps (returns UTC and can be mocked in tests) instead of calling `datetime.now(...)` directly.
- MongoDB does not natively store timezone information. Treat stored datetime values as UTC+0, and re-attach timezone information in the ORM layer. Use `nomad.mongo.fields.UTCDateTimeField` for MongoEngine models so UTC timezone info is consistently restored on reads/writes.
- API responses should be RFC3339 compliant; use the Pydantic field `nomad.models.common.UTCDateTime`.
- Avoid manual timezone handling in application code when ORM/API field types above can enforce it.

## Storage Module

- The existing architecture adopts a two-tier storage approach: the main storage uses a local filesystem; the auxiliary
  storage uses an optional remote filesystem.
- The auxiliary storage `NOMADFileSystem` shall be used as a proxy sitting between user code and the underlying
  filesystem, interactions shall be done via `NOMADFileSystem.target_fs` which may fall back to the local filesystem if
  the auxiliary storage is not configured.
- User code shall be filesystem-agnostic, and should not directly interact with files via for example `os` module.
  Interactions shall be done via `fsspec` abstraction layer.
- `files.FSUtility` provides a set of utility functions for opening files in binary mode, as H5 files, as zip/tar
  archives, etc. Additional utility functions can be added as needed.
- User code (for example, `files.PathObject`) shall handle nominal paths. Absolutely avoid mixing nominal paths and real
  paths in user code.
- The `NOMADFileSystem.real_destination` can be used to resolve the actual path on the corresponding filesystem.