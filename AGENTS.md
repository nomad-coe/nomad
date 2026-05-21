# AGENTS

## Datetime handling

- Prefer `nomad.common.now()` for current timestamps (returns UTC and can be mocked in tests) instead of calling `datetime.now(...)` directly.
- MongoDB does not natively store timezone information. Treat stored datetime values as UTC+0, and re-attach timezone information in the ORM layer. Use `nomad.mongo.fields.UTCDateTimeField` for MongoEngine models so UTC timezone info is consistently restored on reads/writes.
- API responses should be RFC3339 compliant; use the Pydantic field `nomad.models.common.UTCDateTime`.
- Avoid manual timezone handling in application code when ORM/API field types above can enforce it.
