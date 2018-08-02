from collections import namedtuple

S3Config = namedtuple('S3', ['uploads_bucket', 'repository_bucket', 'archive_bucket'])
RabitMQConfig = namedtuple('RabbitMQ', ['host', 'port', 'user', 'password'])
MinioConfig = namedtuple('Minio', ['host', 'port', 'accesskey', 'secret'])
FSConfig = namedtuple('FSConfig', ['tmp'])

s3 = S3Config(
  uploads_bucket='uploads',
  repository_bucket='repository',
  archive_bucket='archive'
)
rabbitmq = RabitMQConfig(
  host = 'localhost',
  port = None,
  user = 'rabbitmq',
  password = 'rabbitmq'
)
minio = MinioConfig(
  host = 'localhost',
  port = 9007,
  accesskey = 'AKIAIOSFODNN7EXAMPLE',
  secret = 'wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY'
)
fs = FSConfig(
  tmp = './infrastructure/data/tmp'
)
