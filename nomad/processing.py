from celery import Celery, chain, chord, group
import nomad.config as config
import nomad.files as files

broker_url = 'pyamqp://%s:%s@localhost//' % (config.rabbitmq.user, config.rabbitmq.password)
backend_url = 'rpc://localhost'
app = Celery('nomad.processing', backend=backend_url, broker=broker_url)


@app.task()
def process(upload):
  mainfiles = [('a', 'pa'), ('b', 'pb'), ('c', 'pc')]
  parsers = group([parse.s(mainfile, parser) for mainfile, parser in mainfiles]).delay()
  return parsers

@app.task()
def parse(mainfile, parser):
  return 'parsed %s with %s' % (mainfile, parser)

if __name__ == '__main__':
  print(~process.s('test'))
