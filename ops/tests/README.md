# Run stress tests

You can run locust in a docker container like this.

```
docker run --rm -p 8089:8089 -v `pwd`:/home -ti locustio/locust -f loadtest_search.py --host https://nomad-lab.eu
```

This will open a HTTP web-interface on 8089. This can be used to control the test.
