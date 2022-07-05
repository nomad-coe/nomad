## Cluster deployment, using Kubernetes and Helm

We use helm charts to describe the deployment of nomad services in a kubernetes cluster.
The NOMAD chart is part of the
[NOMAD source code](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR)
and can be found under `ops/helm/nomad`.

This chart allows to run the nomad app, worker, gui, and proxy in a kubernetes cluster.
The `values.yaml` contains more documentation on the different values.

The chart can be used to run multiple nomad instances in parallel on the same cluster,
by using different URL-path and database names.

The chart does not run any databases and search engines. Those are supposed to run
separately (see also *nomad-full* for an alternative approach) and their hosts, etc.
can be configures via helm values.
