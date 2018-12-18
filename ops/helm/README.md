## Cluster Deployment, Using Kubernetes and Helm

We use helm charts to describe the deployment of nomad services in a kubernetes cluster.

### nomad

This chart allows to run the nomad api, worker, gui, and proxy in a kubernetes cluster.
The `values.yaml` contains more documentation on the different values.

The chart can be used to run multiple nomad instances in parallel on the same cluster,
by using different URL-path and database names.

The chart does not run any databases and search engines. Those are supposed to run
separately (see also *nomad-full* for an alternative approach) and their hosts, etc.
can be configures via helm values.

### rawapi

Similar to *nomad* and similar to *rawapi* in `docker-compose`. Runs rawapi solo for
the *materials project*.

### nomad-full

This chart is under development. It is an attempt to also run all required databases
and search engine in the same kubernetes cluster.