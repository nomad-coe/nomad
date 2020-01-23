# increases the memory limits for flannel ds which causes alot of oom kills and
# further problems down the line
kubectl patch ds -n=kube-system kube-flannel-ds-amd64 -p '{"spec": {"template":{"spec":{"containers": [{"name":"kube-flannel", "resources": {"limits": {"cpu": "250m","memory": "550Mi"},"requests": {"cpu": "100m","memory": "100Mi"}}}]}}}}'