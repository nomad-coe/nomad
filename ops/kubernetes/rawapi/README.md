Deployment requires a secret to log into gitlab's docker registry:

```
kubectl create secret docker-registry gitlab-mpcdf --docker-server=gitlab-registry.mpcdf.mpg.de --docker-username=<your-user-name > --docker-password=<yourpass> --docker-email=<email>
```