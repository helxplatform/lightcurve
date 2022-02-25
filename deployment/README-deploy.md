# Sterling deployment

## Using helm

```
# install
helm --namespace argus-dev upgrade --install evryscope lightcurve --values lightcurve-values.yaml

# check values
helm --namespace argus-dev upgrade --install evryscope lightcurve --values lightcurve-values.yaml --dry-run
```

## Checking S3 Secrets
```
kubectl --namespace argus-dev get secret lightcurve-s3 --template={{.data.S3_ACCESS_KEY_ID}} | base64 -d

kubectl --namespace argus-dev get secret lightcurve-s3 --template={{.data.S3_SECRET_ACCESS_KEY}} | base64 -d
```
