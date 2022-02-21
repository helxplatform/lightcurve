# Sterling deployment

## Using helm

```
# install
helm --namespace argus-dev upgrade --install evryscope lightcurve --values lightcurve-values.yaml

# check values
helm upgrade --install argus-dev lightcurve --values lightcurve-values.yaml --dry-run
```
