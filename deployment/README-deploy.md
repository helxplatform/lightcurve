# Sterling deployment

## Using helm

```
# install
helm upgrade --install dev lightcurve --values lightcurve-values.yaml

# check values
helm upgrade --install dev lightcurve --values lightcurve-values.yaml --dry-run
```
