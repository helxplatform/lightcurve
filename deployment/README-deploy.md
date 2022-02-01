# Sterling deployment

## Using helm

```
helm upgrade --install dev lightcurve --values lightcurve-values.yaml --dry-run

```

## Using kubectl

Deploy both the service and ingress.

Quick summary of the commands:
```
./k8-deploy-service.sh -n argus-dev

kubectl -n argus-dev apply -f lightcurve-ingress.yaml
```

-----

### k8s deployment of the service

Use the script which implements the sed substitutions of the secret S3 values: `k8-deploy-service.sh`

Alternatively, use one of the following manual methods:

1. In the `lightcurve-service.yaml` file set the values for `S3_ACCESS_KEY_ID` and `S3_SECRET_ACCESS_KEY` to allow access to the S3 data files.  These two values are for the **lightcurve-admin** user.

    ```
    kubectl -n argus-dev apply -f lightcurve-service.yaml

    # other
    kubectl -n argus-dev delete deployment lightcurve-deployment
    ```

2. Or use the **sed** utility to transform the yaml file with the required values.

    ```
    sed -e "s/S3_ACCESS_KEY_ID-value/$S3_ACCESS_KEY_ID/g" -e "s/S3_SECRET_ACCESS_KEY-value/$S3_SECRET_ACCESS_KEY/g" lightcurve-service.yaml | kubectl apply -f -
    ```

### k8s ingress creation

Run the kubectl command to apply the ingress.

```
kubectl -n argus-dev apply -f lightcurve-ingress.yaml

# other
kubectl -n argus-dev describe ingress lightcurveservice-ingress
kubectl -n argus-dev delete ingress lightcurveservice-ingress
```
