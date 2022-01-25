if [ -z "${S3_ACCESS_KEY_ID}" ]; then
    echo S3_ACCESS_KEY_ID not set!
elif [ -z "${S3_SECRET_ACCESS_KEY}" ]; then
    echo S3_SECRET_ACCESS_KEY not set!
else
    if [ "$#" -lt 2 ]; then
        echo "-n namespace" required!
    else
        #set -x
        sed -e "s/S3_ACCESS_KEY_ID-value/$S3_ACCESS_KEY_ID/g" \
            -e "s/S3_SECRET_ACCESS_KEY-value/$S3_SECRET_ACCESS_KEY/g" \
            lightcurve-service.yaml | kubectl $@ apply -f -
        #set +x
    fi
fi
