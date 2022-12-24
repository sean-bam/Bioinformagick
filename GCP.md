layout: page
title: "GCP"
permalink: /GCP

**Delete files by date**
```
from google.cloud import storage
from datetime import datetime, timezone

client = storage.Client()
bucket = client.get_bucket("bucket")
blobs = bucket.list_blobs(prefix="path/to/folder", delimiter="/")

today_time = datetime.now(tz=timezone.utc)

for blob in blobs:
    blob_time = blob.time_created
    delta = today_time - blob_time
    if delta.days > 1:
        blob.delete()
```
