layout: page
title: "GCP"
permalink: /GCP

## <a name='TOC'>Google Cloud Platform</a>

1. [Delete files by date](#Deletefilesbydate)


**<a name="Deletefilesbydate">Delete files by date</a>**
```python
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
