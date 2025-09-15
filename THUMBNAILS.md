Remote thumbnails workflow

- Export thumbnails from QuPath using `qupath/thumbnail_export.groovy`. It writes PNG files named `DONOR_ISLET.png` and a local `manifest.csv` with columns: `donor_id,islet_id,filename`.
- Upload the PNGs to a static hosting location (e.g., Amazon S3, CloudFront, Google Cloud Storage, Azure Blob, GitHub raw) that serves files over HTTPS.
- Create a public manifest CSV for the app with columns: `donor_id,islet_id,image_url`. Build `image_url` by concatenating your hosting base URL and the filename.

Example S3 mapping

- Upload `*.png` to `s3://my-bucket/panc_codex_thumbs/`
- Public URL format: `https://my-bucket.s3.amazonaws.com/panc_codex_thumbs/{filename}`
- Generate manifest (bash):

  awk -F, 'NR==1{next} {print $1","$2",https://my-bucket.s3.amazonaws.com/panc_codex_thumbs/"$3}' manifest.csv > manifest_public.csv

Use in the app

- In the Image Browser, choose "Manifest CSV (remote)" and upload `manifest_public.csv`.
- Select donor and islet to view remote thumbnails directly via URLs.

Tips

- Keep donor and islet IDs consistent between your Excel and the manifest.
- If QuPath doesn’t have `islet_id` in measurements, the script falls back to object name; you can set object names with a simple script or manually.
- Prefer 512–1024 px thumbnails to balance speed and visibility.

OneDrive notes

- Create a share link with Anyone-with-link Can View permission for each image or the parent folder.
- The app makes a best-effort to coerce OneDrive/SharePoint links by appending `?download=1` to the URL. For reliability, prefer links that already render the image directly.
- For OneDrive personal (`1drv.ms` or `onedrive.live.com`) and SharePoint (`sharepoint.com`) public links, adding `?download=1` often returns direct content suitable for `st.image`.
- You can generate a manifest by copying the public share URL to the `image_url` column.

QuPath overlays & channels

- The script `qupath/thumbnail_export.groovy` uses the current viewer’s channel visibility and colors. Before running:
  - Turn on the relevant channels: DAPI, INS, GCG, SST, CD3e, CD68, CD163, PDPN, CD31, B3TUBB, LGALS3.
  - Ensure annotations/detections are visible.
- The export draws outlines for objects with classes: `Islet`, `IsletExpanded`, `Nerve`, `Capillary`, `Lymphatic` over the thumbnail.
- `DOWNSAMPLE` controls output scale (8–16 recommended). Increase for smaller, faster thumbnails.
