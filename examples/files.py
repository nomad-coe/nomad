from nomad import infrastructure, files, processing as proc

infrastructure.setup_logging()
infrastructure.setup_mongo()

upload_id = 'NvVyk3gATxCJW6dWS4cRWw'
upload = proc.Upload.get(upload_id)
upload_with_metadata = upload.to_upload_with_metadata()

upload_files = files.PublicUploadFiles(upload_id)
upload_files.repack(upload_with_metadata)

# try:
#     public_upload_files = files.PublicUploadFiles(upload_id)
#     public_upload_files.delete()
# except Exception:
#     pass

# staging_upload_files = files.StagingUploadFiles(upload_id)
# staging_upload_files.pack(upload_with_metadata)
