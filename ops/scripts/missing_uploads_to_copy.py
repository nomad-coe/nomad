import json

if __name__ == '__main__':
    with open('local/missing_calcs_data.json') as f:
        data = json.load(f)

    uploads = {}
    for key in ['no_package', 'no_calcs', 'missing_mainfile']:
        for item in data[key]:
            uploads[item['source_upload_id']] = 'yes'

    with open('local/missing_uploads', 'wt') as f:
        for upload in uploads.keys():
            f.write('%s\n' % upload)
