import csv
import sys

with open(sys.argv[1]) as f:
    reader = csv.reader(f, delimiter=',')
    data = [row for row in reader]

for row in data[1:]:
    if len(row) == 8:
        message, n_entries, entry_id, upload_id, n_uploads, mainfile, parser, exc  = row[0:8]
    if len(row) == 15:
        message, n_entries, entry_id, upload_id, mainfile, _, _, _, _, _, parser = row[0:11]
        n_uploads = 'n.n.'
        exc = None
        if 'process failed' in message:
            continue
    else:
        assert False, "Wrong csv format."
    
    print(f'''
- [ ] ({n_entries} entries/{n_uploads} uploads) *{parser}: {message}*

example: [upload_id: {upload_id}, entry_id: {entry_id}, mainfile: {mainfile}](https://nomad-lab.eu/prod/v1/gui/entry/entry_id/{entry_id})

{"```exc```" if exc else ""}


''')
