'''
A simple example used in the NOMAD webinar API tutorial
'''

from nomad.client import ArchiveQuery, Auth

query = ArchiveQuery(
    owner='user',
    required={
        'run': {
            'system[-1]': '*'
        }
    },
    authentication=Auth(user='yourusername', password='yourpassword'))

print(query)
