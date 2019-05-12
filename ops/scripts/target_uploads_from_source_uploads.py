import sys
from nomad import infrastructure
from nomad.migration import Package
from nomad.processing import Upload

infrastructure.setup_logging()
infrastructure.setup_mongo()

sources = sys.argv[1:]
names = list(package.package_id for package in Package.objects(upload_id__in=sources))
targets = Upload.objects(name__in=names)
for target in targets:
    print(targets.upload_id)
