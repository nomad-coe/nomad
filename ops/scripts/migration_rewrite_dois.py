from nomad import doi, infrastructure, utils
from nomad.datamodel import Dataset

if __name__ == '__main__':
    infrastructure.setup_logging()
    infrastructure.setup_mongo()

    for dataset in Dataset.m_def.m_x('me').objects(doi__exists=True):
        try:
            doi.edit_url(doi=dataset.doi)
        except Exception as e:
            utils.get_logger('__name__').error('could not rewrite doi', exc_info=e)
        else:
            print('Rewrote URL of %s' % dataset.doi)
