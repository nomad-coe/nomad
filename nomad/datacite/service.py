from typing import TYPE_CHECKING

from nomad import datamodel
from nomad.datamodel.data import User
from nomad.datamodel.datamodel import Dataset

from .client import DataCiteClient
from .models import Affiliation, Creator, DoiRequestAttributes, Title, TypesBase
from .utils import generate_target_url, generate_unique_doi_name

if TYPE_CHECKING:
    from nomad.processing.data import Upload


def convert_user_id_to_creator(user_id: str) -> Creator:
    """Generates a DataCite Creator from a NOMAD user ID."""
    user = datamodel.User.get(user_id=user_id)
    return convert_user_to_creator(user)


def convert_user_to_creator(user: User) -> Creator:
    """Generates a DataCite Creator from a NOMAD User."""
    affiliation = Affiliation()
    if user.affiliation is not None:
        affiliation.name += user.affiliation.strip()
    if user.affiliation_address is not None:
        affiliation.name += '; ' + user.affiliation_address.strip()

    creator = Creator(name=user.name, affiliation=[affiliation])
    return creator


def generate_base_attributes() -> DoiRequestAttributes:
    """Generates base attributes for a DOI."""
    data = DoiRequestAttributes()
    data.publisher = 'NOMAD Repository'
    data.types = TypesBase()
    data.types.resourceTypeGeneral = 'Collection'
    data.types.resourceType = 'Materials Science Data'

    return data


def create_doi(attributes: DoiRequestAttributes) -> dict | None:
    client = DataCiteClient()
    if not client.enabled:
        return None

    # Using POST to trigger error if DOI alread exists
    response = client.post_doi(attributes)
    return response


def publish_doi(doi: str):
    client = DataCiteClient()
    if not client.enabled:
        return None

    attributes = DoiRequestAttributes(event='publish')
    response = client.put_doi(doi, attributes)
    return response


def delete_doi(doi: str):
    client = DataCiteClient()
    if not client.enabled:
        return None

    response = client.delete_doi(doi)
    return response


def create_attributes_from_dataset(dataset: Dataset) -> DoiRequestAttributes:
    attributes = generate_base_attributes()
    attributes.titles = [Title(title=dataset.dataset_name)]
    attributes.creators = [convert_user_id_to_creator(dataset.user_id)]
    attributes.publicationYear = dataset.dataset_create_time.year
    attributes.doi = dataset.doi

    return attributes


def create_attributes_from_upload(upload: 'Upload') -> DoiRequestAttributes:
    attributes = generate_base_attributes()
    attributes.titles = [Title(title=upload.upload_name)]
    attributes.creators = [convert_user_id_to_creator(upload.main_author)]
    attributes.publicationYear = upload.upload_create_time.year

    return attributes


def create_doi_for_dataset(dataset: Dataset) -> str | None:
    if dataset.doi is not None:
        raise ValueError('Dataset already has a DOI.')

    attributes = create_attributes_from_dataset(dataset)
    attributes.doi = generate_unique_doi_name()
    attributes.url = generate_target_url(attributes.doi, 'dataset')
    create_doi(attributes)

    return attributes.doi


def create_doi_for_upload(upload: 'Upload') -> str | None:
    if upload.doi is not None:
        raise ValueError('Upload already has a DOI.')

    attributes = create_attributes_from_upload(upload)
    attributes.doi = generate_unique_doi_name()
    attributes.url = generate_target_url(attributes.doi, 'upload')
    create_doi(attributes)

    return attributes.doi


def create_attributes_from_args(
    title: str, publicationYear: int, user: User
) -> DoiRequestAttributes:
    attributes = generate_base_attributes()
    attributes.titles = [Title(title=title)]
    attributes.creators = [convert_user_to_creator(user)]
    attributes.publicationYear = publicationYear
    attributes.doi = generate_unique_doi_name()
    attributes.url = generate_target_url(attributes.doi, 'dataset')

    return attributes


def create_draft_doi_from_args(title: str, year: int, user: User) -> str:
    attributes = create_attributes_from_args(title, year, user)
    create_doi(attributes)

    return attributes.doi


def generate_example_doi_attributes() -> DoiRequestAttributes:
    attributes = generate_base_attributes()
    attributes.titles = [Title(title='Example Dataset')]
    attributes.creators = [Creator(name='John Doe')]
    attributes.publicationYear = 2024
    attributes.doi = generate_unique_doi_name()
    attributes.url = 'https://example.com/dataset'

    return attributes
