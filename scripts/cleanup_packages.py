import os
import sys
import requests
import pandas as pd
from packaging import version


def get_env_var(var_name):
    value = os.getenv(var_name)
    if value is None:
        print(f'Missing {var_name} environment variable')
        sys.exit(1)
    return value


CI_JOB_TOKEN = get_env_var('CI_JOB_TOKEN')
CI_PROJECT_ID = get_env_var('CI_PROJECT_ID')

NB_PACKAGES_PER_PAGE = 100


def fetch_packages():
    packages = []
    hasNext = True
    nextLink = f'https://gitlab.mpcdf.mpg.de/api/v4/projects/{CI_PROJECT_ID}/packages?page=1&per_page={NB_PACKAGES_PER_PAGE}&order_by=created_at&sort=asc'
    job_token = {'JOB-TOKEN': CI_JOB_TOKEN}
    while hasNext and nextLink:
        response = requests.get(nextLink, headers=job_token)
        print(response)

        if response.status_code != 200:
            print('Unable to list Gitlab packages, no cleanup can be done')
            sys.exit(1)

        packages.extend(response.json())
        nextLink = None
        links = response.headers.get('Link', '')
        if 'rel="next"' in links:
            links_parts = links.split(',')
            for part in links_parts:
                if 'rel="next"' in part:
                    nextLink = part[part.find('<') + 1 : part.find('>')]
                    break
        hasNext = nextLink is not None

    return pd.DataFrame(packages)


def find_packages_to_delete(packages: pd.DataFrame):
    """
    Identify development versions of packages that are older than the second latest non-development version.

    This function takes a DataFrame containing package information and identifies development versions
    (e.g., `1.3.2.dev456`) that are older than the second latest non-development version (e.g., `1.3.2`).

    Note: `packaging.version` considers `1.3.2` to be more recent than `1.3.2.dev456`.

    Parameters:
    packages (pd.DataFrame): A DataFrame containing package information with columns 'name' and 'version'.
        Example:
            name    version
            ----    -------
            pkg1    1.3.3
            pkg1    1.3.3.dev123
            pkg1    1.3.2
            pkg1    1.3.2.dev456
            pkg1    1.3.1
            pkg1    1.3.1.dev678

    Returns:
    list: A list of dictionaries representing the packages to delete. Each dictionary contains package information.
        Example:
        [
            {'name': 'pkg1', 'version': '1.3.2.dev456'},
            {'name': 'pkg1', 'version': '1.3.1.dev678'}
        ]
    """
    grouped = packages.groupby('name')

    packages_to_delete = []

    for _, group in grouped:
        # Convert version strings to version objects for comparison
        group['parsed_version'] = group['version'].apply(version.parse)
        sorted_group = group.sort_values(by='parsed_version', ascending=False)
        # Find the non-dev versions
        latest_non_dev = sorted_group[
            ~sorted_group['parsed_version'].apply(lambda x: x.is_prerelease)
        ]
        if len(latest_non_dev) < 2:
            continue
        # Find the second latest non dev version. (eg: given 1.3.3, 1.3.2, 1.3.1, we want: 1.3.2)
        second_latest_non_dev = latest_non_dev.iloc[1]
        second_latest_non_dev_version = second_latest_non_dev['parsed_version']
        # Add dev versions older than the second latest non-dev version to the delete list
        # (eg: 1.3.3, 1.3.3.dev123, 1.3.2, 1.3.2.dev456, 1.3.1, 1.3.1.dev678 -> [1.3.2.dev456, 1.3.1.dev678])
        dev_versions_to_delete = sorted_group[
            (sorted_group['parsed_version'] < second_latest_non_dev_version)
            & (sorted_group['parsed_version'].apply(lambda x: x.is_prerelease))
        ]
        packages_to_delete.extend(dev_versions_to_delete.to_dict('records'))
    return packages_to_delete


def delete_old_packages(packages_to_delete: list[dict]):
    """
    Deletes old packages from a GitLab project registry.

    This function takes a list of dictionaries representing packages to delete. Each dictionary should contain
    the package's 'id', 'name', and 'version'. The function will attempt to delete each package by making a DELETE
    request to the GitLab API.

    Note: This function requires the `CI_PRIVATE_TOKEN` environment variable to be set. The token should have
    maintainer access to the project's registry.

    Parameters:
    packages_to_delete (list[dict]): A list of dictionaries, each containing:
        - 'id' (int): The package ID.
        - 'name' (str): The package name.
        - 'version' (str): The package version.

    Example:
    packages_to_delete = [
        {"id": 123, "name": "pkg1", "version": "1.3.2.dev456"},
        {"id": 124, "name": "pkg1", "version": "1.3.1.dev678"}
    ]

    Environment Variables:
    CI_PRIVATE_TOKEN (str): The private token for authentication with the GitLab API.

    Raises:
    ValueError: If the `CI_PRIVATE_TOKEN` environment variable is not set.
    requests.exceptions.RequestException: If there is an issue with the request to delete a package.
    """
    CI_ACCESS_TOKEN = get_env_var('CI_ACCESS_TOKEN')
    headerToken = {'PRIVATE-TOKEN': CI_ACCESS_TOKEN}

    for package_info in packages_to_delete:
        package_id = package_info['id']
        package_info_label = f"{package_info['name']} - v{package_info['version']}"
        print(f'Deleting package {package_info_label}')
        url = f'https://gitlab.mpcdf.mpg.de/api/v4/projects/{CI_PROJECT_ID}/packages/{package_id}'
        try:
            response = requests.delete(url, headers=headerToken)
            if response.status_code == 204:
                print(f'Package {package_info_label} has been deleted successfully')
            else:
                print(
                    f'/!\\ Unable to delete package {package_info_label}. Status code: {response.status_code}'
                )
        except requests.exceptions.RequestException as e:
            print(f'/!\\ Error: {e}')


if __name__ == '__main__':
    packages = fetch_packages()
    packages_to_delete = find_packages_to_delete(packages)
    df = pd.DataFrame(packages_to_delete)
    print('Packages to be deleted: \n', df[['parsed_version']])
    delete_old_packages(packages_to_delete)
