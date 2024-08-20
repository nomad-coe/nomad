# Updates the values.yaml file based on the nomad configuration. Currently this
# includes the following:
# - populating jupyterhub.singleUser.profileList with NORTH tools

from typing import Any, Dict
import os
from sys import stdout, argv
from ruamel.yaml import YAML
from nomad import config

dir_path = os.path.dirname(os.path.realpath(__file__))
file_path = os.path.join(dir_path, 'values.yaml')

yaml = YAML()
yaml.allow_duplicate_keys = True
with open(file_path, 'r') as file:
    data = yaml.load(file)

# Profiles currently break how NOMAD interacts with JupyterHub, so we're only using
# the extraImages to feed the prePuller
generate_profiles = False

if generate_profiles:
    profile_list = (
        data.setdefault('jupyterhub', {})
        .setdefault('singleuser', {})
        .setdefault('profileList', [])
    )
    for name, tool in config.north.tools.filtered_items():
        profile = next(
            (profile for profile in profile_list if profile['display_name'] == name),
            None,
        )
        if profile is None:
            profile = dict()
            profile_list.append(profile)

        profile.update(
            dict(
                display_name=name,
                description=tool.description,
                kubespawner_override=dict(
                    image=tool.image,
                    image_pull_policy=tool.image_pull_policy,
                ),
            )
        )

        if tool.default_url:
            profile['kubespawner_override']['default_url'] = tool.default_url
        if tool.cmd:
            profile['kubespawner_override']['cmd'] = tool.cmd
        if tool.privileged:
            profile['kubespawner_override']['privileged'] = tool.privileged
            profile['kubespawner_override']['allow_privilege_escalation'] = True
            profile['kubespawner_override']['uid'] = 0

else:
    pre_puller = data.setdefault('jupyterhub', {}).setdefault('prePuller', {})
    extra_images: Dict[str, Any] = {}
    pre_puller['extraImages'] = extra_images

    for name, tool in config.north.tools.filtered_items():
        try:
            image_name, image_tag = tool.image.rsplit(':', 1)
        except ValueError:
            image_name, image_tag = tool.image, 'latest'

        extra_images[name] = dict(name=image_name, tag=image_tag)


if len(argv) == 2:
    with open(argv[1], 'w') as file:
        yaml.dump(data, file)
else:
    yaml.dump(data, stdout)
