from nomad import infrastructure
from nomad.groups import create_user_group


def _create(group_id, group_name, owner, members=None):
    members = members or []
    return create_user_group(
        group_id=group_id, group_name=group_name, owner=owner, members=members
    )


def init_gui_test_groups():
    user0 = infrastructure.user_management.get_user(username='admin').user_id
    user1 = infrastructure.user_management.get_user(username='test').user_id
    user2 = infrastructure.user_management.get_user(username='scooper').user_id
    user3 = infrastructure.user_management.get_user(username='ttester').user_id

    groups = {
        'group0': ('group0', 'Group Admin', user0),
        'group1': ('group1', 'Group Test', user1),
        'group2': ('group2', 'Group Cooper', user2),
        'group3': ('group3', 'Group Tester', user3),
        'group23': ('group23', 'Group 23', user2, [user3]),
    }
    groups = {k: _create(*args) for k, args in groups.items()}

    return groups
