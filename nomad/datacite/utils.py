import secrets

from nomad.config import config


def generate_target_url(doi: str, target_type: str) -> str:
    """Generates the URL the DOI resolves to."""

    return f'{config.gui_url()}/{target_type}/doi/{doi}'


def generate_unique_doi_name() -> str:
    """Generates a unique DOI name of pattern {prefix}/nomad.{random}.

    The prefix is defined in the config. The random part is a base32 Crockford
    string of length 8, split into two groups of 4 characters.
    """
    RANDOM_LENGTH = 8
    RANDOM_SPLIT = 4

    prefix = config.datacite.prefix
    namespace = 'nomad.'
    random_str = generate_random_b32crockford(RANDOM_LENGTH, RANDOM_SPLIT)

    return f'{prefix}/{namespace}{random_str}'


def generate_random_b32crockford(length: int = 8, split: int = 4) -> str:
    """Returns a random base32 Crockford string (lower-case, without check symbol).

    The alphabet includes digits 0-9 and lowercase letters a-z except i, l, o, and u.
    """
    alphabet = '0123456789abcdefghjkmnpqrstvwxyz'  # 32 characters
    number = secrets.randbits(length * 5)  # 5 bits per character

    result = ''
    while number > 0:
        result += alphabet[number & 0b11111]
        number >>= 5
    result = result.ljust(length, '0')  # Pad with zeros if necessary

    if split > 0:
        result = '-'.join(result[i : i + split] for i in range(0, len(result), split))

    return result
