import sys
from nomad.metainfo import Environment
from nomad.metainfo.legacy import LegacyMetainfoEnvironment
import nomad.datamodel.metainfo.common
import nomad.datamodel.metainfo.public
import nomad.datamodel.metainfo.general

m_env = LegacyMetainfoEnvironment()
m_env.m_add_sub_section(Environment.packages, sys.modules['nomad.datamodel.metainfo.common'].m_package)  # type: ignore
m_env.m_add_sub_section(Environment.packages, sys.modules['nomad.datamodel.metainfo.public'].m_package)  # type: ignore
m_env.m_add_sub_section(Environment.packages, sys.modules['nomad.datamodel.metainfo.general'].m_package)  # type: ignore

m_env.m_add_sub_section(Environment.packages, sys.modules['nomad.datamodel.datamodel'].m_package)  # type: ignore
m_env.m_add_sub_section(Environment.packages, sys.modules['nomad.datamodel.dft'].m_package)  # type: ignore
m_env.m_add_sub_section(Environment.packages, sys.modules['nomad.datamodel.ems'].m_package)  # type: ignore
m_env.m_add_sub_section(Environment.packages, sys.modules['nomad.datamodel.encyclopedia'].m_package)  # type: ignore
m_env.m_add_sub_section(Environment.packages, sys.modules['nomad.datamodel.optimade'].m_package)  # type: ignore
