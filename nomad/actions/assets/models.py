from enum import Enum
from typing import Final, Literal

from pydantic import BaseModel, ConfigDict, Field


class ActionAssetPurpose(str, Enum):
    ACTION_START = 'action_start'
    ACTION_SIGNAL = 'action_signal'


ACTION_ASSET_REF_TYPE: Final = 'action_asset_ref'


class ActionAssetRef(BaseModel):
    model_config = ConfigDict(populate_by_name=True)

    nomad_type: Literal['action_asset_ref'] = Field(alias='_nomad_type')
    filename: str
    media_type: str | None = None
    size: int | None = None
    sha256: str | None = None


class ActionAssetUploadResult(ActionAssetRef):
    purpose: ActionAssetPurpose
    action_id: str | None = None
    action_instance_id: str | None = None
    signal_fn_name: str | None = None
