/* eslint-disable quotes */
/*
 * Copyright The NOMAD Authors.
 *
 * This file is part of NOMAD. See https://nomad-lab.eu for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
import '@h5web/app/dist/style-lib.css';
import '@h5web/app/dist/style.css';

import React from 'react'
import H5Web from '../../visualization/H5Web';
import { Paper } from '@material-ui/core'

const NexusCard = (props) => {
    const filepath = props.index.upload_id.substring(0,2)+"/"+props.index.upload_id+"/raw/"+props.index.mainfile

    if (!props.index.parser_name === "parsers/nexus")
      return null
    else
      return (
          <Paper elevation={2} style={{height: '50vh'}}>
            <H5Web filepath={filepath}/>
        </Paper>
      );

}

export default React.memo(NexusCard)