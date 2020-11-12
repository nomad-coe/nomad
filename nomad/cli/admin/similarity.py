#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import os
import json
from nomad import infrastructure
from nomad.datamodel.material import Material, Similarity, DOSSimilarity
from nomad.archive import ArchiveReader, ArchiveWriter
from typing import List


def ingest(input_path: str, batch_size: int, verbose: bool):
    """Used to ingest the given DOS similarity values into MongoDB.

    Args:
        input_path: Path of the msgpack file to ingest.
        batch_size: Batch size for MongoDB bulk ingest.
        verbose: Enable verbose output.
    """
    # Initialize mongo connection
    infrastructure.setup_mongo()

    bulk = Material.m_def.a_mongo.mongo_cls()._get_collection().initialize_ordered_bulk_op()  # pylint: disable=not-callable

    with ArchiveReader(input_path) as reader:
        i = 0
        for material_id in reader:
            material_object = reader[material_id]
            material_dict = material_object.to_dict()

            # Each entry is cycled through the material metainfo definition and
            # the mongo annotation to validate it and only push the annotated
            # data.
            mongo_instance = Material.m_def.a_mongo.mongo_cls(**material_dict, _created=False)  # pylint: disable=not-callable
            material_dict = mongo_instance.to_mongo().to_dict()

            # Add an upsert to the bulk operations
            bulk.find({'_id': material_dict["_id"]}).upsert().update_one({'$set': material_dict})
            i += 1
            if i % batch_size == 0:
                bulk.execute()
                bulk = Material.m_def.a_mongo.mongo_cls()._get_collection().initialize_ordered_bulk_op()  # pylint: disable=not-callable
                if verbose:
                    print("{} inserted".format(i))

    # Final push for the remainder
    bulk.execute()
    if verbose:
        print("{} inserted".format(i))


def update(input_dir: str, output_path: str, verbose: bool):
    """Used to create a compact msgpack file that follows the metainfo schema
    for materials and contains the DOS similarity values.

    Args:
        input_dir: Path of the directory containing the raw similarity filesa.
        output_path: Path of the output msgpack file.
        verbose: Enable verbose output.
    """
    # Find all valid data files in the given directory
    similarity_files: List[str] = []
    for filename in os.listdir(input_dir):
        if filename.endswith(".dat"):
            similarity_files.append(os.path.join(input_dir, filename))
    n_files: int = len(similarity_files)
    if n_files == 0:
        raise ValueError("Could not find similarity files in directory: {}".format(input_dir))
    if verbose:
        print("{} files found".format(n_files))

    # Gather the number of entries to prepare the msgpack file.
    n_entries = 0
    for filepath in similarity_files:
        num_lines = sum(1 for line in open(filepath, "r"))
        n_entries += num_lines
    if verbose:
        print("{} entries found".format(n_entries))

    # Read the data file containing similarities
    # i = 0
    if verbose:
        print("Writing msgpack file...")
    with ArchiveWriter(output_path, n_entries, entry_toc_depth=1) as writer:
        for filepath in similarity_files:
            with open(filepath, "r") as f:
                for line in f:
                    ientry = json.loads(line)

                    for key, value in ientry.items():
                        _, _, imaterial = key.split(":")

                        # Create data according to a metainfo model
                        material = Material()
                        material.material_id = imaterial
                        similarity = material.m_create(Similarity)
                        dos_similarity = similarity.m_create(DOSSimilarity)
                        ids = []
                        values = []

                        for jkey, similarity_value in value.items():
                            _, _, jmaterial = jkey.split(":")
                            ids.append(jmaterial)
                            values.append(similarity_value)

                        dos_similarity.material_ids = ids
                        dos_similarity.values = values

                        # Save as msgpack
                        writer.add(imaterial, material.m_to_dict())

    if verbose:
        print("Finished")
