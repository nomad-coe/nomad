# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
This module contains the objects and setups for running pipelines.
"""
from typing import List
from collections import defaultdict
import networkx as nx
from celery import chain, group, chord
from celery.exceptions import SoftTimeLimitExceeded

from nomad.processing.base import NomadCeleryTask, PROCESS_CALLED, PROCESS_COMPLETED
from nomad.processing.celeryapp import app
import nomad.processing.data
from nomad import config

import json
from kombu.serialization import register


def my_dumps(obj):
    """Custom JSON encoder function for Celery tasks.
    """
    class MyEncoder(json.JSONEncoder):
        def default(self, obj):  # pylint: disable=E0202
            if isinstance(obj, PipelineContext):
                return obj.encode()
            else:
                return json.JSONEncoder.default(self, obj)

    return json.dumps(obj, cls=MyEncoder)


def my_loads(obj):
    """Custom JSON decoder function for Celery tasks.
    """
    def my_decoder(obj):
        if '__type__' in obj:
            if obj['__type__'] == '__pipelinecontext__':
                return PipelineContext.decode(obj)
        return obj
    return json.loads(obj, object_hook=my_decoder)


register("myjson", my_dumps, my_loads, content_type='application/x-myjson', content_encoding='utf-8')


class PipelineContext():
    """Convenience class for storing pipeline execution related information.
    Provides custom encode/decode functions for JSON serialization with Celery.
    """
    def __init__(self, filepath, parser_name, calc_id, upload_id, worker_hostname, re_process=False):
        self.filepath = filepath
        self.parser_name = parser_name
        self.calc_id = calc_id
        self.upload_id = upload_id
        self.worker_hostname = worker_hostname
        self.re_process = re_process

    def encode(self):
        return {
            "__type__": "__pipelinecontext__",
            "filepath": self.filepath,
            "parser_name": self.parser_name,
            "calc_id": self.calc_id,
            "upload_id": self.upload_id,
            "worker_hostname": self.worker_hostname,
            "re_process": self.re_process,
        }

    @staticmethod
    def decode(data):
        return PipelineContext(
            data["filepath"],
            data["parser_name"],
            data["calc_id"],
            data["upload_id"],
            data["worker_hostname"],
            data["re_process"],
        )


class Stage():
    """Stage comprises of a single python function. After this function is
    completed, the stage is completed.
    """
    def __init__(self, name: str, function):
        """
        Args:
            name: Name of the stage. The name is used in resolving stage
                dependencies.
            function: A regular python function that will be executed during
                this stage. The function should not return any values as all
                communication happens through object persistence in MongoDB. The
                function should accept the following arguments:

                    - context: PipelineContext object
                    - stage_name: Name of the stage executing the function
                    - i_stage: The index of this stage in the pipeline
                    - n_stages: Number of stages in this pipeline
        """
        self.name = name
        self._function = function
        self._pipeline: Pipeline = None
        self.index: int = None
        self.dependencies: List[Stage] = []

    def add_dependency(self, name):
        self.dependencies.append(name)

    # def run(self, filepath, parser_name, calc_id, upload_id, worker_hostname, stage_name, i_stage, n_stages):
        # wrapper.delay(self._function.__name__, filepath, parser_name, calc_id, upload_id, worker_hostname, stage_name, i_stage, n_stages)

    def signature(self):
        return wrapper.si(
            self._function.__name__,
            self._pipeline.context,
            self.name,
            self.index,
            len(self._pipeline.stages)
        )


class Pipeline():
    """Pipeline consists of a list of stages. The pipeline is complete when all
    stages are finished.
    """
    def __init__(self, context: PipelineContext):
        """
        Args:
            context: The working context for this pipeline.
        """
        self.context = context
        self.stages: List[Stage] = []

    def add_stage(self, stage: Stage):
        """Adds a stage to this pipeline. The stages are executec in the order
        they are added with this function.

        Args:
            stage: The stage to be added to this pipeline.
        """
        stage._pipeline = self
        stage.index = len(self.stages)
        self.stages.append(stage)


@app.task(
    bind=True, base=NomadCeleryTask, ignore_results=False, max_retries=3,
    acks_late=config.celery.acks_late, soft_time_limit=config.celery.timeout,
    time_limit=config.celery.timeout)
def wrapper(task, function_name, context, stage_name, i_stage, n_stages):
    """This function wraps the function calls made within the stages. This
    wrapper takes care of common tasks like exception handling and timeout
    handling, and also persists the state of the pipelines in to MongoDB.
    """
    # Create the associated Calc object during first stage, otherwise get it
    # from Mongo
    if i_stage == 0:
        calc = nomad.processing.data.Calc.create(
            calc_id=context.calc_id,
            mainfile=context.filepath,
            parser=context.parser_name,
            worker_hostname=context.worker_hostname,
            upload_id=context.upload_id
        )
        # During the first step tell the calculation that it's processing has been
        # started
        if context.re_process:
            calc.current_process = "process"
        else:
            calc.current_process = "re_process"
        calc.process_status = PROCESS_CALLED
        calc.save()
    else:
        calc = nomad.processing.data.Calc.get(context.calc_id)
    logger = calc.get_logger()

    # Get the defined function. If does not exist, log error and fail calculation.
    function = globals().get(function_name, None)
    if function is None:
        calc.fail('Could not find the function associated with the stage.')

    # Try to execute the stage.
    try:
        function(context, stage_name, i_stage, n_stages)
    except SoftTimeLimitExceeded as e:
        logger.error('exceeded the celery task soft time limit')
        calc.fail(e)
    except Exception as e:
        calc.fail(e)
    except SystemExit as e:
        calc.fail(e)

    # After the last stage, save the calculation info
    if i_stage == n_stages - 1:
        calc.process_status = PROCESS_COMPLETED
        calc.save()


@app.task(
    bind=True, base=NomadCeleryTask, ignore_results=False, max_retries=3,
    acks_late=config.celery.acks_late, soft_time_limit=config.celery.timeout,
    time_limit=config.celery.timeout)
def empty_task(task, *args, **kwargs):
    """Empty dummy task used as a callback for Celery chords.
    """
    pass


@app.task(
    bind=True, base=NomadCeleryTask, ignore_results=False, max_retries=3,
    acks_late=config.celery.acks_late, soft_time_limit=config.celery.timeout,
    time_limit=config.celery.timeout)
def upload_cleanup(task, upload_id):
    """Task for cleanin up after processing.
    """
    upload = nomad.processing.data.Upload.get(upload_id)
    upload.processing_finished()
    upload.get_logger().info("Cleanup for upload {} finished.".format(upload_id))


def comp_process(context, stage_name, i_stage, n_stages):
    """Function for processing computational entries: runs parsing and
    normalization.
    """
    # Process calculation
    calc = nomad.processing.data.Calc.get(context.calc_id)
    if context.re_process is True:
        calc.re_process_calc()
        calc.get_logger().warn("Re-processing of calculation {} at path {} finished.".format(context.filepath, context.calc_id))
    else:
        calc.process_calc()
        calc.get_logger().warn("Processing of calculation {} at path {} finished.".format(context.filepath, context.calc_id))


def get_pipeline(context):
    """Used to fetch a pipeline based on a pipeline context. Typically chosen
    simply based on a matched parser name that is stored in the context

    Args:
        context: The context based on which the pipeline is chosen.

    Returns:
        Pipeline: The pipeline to execute for the given context.
    """
    pipeline = Pipeline(context)

    # Phonopy pipeline
    if context.parser_name == "parsers/phonopy":
        stage1 = Stage("comp_process_phonopy", comp_process)
        stage1.add_dependency("comp_process")
        pipeline.add_stage(stage1)
    # DFT pipeline
    else:
        stage1 = Stage("comp_process", comp_process)
        pipeline.add_stage(stage1)

    return pipeline


def run_pipelines(context_generator, upload_id) -> int:
    """Used to start running pipelines based on the PipelineContext objects
    generated by the given generator. Currently the implementation searches all
    the mainfiles in an upload before starting to execute any stages. This
    makes synchronization with Celery much easier.

    Returns:
        The number of pipelines that were started.
    """
    # Resolve all pipelines into disconnected dependency trees and run
    # each tree in parallel.
    stage_dependencies = []
    stages: defaultdict = defaultdict(list)
    stage_names = set()
    n_pipelines = 0
    for context in context_generator:
        pipeline = get_pipeline(context)
        n_pipelines += 1
        for stage in pipeline.stages:

            # Store stage names to be used as nodes
            stage_names.add(stage.name)

            # Store dependencies to be used as edges
            for dependency in stage.dependencies:
                stage_dependencies.append((stage.name, dependency))
            stages[stage.name].append(stage)

    if n_pipelines != 0:
        # Resolve all independent dependency trees
        dependency_graph = nx.DiGraph()
        dependency_graph.add_nodes_from(stage_names)
        dependency_graph.add_edges_from(stage_dependencies)
        dependency_trees = nx.weakly_connected_components(dependency_graph)

        # Form chains for each independent tree.
        chains = []
        for tree_nodes in dependency_trees:
            tree = dependency_graph.subgraph(tree_nodes).copy()
            sorted_nodes = nx.topological_sort(tree)
            groups = []
            for node in reversed(list(sorted_nodes)):
                # Group all tasks for a stage
                tasks = stages[node]
                task_signatures = []
                for task in tasks:
                    task_signatures.append(task.signature())
                task_group = group(*task_signatures)

                # Celery does not allow directly chaining groups. To simulate
                # this behaviour, we instead wrap the groups inside chords
                # which can be chained. The callback is a dummy function that
                # does nothing.
                groups.append(chord(task_group, body=empty_task.si()))

            # Form a chain of stages for this tree
            stage_chain = chain(*groups)
            chains.append(stage_chain)

        # This is is the group that will start executing all independent stage trees parallelly.
        tree_group = group(*chains)

        # After all trees are finished, the upload should perform cleanup
        final_chain = chain(tree_group, upload_cleanup.si(upload_id))
        final_chain.delay()

    return n_pipelines
