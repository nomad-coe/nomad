# -*- coding: utf-8 -*-
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

import pandas as pd
import json
import numpy
from collections import defaultdict
from datetime import date, datetime
import os


class PerovskiteEntryWriter():

    def __init__(self, csv_database_path: str):
        """
        Init method for the PerovskiteDBReader class.
        :param csv_database_path: path to the .csv file ocntaining the database
        :param data_schema_path: path to the excel file ocntaining the database quantities,
        db original data types, and quantities descriptions.

        """
        self.csv_database_path = csv_database_path
        self.df_db = pd.read_csv(self.csv_database_path, skiprows=0, parse_dates=["Ref_publication_date"])

    def read_columns(self):
        '''
        Reading the column names of the csv perovskite database file.
        '''

        self.column_names = self.df_db.columns
        return self.column_names

    def collect_schema_items(self):
        '''
        Generates a list of first level section names and quantities by splitting
        the column name in the first underscore found.
        A manual section name is also added (`Perovskite_deposition`)
        to be splitted in the second underscore found for one exception.
        '''

        self.column_names = self.read_columns()
        self.quantities_list = []
        self.sections = []

        for column in self.column_names:
            if 'Perovskite_deposition' in column:
                self.quantities_list.append('_'.join(column.split('_', 2)[2:]))
                self.sections.append('_'.join(column.split('_', 2)[:2]))
            elif '_' in column:
                self.quantities_list.append(column.split('_', 1)[1])
                self.sections.append(column.split('_', 1)[0])
            else:
                self.quantities_list.append(column)
                self.sections.append(column)

        self.sections = [i.lower() for i in self.sections]

        return self.sections, self.quantities_list

    def collect_schema_items_per_column(self, column):
        '''
        Generates a section and a quantity names by splitting the column name in the first
        underscore found. A manual section name is also added to be splitted in the
        second underscore found for one exception.
        '''

        self.column_names = self.read_columns()
        self.quantity = ''
        self.section = ''

        if 'Perovskite_deposition' in column:
            self.quantity = '_'.join(column.split('_', 2)[2:])
            self.section = '_'.join(column.split('_', 2)[:2])
        elif '_' in column:
            self.quantity = column.split('_', 1)[1]
            self.section = column.split('_', 1)[0]
        else:
            self.quantity = column
            self.section = column

        self.section.lower()

        return self.section, self.quantity

    def entry_writer(self, target_dir):
        '''
        Writes archive entries of the perovskite database
        :param target_dir: path pointing to where to dump the ``archive.json`` files
        '''

        self.sections, self.quantities_list = self.collect_schema_items()

        def create_rec_dd():  # To create an infinite level recursive default dictionary
            return defaultdict(create_rec_dd)

        for row in range(self.df_db.shape[0]):
            target_dict = create_rec_dd()
            for column in range(len(self.quantities_list)):

                # target_dict
                target_dict["data"]
                target_dict["data"]["m_def"] = "nomad.datamodel.metainfo.eln.perovskite_solar_cell_database.PerovskiteSolarCell"
                if not pd.isnull(self.df_db.iloc[row][column]):
                    target_dict["data"][self.sections[column]][self.quantities_list[column]] = self.df_db.iloc[row][column]
                target_dict = dict(target_dict)
                id = target_dict['data']["ref"]['ID']
                sufix = '.archive.json'
            save_path = os.path.join(target_dir, str(id) + sufix)
            with open(save_path, 'w') as fp:
                json.dump(target_dict, fp, cls=MyEncoder)


class MyEncoder(json.JSONEncoder):

    def default(self, obj):  # pylint: disable=E0202
        if isinstance(obj, numpy.integer):
            return int(obj)
        elif isinstance(obj, numpy.floating):
            return float(obj)
        elif isinstance(obj, numpy.bool_):
            return bool(obj)
        elif isinstance(obj, numpy.ndarray):
            return obj.tolist()
        elif isinstance(obj, (datetime, date,)):
            return obj.isoformat()
        else:
            return super(MyEncoder, self).default(obj)
