# Import from Log-Files Created by the PVD software by Roland Mainz
# for PVD Measurementfile v3, 3.0.0
# Version 0.1
# by Pascal Becker
import re
from typing import Dict, List, Any, Tuple
import json
from datetime import datetime
import numpy as np
import pandas as pd
from xarray import Dataset
import os


class Importer:
    def __init__(self):
        self.import_name = "pvd_p_evaporation"
        self.import_type = "processes"
        self.import_long_name = "PVD-P software log"
        self.supported_file_types = "Process logs (*.csv)"

    @staticmethod
    def get_unit(column_name: str) -> str:
        # Eurotherm values
        if column_name.endswith(' Aout'):
            return '%'
        if column_name.endswith(' Pv'):
            return '0.1 nmol cm$^{-1}$ s$^{-1}$'
        if column_name.endswith(' T') or column_name.endswith(' Tsp'):
            return '$^{\\circ}$C'

        # Inficon QCM values
        if column_name[:-1].endswith('Ch_Time_'):
            return 's'
        if column_name[:-1].endswith('C_Thik_'):
            return 'kA$^{\\circ}$'
        if column_name[:-1].endswith('Ratet1_'):
            return 'A$^{\\circ}$ s$^{-1}$'
        if column_name[:-1].endswith('Xlife_'):
            return '%'

        # Pressure
        if column_name[:-1].endswith(' Pressure'):
            return 'mbar'

        # none of the above (default)
        return ''

    def read(self, f, comment: str = None, user: str = None) -> Tuple[List[Dataset], List[None]]:
        header = '{"'
        while True:
            line = f.readline()
            if line[0] == '#':
                if ':' in line:
                    header += line[2:-1].replace(': ', '": "') + '", "'
            else:
                process_data = pd.read_csv(f, sep='\t')
                break

        header = header[:-3] + '}'
        header_dict = json.loads(header)
        header_dict['datetime'] = datetime.strptime(header_dict['Date'] + ' ' + header_dict['Time'],
                                                    '%Y/%m/%d %H:%M:%S')

        # Check if list of substrate_labels in log file matches with selected libraries
        substrate_labels = re.split('[;,]', header_dict['Substrate Number'].replace(' ', ''))

        evaporation_start_index = np.where(process_data['QCM1 SHTSRC_1']
                                            + process_data['QCM1 SHTSRC_2']
                                            + process_data['QCM2 SHTSRC_1']
                                            + process_data['QCM2 SHTSRC_2'] > 0)[0][0]

        datetime_column = []
        process_time_column = []
        for process_time in process_data['Time']:
            datetime_column.append(datetime.strptime(header_dict['Date'] + ' ' + process_time, '%Y/%m/%d %H:%M:%S'))
        evaporation_start_time = datetime_column[evaporation_start_index]
        for process_time in datetime_column:
            if process_time < evaporation_start_time:
                process_time_column.append(-1 * (evaporation_start_time - process_time).seconds)
            else:
                process_time_column.append((process_time - evaporation_start_time).seconds)

        return process_time_column, process_data['Substrate PV'].to_numpy()
        # All dimensions and data prepared for the dataset.
        # Start by making the dimension arrays,

        # x_dim_array = make_dimension_array(         ### will only get you in trouble, uses Kameleont functions
        #     dim_name="x",
        #     dim_values=np.array([0]),
        #     dim_units="mm",
        #     dim_long_name="x position")

        # y_dim_array = make_dimension_array(         ### will only get you in trouble, uses Kameleont functions
        #     dim_name="y",
        #     dim_values=np.array([0]),
        #     dim_units="mm",
        #     dim_long_name="y position"
        # )

        #### uses Kameleont function. Needs to be replaced by something Nexus-like !!! You may want to use the parameters of the read function 'comment' and 'user' for this?!
        # time_dim_array = make_dimension_array(
        #     dim_name="time",
        #     dim_values=np.array(process_time_column),
        #     dim_units="s",
        #     dim_long_name="process time"
        # )

        # pvd_data_das = []
        for column_name in process_data.columns[:-1]:  # last column is empty because column row ends with tab
            if column_name.replace('.', '').replace(' ', '_').lower() != 'time':
                column_unit = self.get_unit(column_name)
                column_data = process_data[column_name].values[:].reshape(1, 1, -1)  ### this reshape is not necessary, if you make your data only depend on the time

                ### uses Kameleont function. Needs to be replaced by something Nexus-like !!!
                ### Generates technically multidemensional arrays depending on x and y.
                ### in reality, data only depends on time, so  x = process_time_column, y = column_data
                # pvd_data_das.append(make_data_array(
                #     data_name=column_name.replace('.', '').replace(' ', '_').lower(),
                #     data=column_data,
                #     x_dim=x_dim_array,
                #     y_dim=y_dim_array,
                #     additional_dims=[time_dim_array],
                #     data_units=column_unit,  # TODO: fill units dependent on file type
                #     data_long_name=column_name
                # ))

        ### This block compiles our data-Arrays into data sets. Don't know, how this would be handled in Nexus files
        # return_ds = []
        # return_json = []
        # # Finally, make the dataset from a list of the data arrays.
        # for library_info in library_infos:
        #     comment = None
        #     if required_params["comment"] != '':
        #         comment = required_params["comment"]
        #     return_ds.append(make_dataset(
        #         data_arrays=pvd_data_das,
        #         data_name=self.import_name,
        #         data_type=self.import_type,
        #         library_id=library_info["identifiers"]["library_id"],
        #         operator=header_dict['operator'],
        #         data_datetime=header_dict['datetime'],
        #         importing_user=user,
        #         comment=comment,
        #         additional_attrs={'process_id': header_dict.get('process ID', '')}
        #     ))

            ### This block generates information for the library_info.json, i.e. adds inforamtion
            ### that a new layer was deposited, which has x-y dimensions of abc and is a solar cell absorber
            # # get substrate dimensions
            # x_length = None
            # y_length = None
            # for layer in library_info["layers"]:
            #     if layer['type'].lower() == 'substrate':
            #         physical_prop = layer["properties"]["physical"]
            #         x_length = (physical_prop["x_length"][0]["value_low"],
            #                     physical_prop["x_length"][0]["value_high"])
            #         y_length = (physical_prop["y_length"][0]["value_low"],
            #                     physical_prop["y_length"][0]["value_high"])
            #         break
            #
            # if x_length and y_length:
            #     set_layer_dim = False
            #     # tyring to assign a substrate holder for the sample
            #     if 23 < x_length[0] < 27 and 23 < y_length[0] < 27:  # inch x inch substrate frame
            #         x_length = (23 - 0.5, 23)  # opening in the substrate holder
            #         y_length = (23 - 0.5, 23)  # estimating 0.5 mm possible shading
            #         set_layer_dim = True
            #     elif 48 < x_length[0] < 51 and 48 < y_length[0] < 51:  # 50 x 50 substrate frame
            #         x_length = (48 - 0.5, 48)  # opening in the substrate holder
            #         y_length = (48 - 0.5, 48)
            #         set_layer_dim = True
            #     # for stripe substrates, the dimensions of the layer depend on it's position in the substrate frame.
            #
            #     temp_knc_path = generate_knc_path(return_ds[-1])
            #     temp_layer = LibraryLayer(library_id=library_info["identifiers"]["library_id"],
            #                               layer_type='absorber',
            #                               name='Perovskite Absorber',
            #                               creation_datetime=header_dict['datetime'],
            #                               position=-1)
            #     if set_layer_dim:
            #         temp_layer.add_property(keys=["physical", "x_length"],
            #                                 value=x_length,
            #                                 unit="mm",
            #                                 source=temp_knc_path
            #                                 )
            #         temp_layer.add_property(keys=["physical", "y_length"],
            #                                 value=y_length,
            #                                 unit="mm",
            #                                 source=temp_knc_path
            #                                 )
            #     return_json.append(temp_layer)

        return None ### actually your things, for the nexus file.   # return_ds, return_json
