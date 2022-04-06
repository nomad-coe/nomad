# Import from Log-Files Created by the PVD software by Roland Mainz
# for PVD Measurementfile v3, 3.0.0
# Version 0.1
# by Pascal Becker

import json
from datetime import datetime
import numpy as np
import pandas as pd


class PVDImporter:
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

    def read(self, f, comment: str = None, user: str = None):
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
        header_dict['datetime'] = datetime.strptime(
            header_dict['Date'] + ' ' + header_dict['Time'], '%Y/%m/%d %H:%M:%S')

        evaporation_start_index = np.where(
            process_data['QCM1 SHTSRC_1'] + process_data['QCM1 SHTSRC_2'] + process_data['QCM2 SHTSRC_1'] + process_data['QCM2 SHTSRC_2'] > 0)[0][0]

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

        return process_time_column, process_data['Substrate PV'].to_numpy(),\
            process_data['Substrate TSP'].to_numpy(), process_data['Vacuum Pressure1'].to_numpy()
