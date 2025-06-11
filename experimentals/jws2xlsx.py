import sys
import olefile as ofio
from struct import unpack
import numpy as np
import os
import pandas as pd
from natsort import natsorted # pip install natsort

#from jwsprocessor
DATAINFO_FMT = '<LLLLLLdddLLLLdddd'
#from jwsprocessor

class DataInfo:
    def __init__(self, data):

        if len(data) < 96:
            raise Exception ("DataInfo should be at least 96 bytes!")

        data = data[:96]
        data_tuple = unpack(DATAINFO_FMT, data)
        self.npoints = data_tuple[5]
        self.x_for_first_point = data_tuple[6]
        self.x_for_last_point = data_tuple[7]
        self.x_increment = data_tuple[8]

def main(fn, channels_num):

    doc = ofio.OleFileIO(fn)

    if doc.exists('Y-Data'):
        DI = DataInfo(doc.openstream('DataInfo').read())
        str = doc.openstream('Y-Data')
        data = str.read()
        fmt = 'f'*DI.npoints
        chunk_size = int(len(data)/channels_num)
        data_chunked = [data[i:i + chunk_size] for i in range(0, len(data), chunk_size)]
        try:
            unpacked_data = [unpack(fmt, data_chunk) for data_chunk in data_chunked]
            #generate xdata
            x = DI.x_for_first_point
            x += DI.x_increment
            xdata = np.arange(DI.x_for_first_point,  # start
                                DI.x_for_last_point + DI.x_increment,  # end+incr.
                                DI.x_increment)  #

            unpacked_data.insert(0, xdata)

            #np.savetxt(jws_file+'.txt', np.c_[unpacked_data].transpose())
            converted_data = np.c_[unpacked_data].transpose()
        except:
            print("Incorrect number of channels.")
    else:
        print("error in file")
    return converted_data

# if __name__ == '__main__':
#     if len(sys.argv) > 1:
#         try:
#             channels_num = int(input("Enter the number of channels:"))
#         except:
#             print('Enter the proper number of channels.')
#             sys.exit(1)


#         infn = sys.argv[1]
#         outfn = sys.argv[2] if len(sys.argv) > 2 else None
#         os.chdir(infn)
#         matches = []
#         for root, dirnames, filenames in os.walk(infn):
#             for filename in filenames:
#                 if filename.endswith(('.jwb','.jws')):
#                     matches.append(os.path.join(root, filename))
#         all_data = []
#         for jws_file in matches:
#             print(jws_file)
#             cur_data = main(jws_file, channels_num)
#             all_data.append(cur_data)
#         np.savetxt('all_data.txt', all_data)

#     else:
#         print ('USAGE:python jws2txt.py input_folder')

channels_num = 3
infn = os. getcwd()
matches = []
filename_list = ["wavelength"]
for root, dirnames, filenames in os.walk(infn):
    for filename in filenames:
        if filename.endswith(('.jwb','.jws')):
            filename_list.append(filename)
            matches.append(os.path.join(root, filename))
            

matches = natsorted(matches) ### sort the value according to number not by char
filename_list[1:] = natsorted(filename_list[1:])


all_mdeg = []
all_ab = []
all_HT = []
for jws_file in matches:
    cur_data = main(jws_file, channels_num)
    all_mdeg.append(cur_data[:,1])
    all_ab.append(cur_data[:,2])
    all_HT.append(cur_data[:,3])
wavelength = [cur_data[:,0]]   
    
all_mdeg = np.array(all_mdeg)
all_mdeg = np.concatenate((np.array(wavelength),all_mdeg), axis = 0).transpose()

all_ab = np.array(all_ab)
all_ab = np.concatenate((np.array(wavelength),all_ab), axis = 0).transpose()

all_HT = np.array(all_HT)
all_HT = np.concatenate((np.array(wavelength),all_HT), axis = 0).transpose()


df1 = pd.DataFrame(all_mdeg)
df1.columns = filename_list

df2 = pd.DataFrame(all_ab)
df2.columns = filename_list

df3 = pd.DataFrame(all_HT)
df3.columns = filename_list
#np.savetxt('all_data.txt', all_data)
# df.to_excel("all_data.xlsx",
#              sheet_name='mdeg')
with pd.ExcelWriter('all_data.xlsx') as writer:  
    df1.to_excel(writer, sheet_name='all_mdeg',index=False)
    df2.to_excel(writer, sheet_name='all_absorbance',index=False)
    df2.to_excel(writer, sheet_name='all_HT',index=False)
