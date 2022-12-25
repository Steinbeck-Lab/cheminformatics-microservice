import os
import sys
import requests
from helpers import*


def main():
    print('NMRShiftDB downloading script has started. Please find downloaded items in the folder "output."\n')
    
    if not os.path.exists('output'):
        os.makedirs('output')
    os.chdir('./output')

    URL = "https://sourceforge.net/projects/nmrshiftdb2/files/data/nmrshiftdb2.nmredata.sd/download"
    response = requests.get(URL)
    open("nmrshiftdb2.nmredata.sd", "wb").write(response.content)

    f = open("nmrshiftdb2.nmredata.sd", "r")
    text = f.read()
    f.close() 

    molecules = text.split('$$$$\n')[:-1]
    print('The total number of molecules found in NMRShiftDB, including duplicates, is: '+ str(len(molecules)) + '\n')

    with_raw = []
    without_raw = []
    for molecule in molecules:
        if "http" in molecule:
            with_raw.append(molecule)
        else:
            without_raw.append(molecule)

    print('The number of molecules found in NMRShiftDB with raw NMR files, with duplicates, is: '+ str(len(with_raw)) + '\n')
    print('The number of molecules found in NMRShiftDB without raw NMR files, with duplicates, is: '+ str(len(without_raw)) + '\n')
    print('Writing NMReData files, this might take a while.')



    if not os.path.exists('without_raw'):
        os.makedirs('without_raw')
    os.chdir('./without_raw')

    for molecule in without_raw:
        write_nmredata(molecule)


    os.chdir("..")
    if not os.path.exists('with_raw'):
        os.makedirs('with_raw')
    os.chdir('./with_raw')

    print('Downloading experimental NMR files, this might take a while. Here you can see the molecules IDs from NMRShiftDB: ')
    for molecule in with_raw:
        download_zips(molecule)
    
    print('NMRShiftDB downloading is finished')
    
    create_spec_folders()
    unzip_spec_files()
    structure_folders()
        

if __name__ == "__main__":
    sys.exit(main())
