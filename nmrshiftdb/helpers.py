"""NMRshiftDB Import Helpers.

This module includes functions used to export NMRShiftDB database as NMReData files and to use those files to detect the locations of the raw NMR files and download them too. 
"""

import os
import wget
import shutil
import zipfile

def get_name(molecule):
    """returns the nmrshiftdb id for each molecule"""
    start = molecule.find('nmrshiftdb2 ') + 12
    end = molecule.find('\n', start)
    name = molecule[start:end]
    
    return name

def write_nmredata(molecule):
    """Writes a molecule in an NMReData file named after its nmrshiftdb id"""
    name = get_name(molecule)
    
    f = open(name+'.nmredata', "w")
    f.write(molecule)
    f.close()
    
    pass

def download_zips(molecule):
    """download all the raw NMR files from the links found in the NMReData file."""
    name = get_name(molecule)
    if not os.path.exists(name):
        os.makedirs(name)
        os.chdir('./'+ name)

        write_nmredata(molecule)

        index = 0
        while molecule.find('http', index) != -1:
            start = molecule.find('http', index)
            end = molecule.find('rawdata', start) + 7
            index = end

            location = molecule[start:end]
            wget.download(location)        

        os.chdir("..")
    pass

def create_spec_folders():
    for molecule in os.listdir("./"): 
        if os.path.isdir(molecule):
            molFolder = os.getcwd() + '/' + molecule
            for spectrum in os.listdir(molFolder):
                if 'zip' in spectrum:
                    name = spectrum[:spectrum.find('.')]
                    specFolder = molFolder + '/' + name
                    os.makedirs(specFolder)
                    shutil.move(molFolder+ '/' +spectrum, specFolder)
                    
    pass

def unzip_spec_files():
    print('unzipping files')
    print('unzipping the following files has failed. Please try to unzip them manually:')
    for molecule in os.listdir("./"): 
        if os.path.isdir(molecule):
            molFolder = os.getcwd() + '/' + molecule
            for spectrum in os.listdir(molFolder):
                specFolder = molFolder + '/' + spectrum
                if os.path.isdir(specFolder):
                    for file in os.listdir(specFolder):
                        if '.zip' in specFolder + '/' +file:
                            try:
                                with zipfile.ZipFile(specFolder + '/' +file, 'r') as zip_ref:
                                    zip_ref.extractall(specFolder)
                                os.remove(specFolder + '/' +file)
                            except:
                                pass
                            
    for path, directories, files in os.walk("."):
        for file in files:
            if 'zip' in file:
                try: 
                    with zipfile.ZipFile(os.path.join(path, file), 'r') as zip_ref:
                        zip_ref.extractall(path)
                    os.remove(os.path.join(path, file))
                except:
                    print('unzipping the following files has failed. Please try to unzip them manually:')
                    print(os.path.join(path, file))
                    
          
    pass

          
def structure_folders():
    print('preparing the folders structure for proper submision to nmrXiv. This might take a while')
    for molecule in os.listdir("./"): 
        if os.path.isdir(molecule):
            print(molecule)
            molFolder = os.getcwd() + '/' + molecule
            for spectrum in os.listdir(molFolder):
                specFolder = molFolder + '/' + spectrum
                if os.path.isdir(specFolder):
                    innerFolder = specFolder
                    while ("acqu" not in os.listdir(innerFolder)) and ("fid" not in os.listdir(innerFolder)):
                        for item in os.listdir(innerFolder):
                            if os.path.isdir(innerFolder + '/' + item):
                                innerFolder +=  '/' + item


                    for f in os.listdir(innerFolder):
                        try:
                            shutil.move(innerFolder + '/'+ f, specFolder) 
                        except:
                            print(innerFolder)
    pass