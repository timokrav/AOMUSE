import logging
from pathlib import Path
import time
from datetime import datetime
import os
from astropy.io import fits
from astropy.table import Table
from tqdm import tqdm
import numpy as np
import json
from pony.orm import *

provider = "mysql"
host = "127.0.0.1"
user = "username"
passwd = "password"
database_name = "DB_name"
rootDir = "/home/dir/to/data"

# FITS filenames to be ignored
ignore_list = ['PIXTABLE', 'MASTER', 'LINES', 'SKY', "STD_", "TRACE"]

##############################
# Actual script operations start from here, there should be no reason to edit these.
# HERE BE DRAGONS
##############################

db = Database()


# Target class describes an entry to the Target table in the database
# The classes have to inherit db.Entity from Pony
class Target(db.Entity):
    #   ----- Attributes -----

    target_name = Required(str, unique=True)  # Required: Cannot be None

    #   ----- Relations -----

    exposures = Set('Exposure')  # One target contains a set of exposures


# Exposure table class
class Exposure(db.Entity):
    #   ----- Attributes -----

    observation_time = Required(str, unique=True)
    insMode = Required(str)
    datacube_header = Optional(LongStr)
    raw_exposure_header = Optional(LongStr)
    raman_image_header = Optional(LongStr)
    pampelmuse_catalog = Optional(LongStr)
    psf_params = Optional(LongStr)
    sources = Optional(LongStr)

    #   ----- Sky parameters -----
    sky_condition_start_time = Optional(float)
    sky_condition_start = Optional(LongStr)
    sky_comment_start = Optional(LongStr)
    sky_condition_end_time = Optional(float)
    sky_condition_end = Optional(LongStr)
    sky_comment_end = Optional(LongStr)

    #   ----- Relations -----

    target = Required('Target')  # One exposure belongs to a target


# This simplifies some header value fetches
def fetch_data(header, keyword):
    try:
        return header[keyword]
    except KeyError:
        print('Keyword ' + keyword + ' not found.')


# This method converts an np.int64 to python native int, because the json library
# at the moment cannot deal with numpy int64
def convert_npint64_to_int(o):
    if isinstance(o, np.int64):
        return int(o)
    raise TypeError


@db_session  # Decorator from Pony to make the function into a db session, thus the database can be modified
def muse_script():
    """
    This function read 5 files (prm, psf, reduced, raman and raw) to store the information of each exposure in a database.

    For the prm file reads the PSFPARS extension and stores it in the 'psfParams' field of the exposure.
    For the psf file reads the 'PM_X' extensions, where X is the id of the source, and stores it in the
    'sources' field of the exposure.
    With other files the scripts takes the headers and saves them as json-like LongStr
    Also, from the primary header it stores the target and instrument mode of the exposure in their own fields    to make it easier the analysis from differents targets and instrument modes.
    """

    start = time.time()
    all_raw_files = []
    all_nightlog_files = []
    all_datacube_files = []
    all_raman_files = []
    all_prm_files = []
    all_psf_files = []

    skip_number = 0

    # FITS file scan and identification
    print("Generating list of fits files...")
    fits_list = list(Path(rootDir).rglob('*.fits'))

    print("Starting fits...")
    logging.info("Starting fits file identification.")
    time.sleep(0.1)  # Makes the terminal output pretty
    for filepath in tqdm(fits_list):
        if any(ignored_word in filepath.parts[-1] for ignored_word in ignore_list):
            skip_number += 1
            continue
        try:
            header = fits.getheader(filepath)
        except:
            continue  # Probably a directory, move on
        try:
            if 'HIERARCH ESO DPR TYPE' in header and header['HIERARCH ESO DPR TYPE'] == "SKY":
                continue  # Sky observation, ignore
            if 'HIERARCH ESO DPR CATG' in header and header['HIERARCH ESO DPR CATG'] == "SCIENCE":
                all_raw_files.append(filepath)
                continue
            if 'HIERARCH ESO PRO CATG' in header and header['HIERARCH ESO PRO CATG'] == "DATACUBE_FINAL":
                if 'PROV2' in header:
                    continue  # This file is a combined cube, we ignore them
                all_datacube_files.append(filepath)
                continue
            if 'HIERARCH ESO PRO CATG' in header and header['HIERARCH ESO PRO CATG'] == "RAMAN_IMAGES":
                all_raman_files.append(filepath)
                continue
        except KeyError:
            pass  # File in question is not a proper ESO FITS product
        if ".prm" in filepath.parts[-1]:
            all_prm_files.append(filepath)
        if ".psf" in filepath.parts[-1]:
            all_psf_files.append(filepath)
        #TODO Add a check for the upcoming psf-per-wavelength files
    logging.info("Finished with fits files.")

    # FITS FZ file scan and identification
    # Sometimes raw products are compressed so we doublecheck them for raw files
    print("Generating list of fz files...")
    fits_list = list(Path(rootDir).rglob('*.fits.fz'))

    print("Starting fz...")
    logging.info("Starting fz file scanning.")
    time.sleep(0.1)  # Makes the terminal output pretty
    for filepath in tqdm(fits_list):
        header = fits.getheader(filepath)
        if 'HIERARCH ESO DPR TYPE' in header and header['HIERARCH ESO DPR TYPE'] == "SKY":
            continue  # Sky observation, ignore
        if 'HIERARCH ESO DPR CATG' in header and header['HIERARCH ESO DPR CATG'] == 'SCIENCE':
            all_raw_files.append(filepath)
            nightlog_path = Path(str(filepath)[:-7] + "NL.txt")
            if Path.exists(nightlog_path):
                all_nightlog_files.append(nightlog_path)
    logging.info("Finished with fz files.")

    # detections.cat file scan and identification
    # pampelmuse source catalogs
    print("Generating list of pampelmuse cat files...")
    logging.info("Generating a list of pampelmuse source catalogs.")
    all_catalog_files = list(Path(rootDir).rglob("*.detections.cat"))


    end = time.time()
    print('{:.3f}'.format(end-start) + " seconds to finish file search")
    print('{:d}'.format(len(all_raw_files)) + " raw science exposures found")
    logging.info('{:d}'.format(len(all_raw_files)) + " raw science exposures found")
    logging.info('{:d}'.format(len(all_nightlog_files)) + " nightlogs found")
    print('{:d}'.format(len(all_nightlog_files)) + " nightlogs found")
    logging.info('{:d}'.format(len(all_datacube_files)) + " reduced datacubes found")
    print('{:d}'.format(len(all_datacube_files)) + " reduced datacubes found")
    logging.info('{:d}'.format(len(all_raman_files)) + " raman files found")
    print('{:d}'.format(len(all_raman_files)) + " raman files found")
    logging.info('{:d}'.format(len(all_prm_files)) + " PRM files found")
    print('{:d}'.format(len(all_prm_files)) + " PRM files found")
    logging.info('{:d}'.format(len(all_psf_files)) + " PSF files found")
    print('{:d}'.format(len(all_psf_files)) + " PSF files found")
    logging.info('{:d}'.format(len(all_catalog_files)) + " catalog files found")
    print('{:d}'.format(len(all_catalog_files)) + " catalog files found")
    logging.info('{:d}'.format(skip_number) + " files skipped")
    print('{:d}'.format(skip_number) + " files skipped")

    # Datacube extraction
    unique_observations = []
    reduced_cube_entries = []

    print("Parsing through found datacubes...")
    for reduced_cube_filename in all_datacube_files:

        cube_parameters = {}
        header = fits.getheader(reduced_cube_filename)

        try:
            target = header['OBJECT']  # Obtain the target name
        except KeyError:
            logging.warning(str(reduced_cube_filename) + " does not have a header key OBJECT, skipping file.")
            print("The header OBJECT does not exist")
            continue
        if Target.exists(target_name=target):  # If the target already exists in the database, get it from the db
            cube_parameters['target'] = Target.get(target_name=target)
        else:
            cube_parameters['target'] = Target(target_name=target)  # Else, create the target

        try:
            cube_parameters['observation_time'] = fetch_data(header, 'DATE-OBS')
            if not cube_parameters['observation_time'] in unique_observations:
                unique_observations.append(cube_parameters['observation_time'])
        except KeyError:
            logging.warning(str(reduced_cube_filename) + " the header DATE-OBS does not exist")
            print("The header DATE-OBS does not exist")
            continue
        except Exception as e:
            logging.warning(str(reduced_cube_filename) + " " + str(e))
            continue

        cube_parameters['instrument_mode'] = fetch_data(header, 'HIERARCH ESO INS MODE')
        cube_parameters['header'] = dict(header)
        try:
            del cube_parameters['header']['COMMENT']  # and then delete the COMMENT key, because it does not have the JSON format and
        except KeyError:  # throws an error when saving it in the database
            pass

        reduced_cube_entries.append(cube_parameters)

    # PSF extraction
    psf_entries = []

    print("Parsing through found prm files...")
    logging.info("Starting prm and psf file parsing.")
    for prm_file in all_prm_files:
        # PRM File
        observation_info = {}
        psfParams = {}  # Dictionary to store the extensions of the prm files
        sources = {}  # Dictionary to store the extensions (each source) of the psf files
        try:
            with fits.open(prm_file) as hduList:
                prm_filepath, prm_filename = os.path.split(prm_file)  # split to get the single file name
                # Notice that single file has the same root name that prm file
                expected_cube_name = prm_filename.replace('.prm', '')
                corresponding_cube_fullpath = prm_filepath + "/" + expected_cube_name
                expected_psf_name = prm_filename.replace('.prm', '.psf')
                corresponding_psf_fullpath = prm_filepath + "/" + expected_psf_name

                observation_info['observation_time'] = 'n/a'
                observation_info['target'] = 'n/a'
                observation_info['instrument_mode'] = 'n/a'
                # Check if the datacube is in the same directory, otherwise look elsewhere for it
                if Path(corresponding_cube_fullpath).exists():
                    observation_info['observation_time'] = fits.getheader(corresponding_cube_fullpath)['DATE-OBS']
                    observation_info['target'] = Target.get(target_name=fits.getheader(corresponding_cube_fullpath)['OBJECT'])
                    observation_info['instrument_mode'] = fits.getheader(corresponding_cube_fullpath)['HIERARCH ESO INS MODE']
                else:
                    for datacube_file in all_datacube_files:
                        if datacube_file.name == expected_cube_name:
                            observation_info['observation_time'] = fits.getheader(datacube_file)['DATE-OBS']
                            observation_info['target'] = Target.get(target_name = fits.getheader(datacube_file)['OBJECT'])
                            observation_info['instrument_mode'] = fits.getheader(datacube_file)['HIERARCH ESO INS MODE']
                            break

                if observation_info['target'] == 'n/a':
                    # PSF without target or observation parameters is kind of useless, skip these
                    logging.warning((prm_file + " Could not find a matching cube for prm file. Skipping file."))
                    print("Couldn't find a match for a prm file " + prm_filename + ".")
                    continue

                try:
                    params = hduList['PSFPARS'].data  # If the prm file does not have the PSFPARS extension, skip
                except KeyError:
                    logging.warning(prm_file + " Could not find PSFPARS in header. Possibly a broken prm file, skipping file.")
                    continue
                try:
                    for i in range(len(params['name'])):  # The column 'name' is a list with the parameter names
                        parameter = params['name'][i]
                        value = params['polyfit'][i].tolist()  # The column 'polyfit' is a list of lists
                        psfParams[parameter] = value  # Store the polyfit of the PSF parameters
                    prmStep = hduList['SPECTRA'].header['CDELT3']
                    prmRestw = hduList['SPECTRA'].header['CRVAL3']
                    prmData = np.arange(hduList['SPECTRA'].header['NAXIS3'])  # Calculate the wavelength
                    prmWavelength = (prmRestw + (prmData * prmStep)) * 10 ** 9
                    psfParams['wavelength'] = prmWavelength.tolist()
                except Exception as e:
                    print(str(prm_file) + " " + str(e))
                    logging.warning(str(prm_file) + " " + str(e))
                    continue
        except FileNotFoundError:
            logging.warning("The file " + str(prm_filename) + " does not exist.")
            print(f"The file {prm_filename} does not exist\n")  # If the prm file does not exists, skip
            continue
        except Exception as e:
            logging.warning(str(prm_filename) + " " + str(e))
            continue

            # PSF File
        try:
            with fits.open(corresponding_psf_fullpath) as hduList:
                for tupla in hduList.info(False):  # Iterate in the list of the extensions
                    if ('PM_' in tupla[1]):  # Looking for the 'PM_X' extensions, where X is the id of the source
                        sourceID = tupla[1]
                        sourceData = {}  # Dictionary that will store the data of the source
                        table = hduList[sourceID].data  # Access to the source extension
                        i = 0
                        for column in table.columns.names:  # Iterate in the list of column names
                            sourceData[column] = table[column].tolist()  # Store the column with its list of values
                            i += 1
                        sources[sourceID] = sourceData  # Store the source data
                if sources == {}:
                    print("Possibly a bad psf file, skipping file.")
                    logging.warning(corresponding_psf_fullpath + " Possibly a bad psf file, PM_* layers are missing")
        except FileNotFoundError:
            logging.warning("The file " + str(corresponding_psf_fullpath) + "does not exist.")
            print(f"The file {corresponding_psf_fullpath} does not exist\n")  # If the psf file does not exists, skip
            continue
        except Exception as e:
            logging.warning(str(prm_filename) + " " + str(e))
            continue

        psf_entries.append((observation_info.copy(), psfParams.copy(), sources.copy()))
    logging.info("Finished prm and psf file parsing.")

    ##############################
    # Catalog parsing starts here
    ##############################

    catalog_entries = []
    print("Parsing through found source catalog files...")
    logging.info("Starting pampelmuse source catalog file parsing.")
    for catalog_file in all_catalog_files:
        observation_info = {}

        try:
            catalog_filepath, catalog_filename = os.path.split(catalog_file)

            expected_cube_name = catalog_filename.replace('.detections.cat', '.fits')
            corresponding_cube_fullpath = catalog_filepath + "/" + expected_cube_name

            observation_info['observation_time'] = 'n/a'
            observation_info['target'] = 'n/a'
            observation_info['instrument_mode'] = 'n/a'
            # Check if the datacube is in the same directory, otherwise look elsewhere for it
            if Path(corresponding_cube_fullpath).exists():
                observation_info['observation_time'] = fits.getheader(corresponding_cube_fullpath)['DATE-OBS']
                observation_info['target'] = Target.get(target_name=fits.getheader(corresponding_cube_fullpath)['OBJECT'])
                observation_info['instrument_mode'] = fits.getheader(corresponding_cube_fullpath)['HIERARCH ESO INS MODE']
            else:
                for datacube_file in all_datacube_files:
                    if datacube_file.name == expected_cube_name:
                        observation_info['observation_time'] = fits.getheader(datacube_file)['DATE-OBS']
                        observation_info['target'] = Target.get(target_name=fits.getheader(datacube_file)['OBJECT'])
                        observation_info['instrument_mode'] = fits.getheader(datacube_file)['HIERARCH ESO INS MODE']
                        break

            if observation_info['target'] == 'n/a':
                # Catalog without target or observation parameters is kind of useless, skip these
                print("Couldn't find a match for a catalog file " + catalog_filename + ".")
                continue

            catalog_table = Table.read(catalog_file, format="ascii.ecsv")
            catalog_array = np.asarray(catalog_table)
            list_of_cat_entries = []
            for i in range(catalog_array.shape[0]):
                entry_dict = {k: v for k, v in zip(catalog_array[i].dtype.names, catalog_array[i])}
                list_of_cat_entries.append(entry_dict)

            catalog_entries.append((observation_info.copy(), list_of_cat_entries.copy()))
        except Exception as e:
            logging.warning(str(catalog_file) + " " + str(e))
            continue

    #######################################
    # Raw exposure extraction starts here
    #######################################

    print("Starting raw exposure extraction")
    logging.warning("Starting raw exposure extraction.")
    raw_fits_entries = []
    for raw_filename in all_raw_files:

        raw_parameters = {}
        header = fits.getheader(raw_filename)

        try:
            target = fetch_data(header, 'OBJECT')  # Obtain the target name
        except KeyError:
            logging.warning(str(raw_filename) + " does not have a header key OBJECT, skipping file.")
            print("The header OBJECT does not exist")
            continue
        if Target.exists(target_name=target):  # If the target already exists in the database, get it from the db
            raw_parameters['target'] = Target.get(target_name=target)
        else:
            raw_parameters['target'] = Target(target_name=target)  # Else, create the target

        try:
            raw_parameters['observation_time'] = header['DATE-OBS']
            if not raw_parameters['observation_time'] in unique_observations:
                unique_observations.append(raw_parameters['observation_time'])
        except KeyError:
            logging.warning(str(raw_filename) + " does not have a header key DATE-OBS, skipping file.")
            print("The header DATE-OBS does not exist")
            continue
        except Exception as e:
            logging.warning(str(raw_filename) + " " + str(e))
            continue

        raw_parameters['instrument_mode'] = fetch_data(header, 'HIERARCH ESO INS MODE')
        raw_parameters['header'] = dict(header)
        try:
            del raw_parameters['header']['COMMENT']  # and then delete the COMMENT key, because it does not have the JSON format and
        except KeyError:  # throws an error when saving it in the database
            pass
        try:
            del raw_parameters['header']['']  # and then delete the '' key, because it does not have the JSON format and
        except KeyError:  # throws an error when saving it in the database
            pass

        raw_fits_entries.append(raw_parameters)
    logging.warning("Finished raw exposure extraction.")

    nightlog_weather_entries = []
    corresponding_rawfile = ""
    # Nightlog extraction
    print("Starting nightlog extraction")
    logging.info("Starting nightlog extraction.")
    for nightlog_filename in all_nightlog_files:
        try:
            nightlog_parameters = {}
            corresponding_rawfile = str(nightlog_filename)[:-6] + "fits.fz"
            raw_header = fits.getheader(corresponding_rawfile)
        except FileNotFoundError:
            logging.warning("Raw file " + corresponding_rawfile + " not found for a nightlog.")
            print("Raw file " + corresponding_rawfile + " not found for a nightlog.")
            continue
        except Exception as e:
            logging.warning(corresponding_rawfile + " " + str(e))
            continue
        try:
            with open(nightlog_filename) as f:
                read_lines = f.readlines()
        except FileExistsError:
            logging.warning("Nightlog file " + nightlog_filename + " not found.")
            print("Nightlog file " + nightlog_filename + " not found.")
            continue
        # Redundant doublecheck that this is indeed a science exposure
        if not ('HIERARCH ESO DPR CATG' in header and header['HIERARCH ESO DPR CATG'] == 'SCIENCE'):
            continue

        obs_time = header["DATE-OBS"]
        try:
            local_start_time = float(raw_header["UTC"])
            integration_time = float(raw_header["EXPTIME"])
        except Exception as e:
            logging.warning(corresponding_rawfile + " " + str(e))
            continue

        local_end_time = (local_start_time + integration_time) % 86400.0
        if local_start_time > 43200.0:
            local_start_time -= 86400.0

        line_extraction = []

        # Read through the nightlog, starting from the weather report
        # Take the hh:mm time and transform into UTC seconds
        for line in read_lines[10:]:
            if line == "---------------------------------------------------\n" or line == "\n":
                break
            line_split = line[:-1].split("\t")
            try:
                weather_time = float(line_split[0].split(":")[0]) * 3600.0 + float(line_split[0].split(":")[1]) * 60.0
            except ValueError:
                continue
            except Exception as e:
                logging.warning(nightlog_filename + " " + str(e))
                continue

            if weather_time > 43200.0:
                weather_time -= 86400.0
            #print(weather_time)
            weather_line = [weather_time,
                            line_split[1],
                            line_split[2]]
            line_extraction.append(weather_line)

        start_weather_time = 0.0
        end_weather_time = 0.0
        observation_start_condition = "None"
        observation_end_condition = "None"
        observation_start_comment = "None"
        observation_end_comment = "None"
        i = -1
        # Loop through the weather comments
        try:
            for line in line_extraction:
                weather_comment_time = line[0]
                i += 1
                if weather_comment_time > 43200:
                    weather_comment_time -= 86400
                # Check we've reached past the exposure start
                if weather_comment_time > local_start_time:
                    break
                start_weather_time = line[0]
                end_weather_time = line[0]
                observation_start_condition = line[1]
                observation_end_condition = line[1]
                observation_start_comment = line[2]
                observation_end_comment = line[2]

            # Loop through the weather comments, but start where we left off previously
            for line in line_extraction[i:]:
                weather_comment_time = line[0]
                if line[0] > 43200:
                    weather_comment_time -= 86400
                # Check we've reached past the exposure end
                if weather_comment_time > local_end_time:
                    break
                end_weather_time = line[0]
                observation_end_condition = line[1]
                observation_end_comment = line[2]
        except Exception as e:
            logging.warning(nightlog_filename + " " + str(e))
            continue

        # Add all the interesting parameters to a dictionary...
        nightlog_parameters["observation_time"] = raw_header["DATE-OBS"]
        nightlog_parameters["sky_condition_start_time"] = start_weather_time
        nightlog_parameters["sky_condition_start"] = observation_start_condition
        nightlog_parameters["sky_comment_start"] = observation_start_comment
        nightlog_parameters["sky_condition_end_time"] = end_weather_time
        nightlog_parameters["sky_condition_end"] = observation_end_condition
        nightlog_parameters["sky_comment_end"] = observation_end_comment

        # ...and throw it into the pile
        nightlog_weather_entries.append(nightlog_parameters)
    logging.info("Finished nightlog extraction.")

    # Raman file extraction
    raman_fits_entries = []
    print("Starting Raman file extraction")
    logging.info("Starting Raman file extraction.")
    for raman_filename in all_raman_files:

        raman_parameters = {}
        header = fits.getheader(raman_filename)

        try:
            target = header['OBJECT']  # Obtain the target name
        except KeyError:
            logging.warning(str(raman_filename + " does not have a header key OBJECT, skipping file."))
            print("The header OBJECT does not exist")
            continue
        if Target.exists(target_name=target):  # If the target already exists in the database, get it from the db
            raman_parameters['target'] = Target.get(target_name=target)
        else:
            raman_parameters['target'] = Target(target_name=target)  # Else, create the target

        try:
            raman_parameters['observation_time'] = fetch_data(header, 'DATE-OBS')
            if not raman_parameters['observation_time'] in unique_observations:
                unique_observations.append(raman_parameters['observation_time'])
        except KeyError:
            logging.warning(str(reduced_cube_filename) + " does not have a header key DATE-OBS, skipping file.")
            print("The header DATE-OBS does not exist")
            continue

        raman_parameters['instrument_mode'] = fetch_data(header, 'HIERARCH ESO INS MODE')
        raman_parameters['header'] = dict(header)
        try:
            del raman_parameters['header']['COMMENT']  # and then delete the COMMENT key, because it does not have the JSON format and
        except KeyError:  # throws an error when saving it in the database
            pass

        raman_fits_entries.append(raman_parameters)
    logging.info("Finished Raman file extraction.")

    # Identification and linking of previous files to a same exposure
    print("Starting the file linking and database update procedure.")
    logging.info("Starting the file linking and database update procedure.")
    new_entries = 0
    modified_entries = 0
    for observation_key in unique_observations:
        observation_dictionary = {}
        for reduced_cube_entry in reduced_cube_entries:
            if observation_key in reduced_cube_entry['observation_time']:
                observation_dictionary['target'] = reduced_cube_entry['target']
                observation_dictionary['observation_time'] = reduced_cube_entry['observation_time']
                observation_dictionary['insMode'] = reduced_cube_entry['instrument_mode']
                observation_dictionary['datacube_header'] = reduced_cube_entry['header']
                reduced_cube_entries.remove(reduced_cube_entry)
                break
        for raw_fits_entry in raw_fits_entries:
            if observation_key in raw_fits_entry['observation_time']:
                observation_dictionary['target'] = raw_fits_entry['target']
                observation_dictionary['observation_time'] = raw_fits_entry['observation_time']
                observation_dictionary['insMode'] = raw_fits_entry['instrument_mode']
                observation_dictionary['raw_exposure_header'] = raw_fits_entry['header']
                raw_fits_entries.remove(raw_fits_entry)
                break
        for nightlog_entry in nightlog_weather_entries:
            if observation_key in nightlog_entry['observation_time']:
                observation_dictionary["sky_condition_start_time"] = nightlog_entry["sky_condition_start_time"]
                observation_dictionary["sky_condition_start"] = nightlog_entry["sky_condition_start"]
                observation_dictionary["sky_comment_start"] = nightlog_entry["sky_comment_start"]
                observation_dictionary["sky_condition_end_time"] = nightlog_entry["sky_condition_end_time"]
                observation_dictionary["sky_condition_end"] = nightlog_entry["sky_condition_end"]
                observation_dictionary["sky_comment_end"] = nightlog_entry["sky_comment_end"]
                nightlog_weather_entries.remove(nightlog_entry)
                break
        for raman_fits_entry in raman_fits_entries:
            if observation_key in raman_fits_entry['observation_time']:
                observation_dictionary['target'] = raman_fits_entry['target']
                observation_dictionary['observation_time'] = raman_fits_entry['observation_time']
                observation_dictionary['insMode'] = raman_fits_entry['instrument_mode']
                observation_dictionary['raman_image_header'] = raman_fits_entry['header']
                raman_fits_entries.remove(raman_fits_entry)
                break
        for psf_entry in psf_entries:
            if observation_key in psf_entry[0]['observation_time']:
                observation_dictionary['target'] = psf_entry[0]['target']
                observation_dictionary['observation_time'] = psf_entry[0]['observation_time']
                observation_dictionary['insMode'] = psf_entry[0]['instrument_mode']
                observation_dictionary['psf_params'] = psf_entry[1]
                observation_dictionary['sources'] = psf_entry[2]
                psf_entries.remove(psf_entry)
                break
        for catalog_entry in catalog_entries:
            if observation_key in catalog_entry[0]['observation_time']:
                observation_dictionary['target'] = catalog_entry[0]['target']
                observation_dictionary['observation_time'] = catalog_entry[0]['observation_time']
                observation_dictionary['insMode'] = catalog_entry[0]['instrument_mode']
                observation_dictionary['pampelmuse_catalog'] = catalog_entry[1]
                catalog_entries.remove(catalog_entry)
                break

        # Modify headers into more JSON-like format
        if 'pampelmuse_catalog' in observation_dictionary:
            observation_dictionary['pampelmuse_catalog'] = json.dumps(observation_dictionary['pampelmuse_catalog'], default=convert_npint64_to_int)
        if 'psf_params' in observation_dictionary:
            observation_dictionary['psf_params'] = json.dumps(observation_dictionary['psf_params'])
        if 'sources' in observation_dictionary:
            observation_dictionary['sources'] = json.dumps(observation_dictionary['sources'])
        if 'raman_image_header' in observation_dictionary:
            observation_dictionary['raman_image_header'] = json.dumps(observation_dictionary['raman_image_header'])
        if 'raw_exposure_header' in observation_dictionary:
            observation_dictionary['raw_exposure_header'] = json.dumps(observation_dictionary['raw_exposure_header'])
        if 'datacube_header' in observation_dictionary:
            observation_dictionary['datacube_header'] = json.dumps(observation_dictionary['datacube_header'])

        # Check if exposure exists, create a new exposure if not and modify existing one if needed
        if not Exposure.exists(observation_time=observation_dictionary['observation_time']):
            Exposure(**observation_dictionary)
            new_entries += 1
        else:
            entry_is_modified = False
            modified_entry_dictionary = {}
            database_entry = Exposure.get(observation_time=observation_dictionary['observation_time'])
            for column, value in observation_dictionary.items():
                # Ignore columns that shouldn't change for an observation
                if (column == 'observation_time' or
                        column == 'target' or
                        column == 'insMode'):
                    continue
                database_value = getattr(database_entry, column)
                if type(value) == str:
                    value = value.strip()
                    database_value = database_value.strip()
                if not value == database_value:
                    modified_entry_dictionary[column] = value
                    entry_is_modified = True
            if entry_is_modified:
                logging.info("Modifying an entry with observation time " + observation_dictionary['observation_time'])
                database_entry.set(**modified_entry_dictionary)
                modified_entries += 1

    print('Finished updating the database, ' + '{:d}'.format(new_entries) + " new entries and "
          + '{:d}'.format(modified_entries) + " modified entries.")
    logging.info('Finished updating the database, ' + '{:d}'.format(new_entries) + " new entries and "
          + '{:d}'.format(modified_entries) + " modified entries.")


# Main part starts here
log_filename = datetime.now().strftime('%d-%m-%Y_%H:%M:%S.log')
logging.basicConfig(filename=log_filename,
                    format='%(asctime)s %(levelname)s %(message)s',
                    datefmt='%d-%m-%Y %H:%M:%S.log',
                    level=logging.INFO)
print("Started a logfile called " + log_filename)
try:
    db.bind(provider="mysql", host="127.0.0.1", user=user, passwd=passwd, db=database_name)
    db.generate_mapping(check_tables=False, create_tables=True)
    logging.info("Connected to SQL database successfully.")
except Exception as e:
    print("Error with SQL connection, check the log for more information.")
    logging.error("SQL connection problem")
    logging.error(e)
    quit(2)

muse_script()
