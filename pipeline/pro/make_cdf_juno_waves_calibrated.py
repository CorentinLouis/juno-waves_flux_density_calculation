#! /usr/bin/python
# -*- coding: utf-8 -*-


_CDFCHECK_BIN = '/Applications/cdf/pds-cdf-1.0.11/bin/cdfcheck'
_CDFLIB_BIN = '/Applications/cdf/cdf/bin/'
_SKTEDITOR_BIN = '/Applications/cdf/skteditor-1.3.1.41/spdfjavaClasses.jar'

#_CDFCHECK_BIN = '/opt/cdf/pds-cdf-1.0.11/bin/cdfcheck'
#_CDFLIB_BIN = '/opt/cdf/cdf38_1-dist/bin/'
#_SKTEDITOR_BIN = '/opt/cdf/skteditor-1.3.1.41/spdfjavaClasses.jar'

from scipy.io import readsav
import numpy
import os
import argparse
from astropy.units import Unit as u
from astropy.time import Time
from datetime import datetime, time, timedelta
import json
import requests
from spacepy import pycdf

def fix_cdf(filename, verbose=True):
    """
    """
    if verbose:
        print('\n\tLOG: CDF fixing')
    if not os.path.isfile(_CDFCHECK_BIN):
        print('WARNING: _CDFCHECK_BIN is not set properly')
        return

    import subprocess
    tmp_filename = filename[:-4] + '_tmp.cdf'

    cur_command = [
        os.path.join(_CDFLIB_BIN, 'cdfconvert'),
        filename,
        tmp_filename
        ]
    cur_p = subprocess.Popen(cur_command, stdout=subprocess.PIPE)
    log_str = cur_p.communicate()[0].decode('ascii')

    os.rename(tmp_filename, filename)
    log_str += 'Moving converted "{}" back to "{}"\n'.format(tmp_filename, filename)

    if verbose:
        print(log_str)
    return

def add_validate(filename, val_dict,verbose=True):
    """
    """
    if verbose:
        print('\n\tLOG: Adding CDF Validate info')

    if isinstance(val_dict, dict):
        val_string = "{}>{}>{}>{}".format(val_dict['name'], val_dict['result'], val_dict['agent'], val_dict['when'])
        from spacepy import pycdf
        c = pycdf.CDF(filename)
        c.readonly(False)
        if 'Validate' in c.attrs.keys():
            tmp_val = c.attrs['Validate'][:]
            tmp_val.append(val_string)
            c.attrs['Validate'] = tmp_val
        else:
            c.attrs['Validate'] = [val_string]
        c.close()
        fix_cdf(filename, verbose=verbose)



def check_cdf(filename, verbose=True):
    """
    """
    if verbose:
        print('\n\tLOG: CDF checking')
    if not os.path.isfile(_CDFCHECK_BIN):
        print('WARNING: _CDFCHECK_BIN is not set properly')
        return
    if not os.path.isfile(_SKTEDITOR_BIN):
        print('WARNING: _SKTEDITOR_BIN is not set properly')
        return
    if not os.path.isdir(_CDFLIB_BIN):
        print('WARNING: _CDFLIB_BIN is not set properly')
        return

    import subprocess
    import socket

    fix_cdf(filename, verbose=verbose)

    verb_flag = ''
    log_str = ""

    if verbose:
        verb_flag = '-v'

    cur_command = ['/bin/bash', str(_CDFCHECK_BIN), verb_flag, filename]
    log_str += ' '.join(cur_command)+'\n'
    cur_p = subprocess.Popen(cur_command, stdout=subprocess.PIPE)
    tmp_log = cur_p.communicate()[0].decode('ascii')
    assert isinstance(tmp_log, str)
    if 'OK' in tmp_log:
        tmp_test = 'pass'
    else:
        tmp_test = 'fail'
    add_validate(filename, {
        'name':"PDS-CDF",
        'result':tmp_test,
        'agent':"{}:{}".format(socket.gethostname().split('.')[0],__file__.split('/')[-1]),
        'when':numpy.datetime64('now').astype(str).split()[0]
    }, verbose=verbose)
    log_str += tmp_log
    print(log_str)

    log_str = ""
    cur_command =[
        'java',
        '-cp',
        '{}:{}/../cdfjava/classes/cdfjava.jar'.format(_SKTEDITOR_BIN, _CDFLIB_BIN),
        'gsfc.spdf.istp.tools.CDFCheck',
        filename]
    log_str += ' '.join(cur_command)+'\n'
    cur_p = subprocess.Popen(cur_command, stdout=subprocess.PIPE)
    log_str += cur_p.communicate()[0].decode('ascii')
    if 'FAILED' in log_str:
        tmp_test = 'fail'
    else:
        tmp_test = 'pass'
    add_validate(filename, {
    'name':"CDF-ISTP",
    'result':tmp_test,
    'agent':"{}:{}".format(socket.gethostname().split('.')[0],__file__.split('/')[-1]),
    'when':numpy.datetime64('now').astype(str).split()[0]
    }, verbose=verbose)

    print(log_str)
    return


def juno_lesia_calibrated_to_cdf(filename, output_path="../../data_n2/spdyn_cdf_data_v02/", data_version=None, verbose=False):
    if data_version == None:
        print("You need to specify a data version")
        return
    elif data_version == 1:
        data_version_name = ''
    else:
        data_version_name = f'_v{data_version:02d}'

    result0=readsav(f"../gain/gain_final_lin_v{data_version:02d}.sav")
    gain=result0['gain_final'].copy()
    # During gain calculation, we only calculating the gain >1 khZ (i.e. channels 16-126)
    # because we use Cassini+Voyager observation (>1 kHz) to obtain the gain
    # For Version =>2 we want to keep the full frequency ranges.
    # Therefore channel [0-15] are calibrated using the same gain than channel=16
    if data_version >= 2:
        gain_under_channel16=numpy.zeros(16)+gain[0]
        gain = numpy.concatenate((gain_under_channel16,gain))
    result1=readsav(f"../background/background_ALLPJ_zlin_v{data_version:02d}.sav")
    background=(result1["background"].copy())*gain
    sigma=(result1["sigma"].copy())*gain
    result=readsav(filename)
    yyyyddd=int(result['yyyyddd'])
    current_day=datetime.strptime(str(yyyyddd),'%Y%j')
    yyyymmdd=current_day.strftime('%Y%m%d')
    t=result['t']
    f=result['f']*u('kHz')
    data=(result['zlincal'].copy()).transpose()
    for ifreq in range (f.shape[0]):
        data[:,ifreq]=data[:,ifreq]+background[ifreq]
    # data[data == 0.0]=-1.0
    # data[data == 0] = data[data >0].min()
    # To be changed?
    # data == 0.0 --> -1
    epoch=numpy.array([current_day+timedelta(seconds=s) for s in t])


    # Checking the existence of a directory for the year
    year_dir=str(yyyyddd)[0:4]
    month_dir=str(yyyymmdd)[4:-2]
    path_dir=output_path+year_dir
    CHECK_FOLDER = os.path.exists(path_dir)
    if not CHECK_FOLDER:
        os.makedirs(path_dir)
        if verbose:
            print("Created folder : ",path_dir)
    path_dir=path_dir+"/"+month_dir
    CHECK_FOLDER = os.path.exists(path_dir)
    if not CHECK_FOLDER:
        os.makedirs(path_dir)
        if verbose:
            print("Created folder : ",path_dir)

    cdfname=f"jno_wav_cdr_lesia_{yyyymmdd}_v{data_version:02d}"
    # next step: ask for permission to remove existing file. For now: remove by default
    if os.path.exists(path_dir+"/"+cdfname+".cdf"):
        os.remove(path_dir+"/"+cdfname+".cdf")

    # Opening CDF object
    pycdf.lib.set_backward(False)  # this is setting the CDF version to be used

    cdf = pycdf.CDF(path_dir+"/"+cdfname+".cdf", '')

    # required settings for ISTP and PDS compliance
    cdf.col_major(True)  # Column Major
    cdf.compress(pycdf.const.NO_COMPRESSION)  # No file level compression

    # SETTING ISTP GLOBAL ATTRIBUTES
    cdf.attrs['TITLE'] = 'LESIA calibrated Juno/Waves data'
    cdf.attrs['Project'] = [
        'OBSPM>Observatoire de Paris',
        'PADC>Paris Astronomical Data Centre'
    ]
    cdf.attrs['Discipline'] = 'Planetary Physics>Waves'
    cdf.attrs['Data_type'] = 'cdr>Calibrated Data Record'
    # to be changed?
    cdf.attrs['Descriptor'] = 'lesia>LESIA calibration'
    cdf.attrs['Data_version'] = str(data_version)
    cdf.attrs['Instrument_type'] = 'Radio and Plasma Waves (space)'
    # to be changed?
    cdf.attrs['Logical_file_id'] = cdfname
    # to be changed?
    cdf.attrs['Logical_source'] = 'jno_wav_cdr_lesia'
    # to be changed?
    cdf.attrs['Logical_source_description'] = 'LESIA calibrated Juno/Waves data'
    cdf.attrs['File_naming_convention'] = 'source_type_descriptor_yyyymmdd_ver'
    cdf.attrs['Mission_group'] = 'Juno'
    # to be changed?
    cdf.attrs['PI_name'] = 'W. S. Kurth'
    # to be changed?
    cdf.attrs['PI_affiliation'] = [
        (
            "Department of Physics and Astronomy, University of Iowa, Iowa City, IA, USA"
        )
    ]
    cdf.attrs['Source_name'] = 'jno_wav>Juno Waves'
    # to be changed?
    cdf.attrs['TEXT'] = (
        "These data are calibrated from the Juno/Waves raw data "
        "Louis, C. K., Zarka, P., Dabidin, K., Lampson, P.-A., Magalhaes, F. P., Boudouma, A., et al. (2021),"
        "Latitudinal beaming of Jupiter's radio emissions from Juno/Waves flux density measurements."
        "Journal of Geophysical Research: Space Physics, 126, e2021JA029435."
        "https://doi.org/10.1029/2021JA029435."
    )
    cdf.attrs['Generated_by'] = [
        "LESIA>Laboratoire d'Etudes Spatiales et d'Instrumentation en Astrophysique"
    ]
    cdf.attrs['Generation_date'] = '{}'.format(
        numpy.datetime64('now').astype(str).split('T')[0].replace('-', '')
    )
    cdf.attrs['DOI'] = "https://doi.org/10.25935/6jg4-mk86"
    cdf.attrs['LINK_TEXT'] = "Juno/Waves calibrated data are available at"
    cdf.attrs['LINK_TITLE'] = 'PADC/MASER repository'
    cdf.attrs['HTTP_LINK'] = 'https://maser.obspm.fr/repository/juno/waves/data/'
    cdf.attrs['MODS'] = [
        (
            f"v{data_version:02d}: First release."
        )
    ]
    cdf.attrs['Parents'] = os.path.basename(filename)
    # to be changed?
    cdf.attrs['Rules_of_use'] = [
        (
            'Access is restricted to the Juno team.'
        ),
        'Contact email: contact.maser@obspm.fr',
        (
            'We kindly request the authors of any communications '
            'and publications using these data to let us know '
            'about them, include minimal citation to the reference '
            'below and appropriate acknowledgements whenever needed.'
        ),
        'References:',
        (
            "Louis, C. K., Zarka, P., Dabidin, K., Lampson, P.-A., Magalhaes, F. P., Boudouma, A., et al. (2021),"
            "Latitudinal beaming of Jupiter's radio emissions from Juno/Waves flux density measurements."
            "Journal of Geophysical Research: Space Physics, 126, e2021JA029435."
            "https://doi.org/10.1029/2021JA029435."
        ),
        'Acknowledgements: see the acknowledgement field.'
    ]
    cdf.attrs['Skeleton_version'] = str(data_version)
    #cdf.attrs['Software_version'] = str(data_version)
    #cdf.attrs['Software_language'] = 'python3'
    cdf.attrs['Time_resolution'] = '1 second'
    # to be changed?
    cdf.attrs['Acknowledgement'] = (
        "The authors acknowledge the Observatoire de Paris, CNES, CNRS "
        "for funding and supporting this work "
        "and the University of Iowa and the Juno/Waves team for providing access "
        "to the Juno/Waves data accessible online from PDS at https://doi.org/10.17189/1519708"
    )

    # SPASE
    # to be changed?
    #cdf.attrs['spase_DatasetResourceID'] = 'spase://NASA/NumericalData/Juno/Waves/cdr_lesia/PT1S'
    cdf.attrs['spase_DatasetResourceID'] = 'spase://NASA/NumericalData/Juno/Waves/lesia_l3a/PT1S'

    # SETTING PDS GLOBAL ATTRIBUTES
    cdf.attrs['PDS_Observation_start_time'] = epoch[0].strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3] + 'Z'
    cdf.attrs['PDS_Observation_stop_time'] = epoch[-1].strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3] + 'Z'
    cdf.attrs['PDS_Observation_target'] = 'Jupiter'
    cdf.attrs['PDS_Observation_type'] = 'Radio'



    # SETTING VESPA GLOBAL ATTRIBUTES
    cdf.attrs['VESPA_dataproduct_type'] = 'DS>Dynamic Spectrum'
    cdf.attrs['VESPA_target_class'] = 'Planet'
    cdf.attrs['VESPA_target_region'] = 'Magnetosphere'
    cdf.attrs['VESPA_feature_name'] = 'Radio emission'
    cdf.attrs['VESPA_instrument_host_name'] = 'Jno>Juno'
    cdf.attrs['VESPA_instrument_name'] = 'Wav>Waves'
    cdf.attrs['VESPA_receiver_name'] = ['LFR>Low Frequency Receiver','HFR>High Frequency Receiver']
    cdf.attrs['VESPA_measurement_type'] = 'phys.flux.density;em.radio;phys.polarization'
    cdf.attrs['VESPA_access_format'] = 'application/x-cdf'
    cdf.attrs['VESPA_bib_reference'] = '10.1029/2021JA029435'

    # sampling step
    cdf.attrs['VESPA_time_sampling_step'] = '1.0'
    cdf.attrs['VESPA_time_sampling_step_unit'] = 's'
    # integration time
    # cdf.attrs['VESPA_time_exp'] =
    # cdf.attrs['VESPA_time_exp_unit'] = 'ms'

    cdf.attrs['VESPA_spectral_range_min'] = str(f[0].value)
    cdf.attrs['VESPA_spectral_range_max'] = str(f[-1].value)
    cdf.attrs['VESPA_spectral_range_unit'] = str(f.unit)


    # sampling step
    # cdf.attrs['VESPA_spectral_sampling_step'] = '75'# * u.kHz
    # cdf.attrs['VESPA_spectral_sampling_step_unit'] = 'kHz'

    # resolution (=bandwidth) (new version resolution will bandwidth/frequency)
    # cdf.attrs['VESPA_spectral_resolution'] = '30'# * u.kHz
    # cdf.attrs['VESPA_spectral_resolution_unit'] = 'kHz'


    # SETTING VARIABLES
    cdf.new('Epoch',
        data=epoch,
        type=pycdf.const.CDF_TIME_TT2000,
        compress=pycdf.const.NO_COMPRESSION)
    cdf['Epoch'].attrs.new('VALIDMIN', data=datetime(2011,8,9,0,0), type=pycdf.const.CDF_TIME_TT2000)
    cdf['Epoch'].attrs.new('VALIDMAX', data=datetime(2025,10,1,0,0), type=pycdf.const.CDF_TIME_TT2000)
    cdf['Epoch'].attrs.new('SCALEMIN', data=current_day, type=pycdf.const.CDF_TIME_TT2000)
    cdf['Epoch'].attrs.new('SCALEMAX', data=current_day+timedelta(days=1), type=pycdf.const.CDF_TIME_TT2000)
    cdf['Epoch'].attrs['CATDESC'] = "SCET time of each spectra."
    cdf['Epoch'].attrs['FIELDNAM'] = "Epoch"
    cdf['Epoch'].attrs.new('FILLVAL', data=-9223372036854775808, type=pycdf.const.CDF_TIME_TT2000)
    cdf['Epoch'].attrs['LABLAXIS'] = "Epoch"
    cdf['Epoch'].attrs['UNITS'] = "ns"
    cdf['Epoch'].attrs['FORM_PTR'] = "CDF_TIME_TT2000"
    cdf['Epoch'].attrs['VAR_TYPE'] = "support_data"
    cdf['Epoch'].attrs['SCALETYP'] = "linear"
    cdf['Epoch'].attrs['MONOTON'] = "INCREASE"
    cdf['Epoch'].attrs['REFERENCE_POSITION'] = "Juno"
    cdf['Epoch'].attrs['SI_CONVERSION'] = "1.0e-9>s"
    cdf['Epoch'].attrs['UCD'] = "time.epoch"
    cdf['Epoch'].attrs['TIME_BASE'] = 'J2000'
    # Is SCET a UTC scale?
    cdf['Epoch'].attrs['TIME_SCALE'] = 'UTC'

    cdf.new('Frequency',
        data=f.value,
        type=pycdf.const.CDF_REAL4,
        compress=pycdf.const.NO_COMPRESSION,
        recVary=False)
    cdf['Frequency'].attrs['CATDESC'] = "Central frequency of each step of the spectral sweep." # Central frequency of each step of the spectral sweep.
    cdf['Frequency'].attrs['DICT_KEY'] = "frequency"
    cdf['Frequency'].attrs['FIELDNAM'] = 'Frequency'
    cdf['Frequency'].attrs.new('FILLVAL', data=-1.0e+31, type=pycdf.const.CDF_REAL4)
    cdf['Frequency'].attrs['FORMAT'] = "E12.4"
    cdf['Frequency'].attrs['LABLAXIS'] = 'Frequency'
    cdf['Frequency'].attrs['UNITS'] = str(f.unit)
    cdf['Frequency'].attrs.new('VALIDMIN', data=f[0].value, type=pycdf.const.CDF_REAL4)
    cdf['Frequency'].attrs.new('VALIDMAX', data=f[-1].value, type=pycdf.const.CDF_REAL4)
    cdf['Frequency'].attrs['VAR_TYPE'] = "support_data"
    cdf['Frequency'].attrs['SCALETYP'] = "log"
    cdf['Frequency'].attrs.new('SCALEMIN', data=f[0].value, type=pycdf.const.CDF_REAL4)
    cdf['Frequency'].attrs.new('SCALEMAX', data=f[-1].value, type=pycdf.const.CDF_REAL4)
    cdf['Frequency'].attrs['UCD'] = "em.freq"
    cdf['Frequency'].attrs['SI_CONVERSION'] = f"{((f.unit/u('Hz')).decompose()).scale}>Hz"

    cdf.new('Data',
        data=data,
        type=pycdf.const.CDF_REAL4,
        compress=pycdf.const.NO_COMPRESSION)
    cdf['Data'].attrs['CATDESC'] = "Calibrated flux density spectrogram measured by Juno/Waves"
    cdf['Data'].attrs['DEPEND_0'] = 'Epoch'
    cdf['Data'].attrs['DEPEND_1'] = 'Frequency'
    cdf['Data'].attrs['DICT_KEY'] = 'electric_field>power'
    cdf['Data'].attrs['DISPLAY_TYPE'] = 'spectrogram'
    cdf['Data'].attrs['FIELDNAM'] = 'Data'
    cdf['Data'].attrs.new('FILLVAL', data=-1.0e+31, type=pycdf.const.CDF_REAL4)
    cdf['Data'].attrs['FORMAT'] = 'E12.4'
    cdf['Data'].attrs['LABLAXIS'] = 'Calibrated flux density'
    cdf['Data'].attrs['UNITS'] = 'V**2 m**-2 Hz**-1'
    cdf['Data'].attrs.new('VALIDMIN', data=1e-25, type=pycdf.const.CDF_REAL4) # 2 decade below the minimal intensity measured in the time range [20160101-20191231]
    cdf['Data'].attrs.new('VALIDMAX', data=1e-1, type=pycdf.const.CDF_REAL4) # 2 decade above the maximal intensity measured in the time range [20160101-20191231]
    cdf['Data'].attrs['VAR_TYPE'] = "Data"
    cdf['Data'].attrs['SCALETYP'] = "log"
    cdf['Data'].attrs.new('SCALEMIN', data=data.min(), type=pycdf.const.CDF_REAL4)
    cdf['Data'].attrs.new('SCALEMAX', data=data.max(), type=pycdf.const.CDF_REAL4)
    cdf['Data'].attrs['UCD'] = "phys.flux.density;em.radio" #;phys.polarization.linear

    cdf.new('Gain',
        data=gain,
        type=pycdf.const.CDF_REAL4,
        compress=pycdf.const.NO_COMPRESSION,
        recVary=False)
    cdf['Gain'].attrs['CATDESC'] = "Calibration gain applied to the Juno Waves FFT filtered data"
    cdf['Gain'].attrs['DEPEND_1'] = 'Frequency'
    # to be changed?
    cdf['Gain'].attrs['DICT_KEY'] = 'power>calibration'
    cdf['Gain'].attrs['DISPLAY_TYPE'] = 'time_series'
    cdf['Gain'].attrs['FIELDNAM'] = 'Gain'
    cdf['Gain'].attrs.new('FILLVAL', data=-1.0e+31, type=pycdf.const.CDF_REAL4)
    cdf['Gain'].attrs['FORMAT'] = 'E12.4'
    cdf['Gain'].attrs['LABLAXIS'] = 'Calibration Gain'
    cdf['Gain'].attrs['UNITS'] = ' '
    cdf['Gain'].attrs.new('VALIDMIN', data=0.00001, type=pycdf.const.CDF_REAL4)
    cdf['Gain'].attrs.new('VALIDMAX', data=10000, type=pycdf.const.CDF_REAL4)
    cdf['Gain'].attrs['VAR_TYPE'] = "support_data"
    cdf['Gain'].attrs['SCALETYP'] = "log"
    cdf['Gain'].attrs.new('SCALEMIN', data=0.00001, type=pycdf.const.CDF_REAL4)
    cdf['Gain'].attrs.new('SCALEMAX', data=10000, type=pycdf.const.CDF_REAL4)
    cdf['Gain'].attrs['UCD'] = "phot.flux.density;obs.calib.flat"

    cdf.new('Background',
        data=background,
        type=pycdf.const.CDF_REAL4,
        compress=pycdf.const.NO_COMPRESSION,
        recVary=False)
    cdf['Background'].attrs['CATDESC'] = "Background calculated from 2016-04-09 to 20190623, multiplied by the Calibraton Gain"
    cdf['Background'].attrs['DEPEND_1'] = 'Frequency'
    cdf['Background'].attrs['DICT_KEY'] = 'power>calibration'
    cdf['Background'].attrs['DISPLAY_TYPE'] = 'time_series'
    cdf['Background'].attrs['FIELDNAM'] = 'Background'
    cdf['Background'].attrs.new('FILLVAL', data=-1.0e+31, type=pycdf.const.CDF_REAL4)
    cdf['Background'].attrs['FORMAT'] = 'E12.4'
    cdf['Background'].attrs['LABLAXIS'] = 'Background x Gain'
    cdf['Background'].attrs['UNITS'] = 'V**2 m**-2 Hz**-1'
    cdf['Background'].attrs.new('VALIDMIN', data=1e-18, type=pycdf.const.CDF_REAL4)
    cdf['Background'].attrs.new('VALIDMAX', data=1e-5, type=pycdf.const.CDF_REAL4)
    cdf['Background'].attrs['VAR_TYPE'] = "support_data"
    cdf['Background'].attrs['SCALETYP'] = "log"
    cdf['Background'].attrs.new('SCALEMIN', data=1e-18, type=pycdf.const.CDF_REAL4)
    cdf['Background'].attrs.new('SCALEMAX', data=1e-5, type=pycdf.const.CDF_REAL4)
    cdf['Background'].attrs['UCD'] = "phot.flux.density;obs.calib"

    cdf.new('Sigma',
        data=sigma,
        type=pycdf.const.CDF_REAL4,
        compress=pycdf.const.NO_COMPRESSION,
        recVary=False)
    cdf['Sigma'].attrs['CATDESC'] = "Standard deviation of the Background calculated from 2016-04-09 to 20190623, multiplied by the Calibration Gain"
    cdf['Sigma'].attrs['DEPEND_1'] = 'Frequency'
    # to be changed?
    cdf['Sigma'].attrs['DICT_KEY'] = 'power>calibration'
    cdf['Sigma'].attrs['DISPLAY_TYPE'] = 'time_series'
    cdf['Sigma'].attrs['FIELDNAM'] = 'Sigma'
    cdf['Sigma'].attrs.new('FILLVAL', data=-1.0e+31, type=pycdf.const.CDF_REAL4)
    cdf['Sigma'].attrs['FORMAT'] = 'E12.4'
    cdf['Sigma'].attrs['LABLAXIS'] = 'Background Standard Deviation x Gain'
    cdf['Sigma'].attrs['UNITS'] = 'V**2 m**-2 Hz**-1'
    cdf['Sigma'].attrs.new('VALIDMIN', data=1e-18, type=pycdf.const.CDF_REAL4)
    cdf['Sigma'].attrs.new('VALIDMAX', data=1e-5, type=pycdf.const.CDF_REAL4)
    cdf['Sigma'].attrs['VAR_TYPE'] = "support_data"
    cdf['Sigma'].attrs['SCALETYP'] = "log"
    cdf['Sigma'].attrs.new('SCALEMIN', data=1e-18, type=pycdf.const.CDF_REAL4)
    cdf['Sigma'].attrs.new('SCALEMAX', data=1e-5, type=pycdf.const.CDF_REAL4)
    cdf['Sigma'].attrs['UCD'] = "phot.flux.density;obs.calib;stat.stdev"
    cdf.close()


    # Checking CDF
    check_cdf(filename=path_dir+"/"+cdfname+".cdf", verbose=verbose)
    print(f"The cdf file '{path_dir}/{cdfname}.cdf' has been produced and saved")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input_path', type=str, help="IDL savefile input path and filename", required=True)
    parser.add_argument('--o', dest='output_path', type=str, default="../../data_n2/spdyn_cdf_data_v02/", help="CDF savefile output path")
    parser.add_argument('-data_version', dest='data_version', type=int, default=2, help="Version of the data")
    parser.add_argument('-v', dest='verbose', help="Verbose flag", action='store_true')
    args = parser.parse_args()
    if os.path.exists(args.input_path):
        juno_lesia_calibrated_to_cdf(args.input_path, args.output_path, args.data_version, verbose=args.verbose)
    else:
        print(f"The file {args.input_path} doesn't exist")
