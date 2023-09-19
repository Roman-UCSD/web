###############################################################################
#
#  Pack ATLAS (BasicATLAS) and PHOENIX (Palantir) model atmospheres into FITS
#  files to be shared with the community
#
###############################################################################

import os, io
import re
import numpy as np
import f90nml
import pickle
import gzip, zipfile
from scipy import constants as scp
import periodictable
import datetime, time

script_dir = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.append(script_dir + '/../../CoolStar/HRFIT/')
sys.path.append(script_dir + '/../../CoolStar/BasicATLAS/')
import HR, atlas

from astropy.io import fits


chemical_symbols = {el.number: el.symbol for el in list(periodictable.elements)[1:]}
standard_abundances = np.loadtxt(script_dir + '/../../CoolStar/BasicATLAS/data/solar.csv', delimiter = ',', unpack = True, dtype = str, usecols = [0, 1, 2, 3, 4])


def get_phoenix_time(model):
    outputs = [model + '/output/structure/fort.6.gz']
    i = 1
    while os.path.isfile(backup := (model + '/output-{}.zip'.format(i))):
        outputs += [backup]
        i += 1

    phoenix_time = 0.0
    iter_count = 0
    start_dates = []
    for output in outputs:
        if output.split('.')[-1] == 'gz':
            f = gzip.open(model + '/output/structure/fort.6.gz', 'r')
            content = f.read().decode('utf-8')
            f.close()
        else:
            with zipfile.ZipFile(output) as z:
                found = False
                for file in z.namelist():
                    if file.split('/')[-1] == 'fort.6':
                        raise ValueError('Uncompressed fort.6 found')
                    if file.split('/')[-1] == 'fort.6.gz':
                        found = True
                        data = io.BytesIO(z.read(file))
                        with gzip.GzipFile(fileobj = data, mode = 'r') as f:
                            content = f.read().decode('utf-8')
                if not found:
                    raise ValueError('fort.6 not found!')

        start_dates += re.findall('This run was started at ([0-9]{4})([0-9]{2})([0-9]{2}):', content)
        start_dates[-1] = time.mktime(datetime.datetime(*np.array(start_dates[-1]).astype(int), 12, 00).timetuple())
        iter_count += len(re.findall('iter *[0-9]+: \(weighted\)', content))
        content = np.array(re.findall('time: used: *([0-9]+) sec', content)).astype(int)
        phoenix_time += (np.max(content) - np.min(content)) / 3600

    start_date = datetime.datetime.fromtimestamp(np.min(start_dates)).strftime('%Y-%m-%d')
    return phoenix_time, iter_count, start_date

def get_atlas_time(model):
    f = open(model + '/output.out', 'r')
    content = f.read()
    f.close()
    synthe_time = re.findall('Finished running SYNTHE in ([0-9]+):([0-9]+):([0-9.]+) s', content)[0]
    synthe_time = int(synthe_time[0]) + int(synthe_time[1]) / 60 + float(synthe_time[2]) / 3600

    meta = HR.atlas.meta(model)
    if model[-1] == '/':
        model = model[:-1]
    f = open(model + '/{}_{}_output.out'.format(grid := os.path.basename(model).split('_')[0], meta['logg']), 'r')
    content = f.read()
    f.close()
    atlas_models = re.findall('{}_{}/models/([0-9]{{3,5}}) successful'.format(grid, meta['logg']), content)
    index = np.where(np.array(atlas_models).astype(int) == meta['teff'])[0][0]
    atlas_time = re.findall('Finished running ATLAS-9 in ([0-9]+):([0-9]+):([0-9.]+) s', content)[index]
    atlas_time = int(atlas_time[0]) + int(atlas_time[1]) / 60 + float(atlas_time[2]) / 3600

    start = 0
    for i in range(index + 1):
        start = content.find('Launcher created', start + 1)
    content = content[start:content.find('successful', start)]
    niter = int(re.findall('([0-9]+) iterations completed', content)[-1])

    return atlas_time, synthe_time, niter

def get_phoenix_meta(model):
    """
    Determine all settings of a PHOENIX run
    """
    # Use the F90NML package to parse the PHOENIX namelist, which would be stored at the bottom of fort.20
    f = gzip.open(model + '/output/structure/fort.20.gz', 'r')
    content = f.read().decode('utf-8')
    f.close()
    content = content[content.find('&PHOENIX'):]
    f = open('temp.nml', 'w')
    f.write(content)
    f.close()
    content = f90nml.read('temp.nml')
    os.remove('temp.nml')
    content = content['PHOENIX']

    # Parse the pickled output
    f = open(model + '/output/structure/fort.pkl', 'rb')
    structure = pickle.load(f)
    f.close()

    output = {}
    output['software'] = 'PHOENIX ' + structure['param']['phxver'][0].decode('utf-8')
    output['teff'] = content['teff']
    output['logg'] = content['logg']
    yscale = float(content['yscale'])
    output['zscale'] = content['zscale']
    output['radius'] = content['r0'] / HR.apc.R_sun.cgs.value

    output['eheu'] = {}
    output['abs_eheu'] = {}
    total_mass = 0.0
    for i, nome in enumerate(content['nome']):
        if nome <= 0:
            continue
        nome = int(np.round(nome / 100))
        mn = float(standard_abundances[1][standard_abundances[0].astype(int) == nome][0].strip())
        abundance = float(standard_abundances[2][standard_abundances[0].astype(int) == nome][0].strip())
        lg_number = content['eheu'][i]
        if nome != 1 and nome != 2:
            lg_number += content['zscale']
        if nome == 2:
            lg_number += content['yscale']
        if chemical_symbols[nome] in ['O', 'Ne', 'Mg', 'Si', 'S', 'Ar', 'Ca', 'Ti']:
            lg_number += content['alpha_scale']
        lg_number = lg_number - np.array(content['eheu'])[np.round(np.array(content['nome']) / 100).astype(int) == 1][0] + 12.0
        number = 10 ** lg_number
        mass = number * mn
        if nome == 2:
            ymass = mass
        if nome == 1:
            xmass = mass
        total_mass += mass
        if nome > 2:
            output['eheu'][nome] = np.round(lg_number - output['zscale'] - abundance, 5)
        else:
            output['eheu'][nome] = np.round(lg_number - abundance, 5)
        output['abs_eheu'][nome] = np.round(lg_number, 5)
    output['X'] = xmass / total_mass
    output['Y'] = ymass / total_mass
    output['Z'] = 1 - xmass / total_mass - ymass / total_mass


    output['vturb'] = np.abs(content['xi'])
    output['mixlen'] = content['mixlng']



    if content['igrains'] == 0:
        output['dust'] = False
    elif content['igrains'] == 2:
        output['dust'] = True
    else:
        raise ValueError('Unknown value of igrains, {} in {}'.format(content['igrains'], model))
    output['clouds'] = content['use_clouds'] > 0
    output['cldalpha'] = content['cloud_covering_factor']
    output['settl'] = content['settleos']

    if content['model'] == 4:
        output['geometry'] = 'SPHERICAL'
    elif content['model'] == 5:
        output['geometry'] = 'PLANE-PARALLEL'
    else:
        raise ValueError('Unknown value of model, {} in {}'.format(content['model'], model))
    output['airorvac'] = 'VAC'

    return output

def reg_search(content, regex, fail_if_not_found = False, return_all = False):
    """
    Search for matches to a regular expression in a file

        content             :     File content or filename
        regex               :     Regular expression
        fail_if_not_found   :     Throw an error when no matches are found (defaults to False)
        return_all          :     Return all matches (otherwise, return the first match only) (defaults to False)

    returns:
        result              :     If "return_all", returns a list of all matched substrings. Otherwise, the first substring only
        content             :     Content of the file that was searched
    """
    if os.path.isfile(content):
        f = open(content, 'r')
        content = f.read()
        f.close()
    result = re.findall(regex, content)
    if len(result) == 0:
        if fail_if_not_found:
            raise ValueError('Broken file!')
        else:
            return False, content
    if not return_all:
        result = result[0]
    return result, content

def adapt_atlas(model, output, overwrite = False):
    if os.path.isfile(output) and (not overwrite):
        return

    # Load meta data
    meta = atlas.meta(model)
    eheu = meta['abun']; del meta['abun']
    settings = atlas.Settings()
    settings.abun = eheu
    settings.zscale = meta['zscale']
    settings.Y = meta['Y']
    meta['X'], meta['Y'], meta['Z'] = settings.mass_fractions()
    del meta['type']
    meta['software'] = 'ATLAS 9 / SYNTHE'
    del meta['res']
    meta['vturb'] = meta['synthe_vturb']
    del meta['synthe_vturb']
    meta['radius'] = 0.0
    result, content = reg_search(model + '/atlas_control.com', '\nCONVECTION *([^ ]+) *([^ ]+) *([^ ]+) *(.+)', return_all = True)
    meta['mixlen'] = float(result[0][1])
    meta['dust'] = False
    meta['clouds'] = False
    meta['settl'] = False
    meta['cldalpha'] = 0.0
    meta['geometry'] = 'PLANE-PARALLEL'
    result, content = reg_search(model + '/synthe_launch.com', '\n(AIR|VAC) +([^ ]+) +([^ ]+) +([^ ]+) +([^ ]+) +([^ ]+) +([^ ]+) +([^ ]+) +([^ ]+) +(.+)', return_all = True)
    meta['airorvac'] = result[0][0]
    meta['date'] = datetime.datetime.fromtimestamp(os.path.getmtime(model + '/synthe.tqs')).strftime('%Y-%m-%d')
    atlas_time, synthe_time, niter = get_atlas_time(model)
    meta['duration'] = atlas_time + synthe_time
    meta['niter'] = niter
    meta['bestiter'] = niter

    # Get absolute abundances
    abs_eheu = settings.abun_solar()
    eheu_column = []
    for elem in abs_eheu:
        if elem != 'H' and elem != 'He':
            abs_eheu[elem] += meta['zscale']
            if elem in eheu:
                abs_eheu[elem] += eheu[elem]
                eheu_column += [eheu[elem]]
            else:
                eheu_column += [00.00]
        if elem == 'He':
            nfracs = settings.atlas_abun()
            solar_he = abs_eheu[elem]
            abs_eheu[elem] = np.log10(nfracs[2] / nfracs[1]) + 12
            eheu_column += [abs_eheu[elem] - solar_he]
        if elem == 'H':
            eheu_column += [00.00]

    structure, units = atlas.read_structure(model)
    structure['radiative_flux'] = structure['radiative_flux'] * 4 * np.pi      # Convert fluxes from Eddington values to conventional
    structure['convective_flux'] = structure['convective_flux'] * 4 * np.pi
    spectrum = atlas.read_spectrum(model)

    # Write out into a FITS file
    tables = []

    columns = {
        'wl': fits.Column(name = 'Wavelength', array = spectrum['wl'], format = 'D', unit = 'AA'),
        'flux': fits.Column(name = 'Total flux density', array = spectrum['flux'] * np.pi, format = 'D', unit = 'ERG/S/CM2/AA'),
        'cont': fits.Column(name = 'Continuum flux density', array = spectrum['cont'] * np.pi, format = 'D', unit = 'ERG/S/CM2/AA'),
    }
    hdr = fits.Header()
    hdr['TABLE'] = 'Emergent spectrum'
    tables += [fits.BinTableHDU.from_columns(columns.values(), header = hdr)]

    replacements = {
        'Mass column density' : 'Column mass density',
        'Density' : 'Mass density',
        'Flux error derivative' : 'Derivative flux error',
        'Rosseland opacity' : 'Rosseland mass opacity',
    }
    columns = {}
    for column in structure:
        name = column[0].upper() + column[1:].replace('_', ' ')
        if name in replacements:
            name = replacements[name]
        if name in ['Layer']:
            continue
        unit = units[column].upper().replace('^', '')
        unit = re.sub('([A-Z]+)-([2-9])', '/\\1\\2', unit)
        unit = re.sub('([A-Z]+)-(1)', '/\\1', unit)
        unit = unit.replace(' ', '')
        if unit == '':
            unit = 'NONE'
        if unit[0] == '/':
            unit = '1' + unit
        if unit == 'BA':
            unit = 'DYN/CM2'
        if unit == 'PERCENT':
            unit = '%'
        columns[column] = fits.Column(name = name, array = structure[column], format = 'D', unit = unit)
    hdr = fits.Header()
    hdr['TABLE'] = 'Atmospheric structure'
    tables += [fits.BinTableHDU.from_columns(columns.values(), header = hdr)]


    columns = {
        'element': fits.Column(name = 'Element', array = list(abs_eheu.keys()), format = 'A2'),
        'eheu': fits.Column(name = 'Relative abundance', array = eheu_column, format = 'D', unit = 'DEX'),
        'abs_eheu': fits.Column(name = 'Absolute abundance', array = list(abs_eheu.values()), format = 'D', unit = 'DEX'),
    }
    hdr = fits.Header()
    hdr['TABLE'] = 'Chemical composition'
    tables += [fits.BinTableHDU.from_columns(columns.values(), header = hdr)]

    hdr = fits.Header()
    for key in meta:
        hdr[key] = meta[key]

    hdul = fits.HDUList([fits.PrimaryHDU(header=hdr)] + tables)
    hdul.writeto(output, overwrite = overwrite)

def adapt_phoenix(model, output, overwrite = False):
    if os.path.isfile(output) and (not overwrite):
        return

    # Load meta data
    meta = get_phoenix_meta(model)
    phoenix_time, iter_count, start_date = get_phoenix_time(model)
    meta['date'] = start_date
    meta['duration'] = phoenix_time

    # Load model for the iteration with the best average flux error
    min_err, max_err, avg_err = HR.phoenix_convergence(model)
    iteration = np.argmin(avg_err) + 1
    structure = HR._phoenix_explorer(model, 'structure', iteration)
    f_structure = HR._phoenix_explorer(model, 'further_structure', iteration)
    spectrum = HR._phoenix_explorer(model, 'spectrum', iteration)
    if len(min_err) != iter_count:
        raise ValueError('Number of iterations mismatch in model {}'.format(model))
    meta['niter'] = len(min_err)
    meta['bestiter'] = iteration

    # Write out into a FITS file
    tables = []

    columns = {
        'wl': fits.Column(name = 'Wavelength', array = spectrum[0], format = 'D', unit = 'AA'),
        'flux': fits.Column(name = 'Flux density', array = spectrum[1] * np.pi, format = 'D', unit = 'ERG/S/CM2/AA'),
    }
    hdr = fits.Header()
    hdr['TABLE'] = 'Emergent spectrum'
    tables += [fits.BinTableHDU.from_columns(columns.values(), header = hdr)]

    columns = {
        'tau': fits.Column(name = 'Rosseland optical depth', array = structure[0], format = 'D', unit = 'NONE'),
        'temp': fits.Column(name = 'Temperature', array = structure[1], format = 'D', unit = 'K'),
        'pgas': fits.Column(name = 'Gas pressure', array = structure[2], format = 'D', unit = 'DYN/CM2'),
        'kappa': fits.Column(name = 'Rosseland mass opacity', array = structure[3], format = 'D', unit = 'CM2/G'),
        'rho': fits.Column(name = 'Mass density', array = structure[4], format = 'D', unit = 'G/CM3'),
        'cmass': fits.Column(name = 'Column mass density', array = f_structure[4], format = 'D', unit = 'G/CM2'),
        'mu': fits.Column(name = 'Mean molecular weight', array = f_structure[6], format = 'D', unit = 'U'),
        'fconv': fits.Column(name = 'Convective flux fraction', array = f_structure[5] / (scp.sigma * 1e3 * meta['teff'] ** 4 / (4 * np.pi)), format = 'D', unit = 'NONE'),
    }
    hdr = fits.Header()
    hdr['TABLE'] = 'Atmospheric structure'
    tables += [fits.BinTableHDU.from_columns(columns.values(), header = hdr)]

    columns = {
        'iter': fits.Column(name = 'Temperature iteration number', array = np.arange(1, len(min_err) + 1), format = 'I'),
        'min': fits.Column(name = 'Minimum error', array = min_err, format = 'D', unit = '%'),
        'avg': fits.Column(name = 'Average error', array = avg_err, format = 'D', unit = '%'),
        'max': fits.Column(name = 'Maximum error', array = max_err, format = 'D', unit = '%'),
    }
    hdr = fits.Header()
    hdr['TABLE'] = 'Radiative flux errors'
    tables += [fits.BinTableHDU.from_columns(columns.values(), header = hdr)]

    eheu = np.array([[chemical_symbols[z], meta['eheu'][z], meta['abs_eheu'][z]] for z in sorted(list(meta['eheu'].keys()))]).T
    columns = {
        'element': fits.Column(name = 'Element', array = eheu[0], format = 'A2'),
        'eheu': fits.Column(name = 'Relative abundance', array = eheu[1].astype(float), format = 'D', unit = 'DEX'),
        'abs_eheu': fits.Column(name = 'Absolute abundance', array = eheu[2].astype(float), format = 'D', unit = 'DEX'),
    }
    del meta['eheu'], meta['abs_eheu']
    hdr = fits.Header()
    hdr['TABLE'] = 'Chemical composition'
    tables += [fits.BinTableHDU.from_columns(columns.values(), header = hdr)]

    hdr = fits.Header()
    for key in meta:
        hdr[key] = meta[key]

    hdul = fits.HDUList([fits.PrimaryHDU(header=hdr)] + tables)
    hdul.writeto(output, overwrite = overwrite)
