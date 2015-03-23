#from Helper.Initializer import initializer

class Diamond(object):

    def __init__(self,name = None):

        self.Specifications = {
            'Name': '',
            'Manufacturer':'',
            'CrystalStructure':'', # eg. 'poly' or 'single'
            'Irradiation':'', # eg. 'no' or 'neutron' or 'proton'
            'Thickness':0, # in um
            'DetectorType':'Pad', # eg. 'Pad' or 'Pixel'
        }

        self.Position = {
                    'xmin': -0.2, # standard position
                    'xmax': 0.15,
                    'ymin': 0.03,
                    'ymax': 0.40
                }

        if name is not None:
            self.SetInfoFromName(name)


    def SetInfoFromName(self,diamondname):
        diamonds = {
            'S30':{
                'crystal':      'single',
                'irradiated':   'neutron',
                'thickness':    539,
                'producer':     'E-6',
                'position': {
                    'xmin': -0.19,
                    'xmax': 0.16,
                    'ymin': 0.01,
                    'ymax': 0.38
                }
            },
            'S125':{
                'crystal':      'single',
                'irradiated':   'proton',
                'thickness':    '',
                'producer':     'E-6',
                'position': {
                    'xmin': -0.23,
                    'xmax': 0.12,
                    'ymin': -0.05,
                    'ymax': 0.32
                }
            },
            'S129':{
                'crystal':      'single',
                'irradiated':   'no',
                'thickness':    528,
                'producer':     'E-6',
                'position': {
                    'xmin': -0.17,
                    'xmax': 0.18,
                    'ymin': -0.14,
                    'ymax': 0.23
                }
            },
            'S129-old-box':{
                'crystal':      'single',
                'irradiated':   'no',
                'thickness':    '',
                'producer':     'II-IV',
                'position': {
                    'xmin': -0.18,
                    'xmax': 0.17,
                    'ymin': -0.09,
                    'ymax': 0.28
                }
            },
            'IIa-1':{
                'crystal':      'single',
                'irradiated':   'no',
                'thickness':    617,
                'producer':     'II-a',
                'position': {
                    'xmin': -0.18,
                    'xmax': 0.17,
                    'ymin': 0.01,
                    'ymax': 0.38
                }
            },
            'IIa-2':{
                'crystal':      'single',
                'irradiated':   'neutron',
                'thickness':    538,
                'producer':     'II-a',
                'position': {
                    'xmin': -0.2,
                    'xmax': 0.15,
                    'ymin': 0.02,
                    'ymax': 0.39
                }
            },
            'IIa-3':{
                'crystal':      'single',
                'irradiated':   'neutron',
                'thickness':    535,
                'producer':     'II-a',
                'position': {
                    'xmin': -0.21,
                    'xmax': 0.14,
                    'ymin': -0.04,
                    'ymax': 0.33
                }
            },
            'IIa-5':{
                'crystal':      'single',
                'irradiated':   'no',
                'thickness':    775,
                'producer':     'II-a',
                'position': {
                    'xmin': -0.2,
                    'xmax': 0.15,
                    'ymin': -0.09,
                    'ymax': 0.28
                }
            },
            '2A87-E':{
                'crystal':      'poly',
                'irradiated':   'neutron',
                'thickness':    500,
                'producer':     'II-IV',
                'position': {
                    'xmin': -0.22,
                    'xmax': 0.13,
                    'ymin': -0.05,
                    'ymax': 0.32
                }
            },
        }
        self.Specifications['Name'] = diamondname
        self.Specifications['Manufacturer'] = diamonds[diamondname]['producer']
        self.Specifications['CrystalStructure'] = diamonds[diamondname]['crystal']
        self.Specifications['Irradiation'] = diamonds[diamondname]['irradiated']
        self.Specifications['Thickness'] = diamonds[diamondname]['thickness']
        self.Position['xmin'] = diamonds[diamondname]['position']['xmin']
        self.Position['xmax'] = diamonds[diamondname]['position']['xmax']
        self.Position['ymin'] = diamonds[diamondname]['position']['ymin']
        self.Position['ymax'] = diamonds[diamondname]['position']['ymax']


    def SetDetectorType(self,type = 'Pad'):
        self.Specifications['DetectorType'] = type