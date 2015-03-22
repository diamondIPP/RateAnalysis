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

        self.Attributes = {

        }

        if name is not None:
            self.SetInfoFromName(name)


    def SetInfoFromName(self,diamondname):
        diamonds = {
            'S30':{
                'crystal':      'single',
                'irradiated':   'neutron',
                'thickness':    539,
                'producer':     'E-6'
            },
            'S125':{
                'crystal':      'single',
                'irradiated':   'proton',
                'thickness':    '',
                'producer':     'E-6'
            },
            'S129':{
                'crystal':      'single',
                'irradiated':   'no',
                'thickness':    528,
                'producer':     'E-6'
            },
            'S129-old-box':{
                'crystal':      'single?',
                'irradiated':   '?',
                'thickness':    '',
                'producer':     'II-IV'
            },
            'IIa-1':{
                'crystal':      'single',
                'irradiated':   'no',
                'thickness':    617,
                'producer':     'II-a'
            },
            'IIa-2':{
                'crystal':      'single',
                'irradiated':   'neutron',
                'thickness':    538,
                'producer':     'II-a'
            },
            'IIa-3':{
                'crystal':      'single',
                'irradiated':   'neutron',
                'thickness':    535,
                'producer':     'II-a'
            },
            'IIa-4':{
                'crystal':      'single',
                'irradiated':   'no',
                'thickness':    775,
                'producer':     'II-a'
            },
            'IIa-5':{
                'crystal':      'single?',
                'irradiated':   '?',
                'thickness':    '',
                'producer':     'II-a'
            },
            '2A87-E':{
                'crystal':      'poly',
                'irradiated':   'neutron',
                'thickness':    500,
                'producer':     'II-IV'
            },
        }
        self.Specifications['Name'] = diamondname
        self.Specifications['Manufacturer'] = diamonds[diamondname]['producer']
        self.Specifications['CrystalStructure'] = diamonds[diamondname]['crystal']
        self.Specifications['Irradiation'] = diamonds[diamondname]['irradiated']
        self.Specifications['Thickness'] = diamonds[diamondname]['thickness']

    def SetDetectorType(self,type = 'Pad'):
        self.Specifications['DetectorType'] = type